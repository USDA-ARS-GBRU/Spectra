#!/usr/bin/env python3
import argparse
import os
import time
from Bio import SeqIO
import csv
import logging
from collections import Counter, defaultdict
import multiprocessing

# Fast reverse complement using translation table
RC_TRANS = str.maketrans("ACGTacgt", "TGCAtgca")
def rc(sequence):
    return sequence.translate(RC_TRANS)[::-1]

# Sliding window utility
def window_tasks(sequence_record, width, spacing, sequence_name, offset=0):
    seq_len = len(sequence_record)
    for i in range(0, seq_len, spacing):
        window = sequence_record[i:i+width]
        if not window:
            continue
        start = i + offset
        # The original script uses min(i+width+offset, seq_len+offset) for end
        # and start + len(windowSeq) for end in windowCount.
        # We'll use the actual window length to be precise.
        window_str = str(window)
        end = start + len(window_str)
        yield (window_str, sequence_name, start, end)

def process_window(task):
    windowSeq, sequence_name, start, end = task
    local_counts = Counter()
    for i in range(len(windowSeq) - GLOBAL_MER_SIZE + 1):
        kmer = windowSeq[i:i + GLOBAL_MER_SIZE]
        if kmer in GLOBAL_KMER_MAP:
            for bin_name in GLOBAL_KMER_MAP[kmer]:
                local_counts[bin_name] += 1

    rows = []
    # To maintain compatibility with multi-pass output, we return rows for all bins
    for bin_name in GLOBAL_BINS:
        rows.append([sequence_name, bin_name, start + 1, end, local_counts.get(bin_name, 0)])
    return rows

def init_worker(kmer_map, mer_size, bins):
    global GLOBAL_KMER_MAP
    global GLOBAL_MER_SIZE
    global GLOBAL_BINS
    GLOBAL_KMER_MAP = kmer_map
    GLOBAL_MER_SIZE = mer_size
    GLOBAL_BINS = bins

def main():
    # CLI arguments
    parser = argparse.ArgumentParser(description="Kmer Mass Query: localize percentile kmers in genomic sequences (optimized)")
    parser.add_argument('-i', '--input', dest='input', required=True, help='Input sequence file (FASTA/FASTQ)')
    parser.add_argument('-f', '--format', default='fasta', help='Input file type [default fasta]')
    parser.add_argument('-q', '--query', required=True, help='Ranked query table file (tsv)')
    parser.add_argument('-w', '--width', type=int, default=3000, help='Window width [default 3000]')
    parser.add_argument('-s', '--spacing', type=int, default=3000, help='Window spacing [default 3000]')
    parser.add_argument('-o', '--output', default='mass_query_report.tsv', help='Output TSV file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose mode')
    parser.add_argument('-c', '--complement', action='store_true', help='Include reverse complements in kmer bins [default False]')
    parser.add_argument('-m', '--mer-size', dest='mer_size', type=int, help='kmer size in query [default 20]', default=20)
    parser.add_argument('-p', '--percentile-step', type=int, default=1, help='Step size for percentile bins [default 1]')
    parser.add_argument('-e', '--percentile-keep', type=int, dest='percentile_keep', default=5, help='Extreme kmer tabulation. Top and bottom N percent kept [default 5]')
    parser.add_argument('-k', '--chunk-size', dest='chunk_size', type=int, help='Max chunk size to work on [default 30000000]', default=30000000)
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads for parallel processing [default 1]')
    parser.add_argument('-x', '--inverse', action='store_true', help='(Deprecated) Tabulate the inverse (now always efficient)')
    parser.add_argument('--minimum-size', dest='minimum_size', type=int, help='Minimum sequence size to include.', default=15000)

    args = parser.parse_args()

    # Logging
    logging.basicConfig(level=logging.INFO if args.verbose else logging.ERROR, format='%(levelname)s: %(message)s')
    logger = logging.getLogger()

    startTime = time.time()

    if not os.path.exists(args.input):
        logger.error(f"Couldn't find input sequence file '{args.input}'")
        return
    elif not os.path.exists(args.query):
        logger.error(f"Couldn't find input query file '{args.query}'")
        return

    # If input chunksize and window size are not compatible, lower chunksize to the next compatible length
    if args.chunk_size % args.width != 0:
        args.chunk_size -= args.chunk_size % args.width

    # Count total kmers
    logger.info("Counting kmers in query file...")
    with open(args.query) as f:
        f.readline() # skip header
        tableLength = sum(1 for _ in f)

    logger.info(f"Query has {tableLength:,} kmers")

    # Identify bins of interest
    step = args.percentile_step
    bin_edges = list(range(0, 100, step))
    interest_bins = []
    for b in bin_edges:
        if not (args.percentile_keep - 1 < b < 100 - args.percentile_keep):
            interest_bins.append(b)

    # Map kmers to bins
    logger.info("Loading kmers into bins...")
    kmer_map = defaultdict(list)
    bin_names = []
    for b in interest_bins:
        start_idx = int(b / 100 * tableLength)
        end_idx = int(min(100, b + step) / 100 * tableLength)
        bin_name = f"pct{b + step:03d}"
        bin_names.append(bin_name)

        with open(args.query) as f:
            f.readline() # skip header
            for idx, line in enumerate(f):
                if idx < start_idx: continue
                if idx >= end_idx: break
                kmer = line.strip().split("\t")[0].upper()
                kmer_map[kmer].append(bin_name)
                if args.complement:
                    kmer_map[rc(kmer)].append(bin_name)

    logger.info(f"Loaded {len(kmer_map):,} unique kmers across {len(bin_names)} bins")

    # Prepare for parallel processing
    pool = multiprocessing.Pool(processes=args.threads, initializer=init_worker, initargs=(kmer_map, args.mer_size, bin_names))

    # Prepare output
    with open(args.output, "w", newline="") as fileOutput:
        tsvWriter = csv.writer(fileOutput, delimiter="\t")
        tsvWriter.writerow(["Sequence", "Bin", "Start", "End", "Count"])

        # Scan genome once
        for record in SeqIO.parse(args.input, args.format):
            sequence_name = record.id
            sequenceLength = len(record)
            if sequenceLength < args.minimum_size:
                continue

            logger.info(f"Processing sequence {sequence_name} ({sequenceLength:,} bp)")

            if sequenceLength > args.chunk_size:
                for i in range(0, sequenceLength, args.chunk_size):
                    sub_seq = str(record.seq[i:i + args.chunk_size]).upper()
                    tasks = window_tasks(sub_seq, args.width, args.spacing, sequence_name, offset=i)
                    for result_rows in pool.imap(process_window, tasks):
                        for row in result_rows:
                            tsvWriter.writerow(row)
                    del sub_seq
            else:
                seq_str = str(record.seq).upper()
                tasks = window_tasks(seq_str, args.width, args.spacing, sequence_name)
                for result_rows in pool.imap(process_window, tasks):
                    for row in result_rows:
                        tsvWriter.writerow(row)
                del seq_str

    pool.close()
    pool.join()
    logger.info(f"Execution time in seconds: {time.time() - startTime:.2f}")

if __name__ == "__main__":
    main()