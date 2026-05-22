#!/usr/bin/env python3
import argparse
import os
import time
import math
from Bio import SeqIO
import csv
import logging
from collections import Counter
import multiprocessing

# Fast reverse complement using translation table
RC_TRANS = str.maketrans("ACGTacgt", "TGCAtgca")
def rc(sequence):
    return sequence.translate(RC_TRANS)[::-1]

# Global variables for workers (shared via Copy-on-Write)
GLOBAL_KMER_MAP = None
GLOBAL_MER_SIZE = None
GLOBAL_BIN_NAMES = None

# Sliding window utility
def window_tasks(sequence_record, width, spacing, sequence_name, offset=0):
    seq_len = len(sequence_record)
    for i in range(0, seq_len, spacing):
        window = sequence_record[i:i+width]
        if not window:
            continue
        start = i + offset
        window_str = str(window)
        end = start + len(window_str)
        yield (window_str, sequence_name, start, end)

def process_window(task):
    window_seq, sequence_name, start, end = task
    local_counts = Counter()
    for i in range(len(window_seq) - GLOBAL_MER_SIZE + 1):
        kmer = window_seq[i:i + GLOBAL_MER_SIZE]
        bin_ids = GLOBAL_KMER_MAP.get(kmer)
        if bin_ids is not None:
            if isinstance(bin_ids, int):
                local_counts[bin_ids] += 1
            else:
                for bin_id in bin_ids:
                    local_counts[bin_id] += 1

    rows = []
    for bin_id, bin_name in enumerate(GLOBAL_BIN_NAMES):
        count = local_counts.get(bin_id, 0)
        rows.append([sequence_name, bin_name, start + 1, end, count])
    return rows

def init_worker(kmer_map, mer_size, bin_names):
    global GLOBAL_KMER_MAP
    global GLOBAL_MER_SIZE
    global GLOBAL_BIN_NAMES
    GLOBAL_KMER_MAP = kmer_map
    GLOBAL_MER_SIZE = mer_size
    GLOBAL_BIN_NAMES = bin_names

def main():
    # CLI arguments
    parser = argparse.ArgumentParser(description="Kmer Mass Query: localize percentile kmers in genomic sequences")
    parser.add_argument('-i', '--input', dest='input', required=True, help='Input sequence file (FASTA/FASTQ)')
    parser.add_argument('-f', '--format', default='fasta', help='Input file type [default fasta]')
    parser.add_argument('-q', '--query', required=True, help='Ranked query table file (tsv)')
    parser.add_argument('-w', '--width', type=int, default=3000, help='Window width [default 3000]')
    parser.add_argument('-s', '--spacing', type=int, default=3000, help='Window spacing [default 3000]')
    parser.add_argument('-o', '--output', default='mass_query_report.tsv', help='Output TSV file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose mode')
    parser.add_argument('-c', '--complement', action='store_true', help='Include reverse complements in kmer sets [default False]')
    parser.add_argument('-m', '--mer-size', dest='mer_size', type=int, help='kmer size in query [default 20]', default=20)
    parser.add_argument('--percentile-low', dest='percentile_low', type=float, default=5, help='Bottom N percent of kmers to keep [default 5]')
    parser.add_argument('--percentile-high', dest='percentile_high', type=float, default=5, help='Top N percent of kmers to keep [default 5]')
    parser.add_argument('-e', '--percentile-keep', type=int, dest='percentile_keep', default=None, help='Deprecated: use --percentile-low and --percentile-high instead')
    parser.add_argument('-k', '--chunk-size', dest='chunk_size', type=int, help='Max chunk size to work on [default 30000000]', default=30000000)
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads for parallel processing [default 1]')
    parser.add_argument('--minimum-size', dest='minimum_size', type=int, help='Minimum sequence size to include.', default=15000)

    args = parser.parse_args()

    # Logging
    logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')
    logger = logging.getLogger()
    if args.verbose:
        logger.setLevel(logging.INFO)

    start_time = time.time()

    if not os.path.exists(args.input):
        logger.error(f"Couldn't find input sequence file '{args.input}'")
        return
    elif not os.path.exists(args.query):
        logger.error(f"Couldn't find input query file '{args.query}'")
        return

    # If input chunksize and window size are not compatible, lower chunksize to the next compatible length
    if args.chunk_size % args.width != 0:
        args.chunk_size -= args.chunk_size % args.width

    if args.percentile_keep is not None:
        args.percentile_low = args.percentile_keep
        args.percentile_high = args.percentile_keep

    # Count total kmers
    logger.info("Processing query file...")
    try:
        with open(args.query) as f:
            f.readline() # skip header
            table_length = sum(1 for _ in f)
    except Exception as e:
        logger.error(f"Error reading query file: {e}")
        return

    logger.info(f"Query has {table_length:,} kmers")

    # Define thresholds
    low_cutoff = int(args.percentile_low / 100 * table_length)
    high_cutoff = int((100 - args.percentile_high) / 100 * table_length)

    # Label "low" and "high"
    bin_names = ["low", "high"]
    # 0 = low, 1 = high

    # Map kmers to "low" and "high" labels in a single pass
    logger.info("Loading kmers into sets...")
    kmer_map = {}

    def add_to_map(k, b_id):
        existing = kmer_map.get(k)
        if existing is None:
            kmer_map[k] = b_id
        elif isinstance(existing, int):
            if existing != b_id:
                kmer_map[k] = [existing, b_id]
        else:
            if b_id not in existing:
                existing.append(b_id)

    try:
        with open(args.query) as f:
            f.readline() # skip header
            for idx, line in enumerate(f):
                if idx < low_cutoff:
                    kmer = line.strip().split("\t")[0].upper()
                    add_to_map(kmer, 0)
                    if args.complement:
                        add_to_map(rc(kmer), 0)
                elif idx >= high_cutoff:
                    kmer = line.strip().split("\t")[0].upper()
                    add_to_map(kmer, 1)
                    if args.complement:
                        add_to_map(rc(kmer), 1)
    except Exception as e:
        logger.error(f"Error mapping kmers: {e}")
        return

    logger.info(f"Loaded {len(kmer_map):,} unique kmers across {len(bin_names)} sets")

    # Prepare for parallel processing using global variables for Copy-on-Write sharing
    global GLOBAL_KMER_MAP
    global GLOBAL_MER_SIZE
    global GLOBAL_BIN_NAMES
    GLOBAL_KMER_MAP = kmer_map
    GLOBAL_MER_SIZE = args.mer_size
    GLOBAL_BIN_NAMES = bin_names

    # Calculate genome assembly quality metric
    # Metric: accumulation of extreme k-mers relative to random k-mers
    # We'll compute this by sequence later if needed, but here we can prepare global stats
    logger.info("Computing assembly quality metric components...")

    # Use 'fork' to share memory efficiently on Unix-like systems
    try:
        mp_context = multiprocessing.get_context('fork')
    except ValueError:
        # 'fork' not available (e.g. on Windows), fallback to default
        mp_context = multiprocessing.get_context()
        logger.warning("Multiprocessing 'fork' not available. Memory usage may be higher.")

    if mp_context.get_start_method() == 'fork':
        pool = mp_context.Pool(processes=args.threads)
    else:
        pool = mp_context.Pool(processes=args.threads, initializer=init_worker, initargs=(kmer_map, args.mer_size, bin_names))

    # Prepare output
    try:
        quality_metrics = []
        with open(args.output, "w", newline="") as file_output:
            tsv_writer = csv.writer(file_output, delimiter="\t")
            tsv_writer.writerow(["Sequence", "Bin", "Start", "End", "Count"])

            # Scan genome once
            for record in SeqIO.parse(args.input, args.format):
                sequence_name = record.id
                sequence_length = len(record)
                if sequence_length < args.minimum_size:
                    continue

                logger.info(f"Processing sequence {sequence_name} ({sequence_length:,} bp)")
                seq_total_extreme_low = 0
                seq_total_extreme_high = 0
                seq_windows = 0

                if sequence_length > args.chunk_size:
                    for i in range(0, sequence_length, args.chunk_size):
                        sub_seq = str(record.seq[i:i + args.chunk_size]).upper()
                        tasks = window_tasks(sub_seq, args.width, args.spacing, sequence_name, offset=i)
                        for result_rows in pool.imap(process_window, tasks):
                            seq_windows += 1
                            for row in result_rows:
                                tsv_writer.writerow(row)
                                bin_name = row[1]
                                count = row[4]
                                # Labels are "low" or "high"
                                if bin_name == "low":
                                    seq_total_extreme_low += count
                                else:
                                    seq_total_extreme_high += count
                        del sub_seq
                else:
                    seq_str = str(record.seq).upper()
                    tasks = window_tasks(seq_str, args.width, args.spacing, sequence_name)
                    for result_rows in pool.imap(process_window, tasks):
                        seq_windows += 1
                        for row in result_rows:
                            tsv_writer.writerow(row)
                            bin_name = row[1]
                            count = row[4]
                            if bin_name == "low":
                                seq_total_extreme_low += count
                            else:
                                seq_total_extreme_high += count
                    del seq_str

                # Calculate sequence-level metric
                # Extreme k-mer density (extreme k-mers per bp)
                if sequence_length > 0:
                    low_density = seq_total_extreme_low / sequence_length
                    high_density = seq_total_extreme_high / sequence_length
                    quality_metrics.append({
                        'Sequence': sequence_name,
                        'Length': sequence_length,
                        'ExtremeLowCount': seq_total_extreme_low,
                        'ExtremeHighCount': seq_total_extreme_high,
                        'LowDensity': low_density,
                        'HighDensity': high_density
                    })

        # Write quality metrics to a separate file
        metrics_output = args.output.replace('.tsv', '_metrics.tsv')
        with open(metrics_output, 'w', newline='') as f:
            # Include metadata about thresholds used
            f.write(f"# Low_Percentile_Threshold: {args.percentile_low}\n")
            f.write(f"# High_Percentile_Threshold: {args.percentile_high}\n")
            fieldnames = ['Sequence', 'Length', 'ExtremeLowCount', 'ExtremeHighCount', 'LowDensity', 'HighDensity']
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for m in quality_metrics:
                writer.writerow(m)
        logger.info(f"Assembly quality metrics written to {metrics_output}")
    except Exception as e:
        logger.error(f"Error during mass query processing: {e}")
    finally:
        pool.close()
        pool.join()

    logger.info(f"Execution time in seconds: {time.time() - start_time:.2f}")

if __name__ == "__main__":
    main()
