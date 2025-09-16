#!/usr/bin/env python3
import argparse
import os
import time
from Bio import SeqIO, Seq
import csv
import sys
import logging
from collections import Counter, namedtuple
import sys

def rc(sequence):
    sequence = Seq.Seq(sequence)
    return str(sequence.reverse_complement())


# Sliding window utilities
# A tuple with the elements for generating a sequence start/end paired with a percentile and its count
WindowTask = namedtuple("WindowTask", ["seq", "queries", "start", "end", "headers"])
def window_tasks(sequence_record, queries, width, spacing, headers, offset=0):
    seq_len = len(sequence_record)

    for i in range(0, seq_len, spacing):
        window = sequence_record[i:i+width]
        if not window:
            continue
        yield WindowTask(window, queries, i+offset, min(i+width+offset, seq_len+offset), headers)

def windowCount(seq, mer_size):
    windowSeq, queries, start, end, headers = seq
    windowSeq = str(windowSeq)
    counts = Counter(windowSeq[i:i + mer_size] for i in range(len(windowSeq) - 2))
    return headers + [start + 1, min(start + len(windowSeq), end), sum([counts.get(q, 0) for q in queries])]

def inverseWindowCount(seq, mer_size):
    windowSeq, queries, start, end, headers = seq
    windowSeq = str(windowSeq)
    counts = Counter(windowSeq[i:i + mer_size] for i in range(len(windowSeq) - 2))

    return headers + [start + 1, min(start + len(windowSeq), end), sum([counts.get(count, 0) for count in counts if count in queries])]

# CLI arguments
parser = argparse.ArgumentParser(description="Kmer Mass Query: localize percentile kmers in genomic sequences (multi-pass)")
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
parser.add_argument('-x', '--inverse', action='store_true',
                    help='Tabulate the inverse (for large kmer lists with-respect-to window sizes)')

args = parser.parse_args()

# Logging
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()
if args.verbose:
    logger.setLevel(logging.INFO)

startTime = time.time()

if not os.path.exists(args.input):
    logger.error(f"Couldn't find input sequence file '{args.input}'")
    exit()
elif not os.path.exists(args.query):
    logger.error(f"Couldn't find input query file '{args.query}'")
    exit()

# If input chunksize and window size are not compatible, lower chunksize to the next compatible length
if args.chunk_size%args.width!=0:
    args.chunksize-=args.chunk_size%args.width

callableProcess = inverseWindowCount if args.inverse else windowCount

# Count total kmers
logger.info("Counting kmers in query file...")
tableLength = sum(1 for _ in open(args.query)) - 1  # minus header
logger.info(f"Query has {tableLength:,} kmers")
# Multi-pass over percentile bins
step = args.percentile_step
bin_edges = list(range(0, 100, step))

# Infer k from first kmer in query file
with open(args.query) as f:
    f.readline()  # header
    first_line = f.readline().strip().split("\t")
    k = len(first_line[0])

# Prepare output
with open(args.output, "w", newline="") as fileOutput:
    tsvWriter = csv.writer(fileOutput, delimiter="\t")
    header = ["Sequence", "Bin", "Start", "End", "Count"]
    tsvWriter.writerow(header)
    # Parse genome once per bin
    for b in bin_edges:
        if args.percentile_keep-1 < b < 100-args.percentile_keep:
            pass
        else:
            start_idx = int(b / 100 * tableLength)
            end_idx = int(min(100, b + step) / 100 * tableLength)
            bin_name = f"pct{b + step:03d}"
            logger.info(f"Loading kmers for bin {bin_name} (lines {start_idx}-{end_idx})")
            kmers = set()
            with open(args.query) as f:
                f.readline()  # skip header
                for idx, line in enumerate(f):
                    if idx < start_idx:
                        continue
                    if idx >= end_idx:
                        break
                    kmer = line.strip().split("\t")[0]
                    kmers.add(kmer)
                    if args.complement:
                        kmers.add(rc(kmer))
            logger.info(f"Bin {bin_name} has {len(kmers):,} kmers")
            # Scan genome for this bin
            for record in SeqIO.parse(args.input, args.format):
                sequence_name = record.id
                seq_record = str(record.seq).upper()
                headers = [sequence_name, bin_name]
                sequenceLength = len(seq_record)
                if sequenceLength > args.chunk_size:
                    logging.info(f"Sequence {sequence_name} is large. Breaking into smaller segments")
                    indices = range(0, sequenceLength, args.chunk_size)
                    for sequenceIndex in indices:
                        sub_seq = seq_record[sequenceIndex:sequenceIndex + args.chunk_size]
                        for task in window_tasks(sub_seq, kmers, args.width, args.spacing, headers, offset=sequenceIndex):
                            row = callableProcess(task, args.mer_size)
                            tsvWriter.writerow(row)
                        del sub_seq
                else:
                    for task in window_tasks(seq_record, kmers, args.width, args.spacing, headers, offset=0):
                        row = callableProcess(task, args.mer_size)
                        tsvWriter.writerow(row)
                logging.info(f"{bin_name} kmers on sequence {sequence_name} windows written to output file")

logger.info(f"Execution time in seconds: {time.time() - startTime:.2f}")