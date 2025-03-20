#!/usr/bin/env python3

# TODO: create scalability for n-mer

import os
import time
from Bio import SeqIO
import csv
import re
import multiprocessing as mp
import logging

logging.basicConfig(level=logging.INFO)

# current, most-efficient window counter
def mpWindowCount(seq, queries, start, end, names):
    return names + [start + 1, end] + [seq.count_overlap(a) for a in queries]

# broken function - returns non-overlapping chunks?
def pdWindowCount(seq, queries, start, end, names):
    return names + [start + 1, end] + [len(re.findall(f'(?={a})', str(seq))) for a in queries]

# original, slow window counter
def windowCount(seq, queries, proportions):
    counts = [0] * len(queries)
    for i in range(0, len(seq) - 2):
        if seq[i:i + 3] in queries:
            counts[queries.index(seq[i:i + 3])] += 1
    if proportions:
        total = sum(counts)
        if total:
            return [a / total for a in counts]
        else:
            return [0] * 64
    else:
        return counts

def execute(args):
    startTime = time.time()
    if not os.path.exists(args.input_sequence):
        logging.error(f"Couldn't find input file '{args.input_sequence}'")
        exit()

    if args.nt > mp.cpu_count():
        threads = mp.cpu_count()
    else:
        threads = args.nt

    sequences = {}
    try:
        sequences = SeqIO.to_dict(SeqIO.parse(args.input_sequence, args.sequence_format))
        if len(sequences.keys()) == 0:
            logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has incorrectly formatted sequences")
            exit()
    except ValueError:
        logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}'")
        exit()

    with open(args.output, 'w', newline='') as fileOutput:
        tsvWriter = csv.writer(fileOutput, delimiter='\t')

        bases = ["A", "C", "G", "T"]
        queries = [f"{a}{b}{c}" for c in bases for b in bases for a in bases]

        tsvWriter.writerow(["Library", "Sequence", "Start", "End"] + queries)

        for sequence in sequences:
            goldilocksSize = len(sequences[sequence])

            pool = mp.Pool(processes=threads)
            headers = sequence.split("_") if args.libraries else [os.path.basename(args.input_sequence), sequence]
            rows = [pool.apply(mpWindowCount, [sequences[sequence].seq.upper()[i:i+args.width], queries, i, i+args.width, headers]) for i in range(0, len(sequences[sequence]), args.spacing)]
            pool.join()
            tsvWriter.writerows(rows)

            logging.info(f"Sequence {sequence} windows written to output file")
    logging.info(f'Execution time in seconds: {time.time() - startTime}')
