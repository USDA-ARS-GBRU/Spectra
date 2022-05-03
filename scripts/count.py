#!/usr/bin/env python3

# TODO: create scalability for n-mer

import os

import numpy as np
from Bio import SeqIO
import csv
from math import ceil, floor
import re
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
    if not os.path.exists(args.input_sequence):
        logging.error(f"Couldn't find input file '{args.input_sequence}'")
        exit()

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
        queries = [a + b + c for c in bases for b in bases for a in bases]

        tsvWriter.writerow(["Library", "Sequence", "Start", "End"] + queries)

        for sequence in sequences:
            i = 0

            while i < len(sequences[sequence]):
                sequences[sequence]
                headers = sequence.split("_") if args.libraries else [os.path.basename(args.input_sequence), sequence]
                rowText = mpWindowCount(sequences[sequence].seq.upper()[i:i + args.width], queries, i, i + args.width, headers)
                tsvWriter.writerow(rowText)
                i += args.spacing
            headers = sequence.split("_") if args.libraries else [os.path.basename(args.input_sequence), sequence]
            rowText = mpWindowCount(sequences[sequence].seq.upper()[i:i + args.width], queries, i, i + args.width, headers)
            tsvWriter.writerow(rowText)

            logging.info(f"Sequence {sequence} windows written to output file")
