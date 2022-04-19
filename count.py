#!/usr/bin/env python

# TODO: create scalability for n-mer

import os
from Bio import SeqIO
import csv
from math import ceil, floor

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
        print("Couldn't find input file '"+args.input_sequence)
        exit()
    try:
        sequences = SeqIO.to_dict(SeqIO.parse(args.input_sequence, args.sequence_format))
    except ValueError:
        print("Sequence file '"+args.input_sequence+"' could not be loaded in format '"+args.sequence_format+"'")
        exit()

    if len(sequences.keys()) == 0:
        print("Sequence file '"+args.input_sequence+"' could not be loaded in format '"+args.sequence_format+"' or has incorrectly formatted sequences")
        exit()

    with open(args.output, 'w', newline='') as fileOutput:
        tsvWriter = csv.writer(fileOutput, delimiter='\t')

        window_lead = floor(args.width / 2)
        window_tail = ceil(args.width / 2)

        bases = ["A", "C", "G", "T"]
        queries = [a + b + c for c in bases for b in bases for a in bases]

        tsvWriter.writerow(["Library", "Sequence", "Start", "End"] + queries)

        for sequence in sequences:
            i = window_lead
            while i + window_tail < len(sequences[sequence]):
                rowCount = windowCount(str(sequences[sequence][i - window_lead:i + window_tail].seq).upper(), queries, args.proportions)
                headers = sequence.split("_") if args.libraries else [os.path.basename(args.input_sequence), sequence]
                rowText = headers + [i - window_lead +1, i + window_tail] + rowCount
                tsvWriter.writerow(rowText)
                i += args.spacing
            rowCount = windowCount(str(sequences[sequence][i - window_lead:i + window_tail].seq).upper(), queries, args.proportions)
            headers = sequence.split("_") if args.libraries else [os.path.basename(args.input_sequence), sequence]
            rowText = headers + [i - window_lead + 1, i + window_tail] + rowCount
            tsvWriter.writerow(rowText)