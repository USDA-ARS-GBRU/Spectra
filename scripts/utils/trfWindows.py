#!/usr/bin/env python3

# given an input trf-tsv, generate a windowed repeat content tabulation

import logging
import argparse
import os
import time

import pandas
from Bio import SeqIO, Seq
import re
import pandas as pd
import csv

# Frame shift checks the outer bounds fo each window, then returns a content proportion of how many bases are within a repeat
def frameShift(values):
    if values[3] > values[5]:
        values[3] = values[5]
    totalLength = values[3] - values[2] + 1
    return [values[0], values[1], values[2], values[3], (sum(values[4])) / totalLength]
    #return rangeTable['start']-rangeTable['end']

# Set up logger
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

# Set up top level module argparser
parser = argparse.ArgumentParser(description='trfWindows: a tool for processing trf results for sliding-frame repeat content')
parser.add_argument('-i', '--input', dest='input_table', type=str, help='Input table file', required=True)
parser.add_argument('-f', '--fasta', dest='input_fasta', type=str, help='Input fasta file', required=True)
parser.add_argument('-o', '--output', dest='output_table', type=str, help='Output table file', required=True)
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
parser.add_argument('-w', '--width', dest='width', type=int, help='Window width', default='30000')
parser.add_argument('-s', '--spacing', dest='spacing', type=int, help='Window spacing', default='30000')
args = parser.parse_args()

if args.verbose:
    logger.setLevel(logging.INFO)

if not os.path.exists(args.input_table):
    logging.error(f"Couldn't find input file '{args.input_table}'")
    exit()
elif not os.path.exists(args.input_fasta):
    logging.error(f"Couldn't find input file '{args.input_fasta}'")
    exit()

sequences = {}
try:
    sequences = SeqIO.index(args.input_fasta, 'fasta')
    if len(sequences.keys()) == 0:
        logging.error(
            f"Sequence file '{args.input_fasta}' could not be loaded in fasta format or has incorrectly formatted sequences")
        exit()
except ValueError:
    logging.error(f"Sequence file '{args.input_fasta}' could not be loaded in fasta format")
    exit()

table = pd.read_csv(args.input_table, delimiter='\t').groupby('sequence')
with open(args.output_table, 'w', newline='') as fileOutput:
    tsvWriter = csv.writer(fileOutput, delimiter='\t')
    tsvWriter.writerow(['Library', 'Sequence', 'Start', 'End', 'Proportion'])
    for group in table:
        sequenceLength = len(sequences[group[0]])

        sequenceHits = [0] * sequenceLength
        for row in group[1].iterrows():
            sequenceHits[row[1]['start']:row[1]['end']] = [1] * (row[1]['end'] - row[1]['start'] + 1)

        # frameShift needs: sequence name, start, end, subset of table with repeats
        toProcess = [[args.input_fasta, group[0], i, i + args.width - 1, sequenceHits[i:(i + args.width)], sequenceLength] for i in range(1, sequenceLength+args.spacing, args.spacing) if sequenceLength >= i]
        #toProcess = [[args.input_fasta, group[0], i, i + args.width - 1, group[1].loc[group[1]['start'] >= i].loc[group[1]['end'] < i+args.width], sequenceLength] for i in range(1, sequenceLength+args.spacing, args.spacing) if sequenceLength >= i]
        rows = map(frameShift, toProcess)
        tsvWriter.writerows(rows)
