#!/usr/bin/env python3

# given an input trf-tsv, generate a windowed repeat content tabulation

import logging
import argparse
import os
import time
from Bio import SeqIO
import pandas as pd
import csv

# Frame shift checks the outer bounds of each window, then returns a content proportion of how many bases are within a repeat
def frame_shift(values):
    # values: [input_fasta, seq_name, start, end, sequenceHits_subset, sequenceLength]
    start = values[2]
    end = values[3]
    seq_len = values[5]
    if end > seq_len:
        end = seq_len

    actual_length = end - start + 1
    if actual_length <= 0:
        proportion = 0
    else:
        proportion = sum(values[4]) / actual_length

    return [values[0], values[1], start, end, proportion]

def main():
    # Set up top level module argparser
    parser = argparse.ArgumentParser(description='trfWindows: a tool for processing trf results for sliding-frame repeat content')
    parser.add_argument('-i', '--input', dest='input_table', type=str, help='Input table file', required=True)
    parser.add_argument('-f', '--fasta', dest='input_fasta', type=str, help='Input fasta file', required=True)
    parser.add_argument('-o', '--output', dest='output_table', type=str, help='Output table file', required=True)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
    parser.add_argument('-w', '--width', dest='width', type=int, help='Window width', default=30000)
    parser.add_argument('-s', '--spacing', dest='spacing', type=int, help='Window spacing', default=30000)
    args = parser.parse_args()

    # Set up logger
    logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')
    logger = logging.getLogger()
    if args.verbose:
        logger.setLevel(logging.INFO)

    if not os.path.exists(args.input_table):
        logger.error(f"Couldn't find input file '{args.input_table}'")
        return
    elif not os.path.exists(args.input_fasta):
        logger.error(f"Couldn't find input file '{args.input_fasta}'")
        return

    try:
        sequences = SeqIO.index(args.input_fasta, 'fasta')
        if not sequences:
            logger.error(f"Sequence file '{args.input_fasta}' could not be loaded or is empty")
            return
    except Exception as e:
        logger.error(f"Error loading sequence file '{args.input_fasta}': {e}")
        return

    try:
        table_groups = pd.read_csv(args.input_table, delimiter='\t').groupby('sequence')
    except Exception as e:
        logger.error(f"Error reading input table '{args.input_table}': {e}")
        return

    try:
        with open(args.output_table, 'w', newline='') as file_output:
            tsv_writer = csv.writer(file_output, delimiter='\t')
            tsv_writer.writerow(['Library', 'Sequence', 'Start', 'End', 'Proportion'])

            for seq_name, group in table_groups:
                if seq_name not in sequences:
                    logger.warning(f"Sequence '{seq_name}' found in table but not in FASTA")
                    continue

                sequence_length = len(sequences[seq_name])
                sequence_hits = [0] * (sequence_length + 1)

                for _, row in group.iterrows():
                    # Ensure indices are within bounds
                    hit_start = max(0, int(row['start']))
                    hit_end = min(sequence_length, int(row['end']))
                    for i in range(hit_start, hit_end + 1):
                        sequence_hits[i] = 1

                for i in range(1, sequence_length + 1, args.spacing):
                    end_idx = min(i + args.width - 1, sequence_length)
                    hits_subset = sequence_hits[i:end_idx + 1]

                    row_data = frame_shift([args.input_fasta, seq_name, i, end_idx, hits_subset, sequence_length])
                    tsv_writer.writerow(row_data)
    except Exception as e:
        logger.error(f"Error writing output table '{args.output_table}': {e}")

if __name__ == "__main__":
    main()
