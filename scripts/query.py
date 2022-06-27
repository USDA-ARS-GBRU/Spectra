#!/usr/bin/env python3

# TODO: create scalability for n-mer

import os
import time
from Bio import SeqIO
import csv
import logging
import spectral

logging.basicConfig(level=logging.INFO)

def execute(args):
    startTime = time.time()
    if not os.path.exists(args.input_sequence):
        logging.error(f"Couldn't find input file '{args.input_sequence}'")
        exit()
        
    if "," in args.query:
        queries = args.query.split(',')
    else:
        queries = [args.query]
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
        tsvWriter.writerow(["Library", "Sequence", "Start", "End"] + queries)

        for sequence in sequences:
            headers = sequence.split("_") if args.libraries else [os.path.basename(args.input_sequence), sequence]
            toProcess = [[sequences[sequence][i:i+args.width].seq.upper(), queries, i, i + args.width, headers] for i in range(0, len(sequences[sequence])+args.spacing, args.spacing) if len(sequences[sequence][i:i+args.width].seq) > 0]
            rows = map(spectral.windowCount, toProcess)
            tsvWriter.writerows(rows)

            logging.info(f"Sequence {sequence} windows written to output file")
    logging.info(f'Execution time in seconds: {time.time() - startTime}')


