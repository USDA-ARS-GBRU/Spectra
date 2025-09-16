#!/usr/bin/env python3

# TODO: create scalability for n-mer

import os
import time
from Bio import SeqIO, Seq
import csv
import logging
import spectral
import pandas as pd
import itertools
from collections import namedtuple, Counter

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

WindowTask = namedtuple("WindowTask", ["seq", "queries", "start", "end", "headers"])
def window_tasks(sequence_record, queries, width, spacing, headers, offset=0):
    seq_str = str(sequence_record.seq).upper()
    seq_len = len(seq_str)

    for i in range(0, seq_len, spacing):
        window = seq_str[i:i+width]
        if not window:
            continue
        yield WindowTask(window, queries, i+offset, min(i+width+offset, seq_len+offset), headers)


def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    startTime = time.time()
    if not os.path.exists(args.input_sequence):
        logging.error(f"Couldn't find input file '{args.input_sequence}'")
        exit()

    sequences = SeqIO.index(args.input_sequence, args.sequence_format)
    try:
        if len(sequences.keys()) == 0:
            logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has incorrectly formatted sequences")
            exit()
    except ValueError:
        logging.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}'")
        exit()

    queries_all=spectral.setMers(args.mer_size)
    rep_map = {q: spectral.canonical(q) for q in queries_all}
    #queries = ["".join(a) for a in list(itertools.product(bases, repeat=args.mer_size))]
    if args.complement:
        queries = spectral.mapCanonicalMers(queries_all)
    else:
        queries = {queries_all[a]: [a] for a in range(len(queries_all))}

    with open(args.output, 'w', newline='') as fileOutput:
        tsvHeaders = ["Library", "Sequence", "Start", "End"] + list(queries.keys())
        tsvWriter = csv.writer(fileOutput, delimiter='\t')
        tsvWriter.writerow(tsvHeaders)
        callableProcess = spectral.windowCount if args.overlap else spectral.windowCountNoOverlap

        for sequence_name, seq_record in sequences.items():
            headers = sequence_name.split("_") if args.libraries else [os.path.basename(args.input_sequence),
                                                                       sequence_name]
            sequenceLength = len(seq_record)
            if sequenceLength >= args.chunk_size:
                logging.info(f"Sequence {sequence_name} is large. Breaking into smaller segments")
                indices = range(0, sequenceLength, args.chunk_size)
                for sequenceIndex in indices:
                    sub_seq = seq_record[sequenceIndex:sequenceIndex + args.chunk_size]
                    for task in window_tasks(sub_seq, queries_all, args.width, args.spacing, headers, offset=sequenceIndex):
                        row = callableProcess(task)
                        if args.complement:
                            tsvWriter.writerow(spectral.collapseRC(row, queries))
                        else:
                            tsvWriter.writerow(row)
                logging.info(f"Sequence {sequence_name} windows written to output file")
            else:
                for task in window_tasks(seq_record, queries_all, args.width, args.spacing, headers):
                    row = callableProcess(task)
                    if args.complement:
                        tsvWriter.writerow(spectral.collapseRC(row, queries))
                    else:
                        tsvWriter.writerow(row)
                logging.info(f"Sequence {sequence_name} windows written to output file")

    #if args.complement:
    #    logging.info("Simplifying forward and r-c counts")
    #    spectra = pd.read_csv(args.output, delimiter='\t')
    #    spectra = spectral.simplify(spectra, dim=len(queries))
    #    spectra.to_csv(args.output, sep='\t', index=False)

    logging.info(f'Execution time in seconds: {time.time() - startTime}')
