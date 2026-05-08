#!/usr/bin/env python3

import os
import time
from Bio import SeqIO
import csv
import logging
import spectral
import pandas as pd
import multiprocessing
from collections import namedtuple

logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')
logger = logging.getLogger()

WindowTask = namedtuple("WindowTask", ["seq", "queries", "start", "end", "headers"])

def window_tasks_generator(sequences, queries, width, spacing, libraries, input_sequence):
    for sequence_id in sequences:
        headers = sequence_id.split("_") if libraries else [os.path.basename(input_sequence), sequence_id]
        seq_record = sequences[sequence_id]
        seq_str = str(seq_record.seq).upper()
        seq_len = len(seq_str)

        for i in range(0, seq_len, spacing):
            window = seq_str[i:i+width]
            if not window:
                continue
            yield WindowTask(window, queries, i, min(i+width, seq_len), headers)

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    start_time = time.time()
    if not os.path.exists(args.input_sequence):
        logger.error(f"Couldn't find input file '{args.input_sequence}'")
        return

    if "," in args.query:
        queries = [a.upper() for a in args.query.split(',')]
    else:
        queries = [args.query.upper()]

    try:
        if args.memory:
            sequences = SeqIO.index(args.input_sequence, args.sequence_format)
        else:
            sequences = SeqIO.to_dict(SeqIO.parse(args.input_sequence, args.sequence_format))

        if not sequences:
            logger.error(f"Sequence file '{args.input_sequence}' could not be loaded in format '{args.sequence_format}' or has no sequences")
            return
    except Exception as e:
        logger.error(f"Error loading sequence file '{args.input_sequence}': {e}")
        return

    if args.complement:
        new_queries = []
        for query in queries:
            query_rc = spectral.rc(query)
            if query not in new_queries:
                new_queries.append(query)
            if query_rc not in new_queries:
                new_queries.append(query_rc)
        queries = new_queries

    try:
        with open(args.output, 'w', newline='') as file_output:
            tsv_writer = csv.writer(file_output, delimiter='\t')
            tsv_writer.writerow(["Library", "Sequence", "Start", "End"] + queries)

            callable_process = spectral.window_count if args.overlap else spectral.window_count_no_overlap
            tasks = window_tasks_generator(sequences, queries, args.width, args.spacing, args.libraries, args.input_sequence)

            if args.threads > 1:
                pool = multiprocessing.Pool(processes=args.threads)
                for row in pool.imap(callable_process, tasks):
                    tsv_writer.writerow(row)
                pool.close()
                pool.join()
            else:
                for task in tasks:
                    tsv_writer.writerow(callable_process(task))
    except Exception as e:
        logger.error(f"Failed to write to {args.output}: {e}")
        return

    if args.complement:
        logger.info("Simplifying forward and r-c counts")
        try:
            spectra = pd.read_csv(args.output, delimiter='\t')
            spectra = spectral.simplify(spectra, dim=len(queries))
            spectra.to_csv(args.output, sep='\t', index=False)
        except Exception as e:
            logger.error(f"Error during simplification: {e}")

    logger.info(f'Execution time in seconds: {time.time() - start_time}')
