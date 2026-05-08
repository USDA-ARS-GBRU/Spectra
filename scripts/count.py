#!/usr/bin/env python3

import os
import time
from Bio import SeqIO
import csv
import logging
import spectral
import multiprocessing
from collections import namedtuple

logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')
logger = logging.getLogger()

WindowTask = namedtuple("WindowTask", ["seq", "queries", "start", "end", "headers"])

def window_tasks_generator(input_sequence, sequence_format, queries_all, width, spacing, libraries, minimum_size, chunk_size):
    for seq_record in SeqIO.parse(input_sequence, sequence_format):
        sequence_name = seq_record.id
        headers = sequence_name.split("_") if libraries else [os.path.basename(input_sequence), sequence_name]
        sequence_length = len(seq_record)

        if sequence_length < minimum_size:
            continue

        # Process in chunks to save memory
        for chunk_start in range(0, sequence_length, chunk_size):
            chunk_end = min(chunk_start + chunk_size, sequence_length)
            sub_seq = str(seq_record.seq[chunk_start:chunk_end]).upper()

            for i in range(0, len(sub_seq), spacing):
                window = sub_seq[i:i+width]
                if not window or (len(window) < width and (chunk_start + i + len(window) < sequence_length)):
                    # Skip incomplete windows unless it's the very end of the sequence
                    continue

                start = chunk_start + i
                end = start + len(window)
                yield WindowTask(window, queries_all, start, end, headers)

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    start_time = time.time()
    if not os.path.exists(args.input_sequence):
        logger.error(f"Couldn't find input file '{args.input_sequence}'")
        return

    # Pre-calculate queries
    queries_all = spectral.set_mers(args.mer_size)
    if args.complement:
        queries = spectral.map_canonical_mers(queries_all)
    else:
        queries = {queries_all[a]: [a] for a in range(len(queries_all))}

    # Prepare for output
    tsv_headers = ["Library", "Sequence", "Start", "End"] + list(queries.keys())

    # Process windows in parallel
    callable_process = spectral.window_count if args.overlap else spectral.window_count_no_overlap

    try:
        with open(args.output, 'w', newline='') as file_output:
            tsv_writer = csv.writer(file_output, delimiter='\t')
            tsv_writer.writerow(tsv_headers)

            tasks = window_tasks_generator(
                args.input_sequence,
                args.sequence_format,
                queries_all,
                args.width,
                args.spacing,
                args.libraries,
                args.minimum_size,
                args.chunk_size
            )

            if args.threads > 1:
                pool = multiprocessing.Pool(processes=args.threads)
                # Use imap to maintain order and be memory efficient
                for row in pool.imap(callable_process, tasks):
                    if args.complement:
                        tsv_writer.writerow(spectral.collapse_rc(row, queries, dim=len(queries_all)))
                    else:
                        tsv_writer.writerow(row)
                pool.close()
                pool.join()
            else:
                for task in tasks:
                    row = callable_process(task)
                    if args.complement:
                        tsv_writer.writerow(spectral.collapse_rc(row, queries, dim=len(queries_all)))
                    else:
                        tsv_writer.writerow(row)
    except Exception as e:
        logger.error(f"Failed to write to {args.output}: {e}")
        return

    logger.info(f'Execution time in seconds: {time.time() - start_time}')
