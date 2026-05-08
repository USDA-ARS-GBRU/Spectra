#!/usr/bin/env python3

import os
import csv
import logging

logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')
logger = logging.getLogger()

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    try:
        with open(args.output, 'w', newline='') as output_file:
            headers = None
            tsv_writer = csv.writer(output_file, delimiter='\t')
            for spectra_input in args.input_tsvs:
                if not os.path.exists(spectra_input):
                    logger.error(f"Could not find input file '{spectra_input}', skipping")
                    continue

                with open(spectra_input, 'r') as input_file:
                    tsv_reader = csv.reader(input_file, delimiter='\t')
                    try:
                        current_headers = next(tsv_reader)
                    except StopIteration:
                        continue

                    if headers is None:
                        headers = current_headers
                        tsv_writer.writerow(headers)

                    for row in tsv_reader:
                        tsv_writer.writerow(row)
                logger.info(f"{spectra_input} collated")
    except Exception as e:
        logger.error(f"Error during collation: {e}")
