#!/usr/bin/env python3

import os
import csv
import logging
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    with open(args.output, 'w') as outputFile:
        headers = 0
        tsvWriter = csv.writer(outputFile, delimiter='\t')
        for spectraInput in args.input_tsvs:
            if not os.path.exists(spectraInput):
                logging.error(f"Could not find input file '{spectraInput}', skipping this file")
            else:
                with open(spectraInput, 'r') as inputFile:
                    tsvReader = csv.reader(inputFile, delimiter='\t')
                    if not headers:
                        headers = next(tsvReader)
                        tsvWriter.writerow(headers)
                    else:
                        next(tsvReader)
                    for spectraWindow in tsvReader:
                        tsvWriter.writerow(spectraWindow)
                logging.info(f"{spectraInput} collated")
