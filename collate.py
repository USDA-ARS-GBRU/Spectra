#!/usr/bin/env python

import os
import csv

def execute(args):
    with open(args.output, 'w') as outputFile:
        headers = 0
        tsvWriter = csv.writer(outputFile, delimiter='\t')
        for spectraInput in args.input_tsvs:
            if not os.path.exists(spectraInput):
                print("Could not find input file '"+spectraInput+"', skipping this file")
            else:
                with open(spectraInput, 'r') as inputFile:
                    tsvReader = csv.reader(inputFile, delimiter='\t')
                    if not headers:
                        headers = next(tsvReader)
                        tsvWriter.writerow(headers)
                    for spectraWindow in tsvReader:
                        tsvWriter.writerow(spectraWindow)
