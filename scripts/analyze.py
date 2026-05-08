#!/usr/bin/env python3

import os
import logging
import pandas as pd
import spectral
import matplotlib.pyplot as plt
import numpy as np
import ruptures as rpt

logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')
logger = logging.getLogger()

def gff_writer(results, output):
    try:
        with open(f"{output}_bins.gff", 'w') as out_file:
            for i, (_, line) in enumerate(results.iterrows()):
                strand = '+' if i % 2 == 0 else '-'
                out_file.write(f"{line['Sequence']}\tSpectra-bins\tbin-region\t{line['Start']}\t{line['End']}\t.\t{strand}\t.\tBin ID:{line['Bin']}, size: {line['Length']}\n")
    except Exception as e:
        logger.error(f"Error writing GFF file: {e}")

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    if not os.path.exists(args.input_tsv):
        logger.error(f"Could not find input file '{args.input_tsv}'")
        return

    try:
        spectra = pd.read_csv(args.input_tsv, delimiter='\t')
    except Exception as e:
        logger.error(f"Error reading {args.input_tsv}: {e}")
        return

    index_length = 4
    spectra_dimensions = len(spectra.columns) - index_length

    if args.frequency:
        if args.penalty > 1:
            args.penalty = 0.5
        spectra = spectral.count_to_frequency(spectra, index=index_length, dim=spectra_dimensions)

    breakpoints = spectral.get_breakpoints(spectra, penalty=args.penalty, min_size=args.size, index=index_length, dim=spectra_dimensions)
    spectra = spectral.apply_breakpoints(spectra, breakpoints)

    if args.output_prefix:
        try:
            spectra.to_csv(f"{args.output_prefix}.tsv", sep='\t', index=False)
        except Exception as e:
            logger.error(f"Error writing output TSV: {e}")

    # After applying breakpoints, 'Block' column might have been added, adjusting dimensions
    # Actually spectral.get_breakpoint_frequencies handles it if we pass correct dim
    results = spectral.get_breakpoint_frequencies(spectra, args.frequency, index=index_length, dim=spectra_dimensions)

    if args.output_prefix:
        try:
            results.to_csv(f"{args.output_prefix}_bins.tsv", index=False, sep='\t')
            gff_writer(results, args.output_prefix)
        except Exception as e:
            logger.error(f"Error writing results: {e}")
    else:
        for _, line in results.iterrows():
            print(f"{line['Library']}, {line['Sequence']}, {line['Start']}, {line['End']}")
