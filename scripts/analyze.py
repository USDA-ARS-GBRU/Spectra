#!/usr/bin/env python3

# provide info about a given spectra file

import ruptures as rpt
import logging
import os
import pandas as pd
import matplotlib.pyplot as plt
import spectral
import numpy as np
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

# plots breakpoints calculated through ruptures
def plotBreakpoints(spectra, index=4, dim=64, penalty=1000000, min_size=5, output=""):
    spectra = spectra.groupby(['Library', 'Sequence'])
    for group in spectra:
        data = group[1].iloc[0:len(group[1]), index:index + dim].to_numpy()
        if len(data) > min_size * 2:
            dataAlgo = rpt.KernelCPD(min_size=min_size).fit(data).predict(pen=penalty)
            fig, ax_arr = rpt.display(np.ndarray(shape=1), dataAlgo, figsize=(30, 1))
            plt.savefig(f'breakpoints_breakdown_{output}_{group[0][0]}_{group[0][1]}.png')
            plt.close()

def gffWriter(results, output):
    with open(f"{output}_bins.gff", 'w') as outFile:
        i = 0
        for line in results.iterrows():
            outFile.write(f"{line[1]['Sequence']}\tSpectra-bins\tbin-region\t{line[1]['Start']}\t{line[1]['End']}\t.\t{'+' if i%2==0 else '-'}\t.\tBin ID:{line[1]['Bin']}, size: {line[1]['Length']}\n")
            i += 1

# function to iterate over the windows and convert them to frequencies
# raw counts misrepresent kmer diversity in gappy alignments
def padWindows(spectra, tally, mer):
    for index, row in spectra.iterrows():
        spectra.loc[index, mer] = [a/tally[index] for a in row[mer]]
    return spectra

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    if not os.path.exists(args.input_tsv):
        logging.error(f"Could not find input file '{args.input_tsv}'")
        exit()

    spectra = pd.read_csv(args.input_tsv, delimiter='\t')
    indexLength = 4
    spectraDimensions = len(spectra.columns) - indexLength
    if args.frequency:
        # penalty is temporarily locked to 1 or lower for frequency-based breakpoints
        if args.penalty > 1:
            args.penalty = 0.5
        spectra = spectral.countToFrequency(spectra, index=indexLength, dim=spectraDimensions)
    breakpoints = spectral.getBreakpoints(spectra, penalty=args.penalty, min_size=args.size, index=indexLength, dim=spectraDimensions)
    spectra = spectral.applyBreakpoints(spectra, breakpoints)

    if args.output_prefix:
        spectra.to_csv(f"{args.output_prefix}.tsv", sep='\t', index=False)

    spectraDimensions -= 1
    results = spectral.getBreakpointFrequencies(spectra, args.frequency, index=indexLength, dim=spectraDimensions)
    if args.output_prefix:
        results.to_csv(f"{args.output_prefix}_bins.tsv", index=False, sep='\t')
        gffWriter(results,args.output_prefix)
    else:
        for line in results.iterrows():
            print(f"{line[1][0]}, {line[1][1]}, {line[1][4]}, {line[1][5]}")
    #plotBreakpoints(spectra, penalty=args.penalty, min_size=args.size, index=indexLength, dim=spectraDimensions, output=args.output_tsv if args.output_tsv else '')
    #