#!/usr/bin/env python3
import os
import numpy as np
import scipy.stats as sp
import pandas as pd
import logging
import spectral
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()


def groupByRows(group, rows, dim=64):
    subGroups = group.groupby(group.index // rows)
    for subGroup in subGroups:
        subGroup[1][subGroup[1].keys()[4:dim+4]].agg('sum')
    return subGroups

def normalize(row, frequencies, chiValue=1, dim=64):
    if row[4:dim+4].sum() != 0:
        localFrequencies = np.array(np.divide(row[4:dim+4], row['End'] - row['Start'] - 1))
        if np.sum(localFrequencies) < 1:
            localFrequencies = np.divide(localFrequencies, np.sum(localFrequencies))
        chiResult = sp.chisquare(localFrequencies, f_exp=frequencies)
        return 1 if chiResult[0] > chiValue else 0
    else:
        return 1

def filterNormal(spectra, frequencies, file, chiValue=1, dim=64):
    spectra = spectra.assign(Normal=0)
    frequencies = np.array(np.divide(list(frequencies.values()), np.sum(list(frequencies.values()))))
    for row in spectra.iterrows():
        spectra.iloc[row[0], 68] = normalize(row[1], frequencies, chiValue)
    resultsOut = spectra.copy()

    resultsOut.iloc[spectra['Normal'] == 1, 4:dim+4] = 0
    del resultsOut['Normal']
    resultsOut.to_csv(f"normal_{file}", sep='\t', index=False)
    spectra.iloc[spectra['Normal'] == 0, 4:dim+4] = 0
    del spectra['Normal']
    spectra.to_csv(f"outlier_{file}", sep='\t', index=False)

def reduceFrequencies(spectra, frequencies):
    width = spectra['End'][0] - spectra['Start'][0] - 1
    for mer in frequencies:
        count = width * frequencies[mer]
        spectra[mer] = spectra[mer].apply(lambda x: round((x - count)/count))

    return spectra

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    if not os.path.exists(args.input_tsv):
        logging.error(f"Could not find input file '{args.input_tsv}'")
        exit()

    spectra = pd.read_csv(args.input_tsv, delimiter='\t')
    if args.resize_window:
        targetFactor = int(args.resize_window)
        if targetFactor > 1:
            spectraGroups = spectra.groupby(['Library', 'Sequence'])
            spectra = pd.DataFrame()
            for group in spectraGroups:
                for index in range(0, len(group[1]), targetFactor):
                    subset = group[1].iloc[index:index+targetFactor]
                    newSubset = subset.sum().to_frame().transpose()
                    newSubset['Library'][0] = group[0][0]
                    newSubset['Sequence'][0] = group[0][1]
                    newSubset['Start'][0] = subset['Start'].min()
                    newSubset['End'][0] = subset['End'].max()
                    spectra = pd.concat([spectra, newSubset], ignore_index=True)
        spectra.reindex()
    frequencies = {}
    if args.weighted_filter or args.weighted_normalization or args.verbose:
        frequencies = spectral.getGlobalFrequencies(spectra)

    if args.print:
        logging.info('Reported frequencies (JSON/pydict):')
        logging.info(frequencies)

    if args.weighted_filter:
        filterNormal(spectra, frequencies, args.input_tsv)

    if args.weighted_normalization:
        spectra = reduceFrequencies(spectra, frequencies)

    if args.convert and args.frequencies:
        spectra = spectral.frequencyToCount(spectra)
    elif args.convert:
        spectra = spectral.countToFrequency(spectra)

    if args.simplify:
        spectra = spectral.simplify(spectra)

    if args.output:
        spectra.to_csv(args.output, sep='\t', index=False)
