#!/usr/bin/env python3
import os
import numpy as np
import scipy.stats as sp
import pandas as pd
import multiprocessing as mp

def groupByRows(group, rows):
    subGroups = group.groupby(group.index // rows)
    for subGroup in subGroups:
        subGroup[1][subGroup[1].keys()[4:]].agg('sum')
    return subGroups

def getGlobalFrequencies(spectra):
    counts = np.array(spectra.sum())
    return dict(zip(spectra.columns[4:], np.divide(counts[4:], counts[3] - counts[2] - len(spectra.Library))))

def normalize(row, frequencies, chiValue=1):
    if row[4:68].sum() != 0:
        localFrequencies = np.array(np.divide(row[4:68], row['End'] - row['Start'] - 1))
        if np.sum(localFrequencies) < 1:
            localFrequencies = np.divide(localFrequencies, np.sum(localFrequencies))
        chiResult = sp.chisquare(localFrequencies, f_exp=frequencies)
        return 1 if chiResult[0] > chiValue else 0
    else:
        return 1

def filterNormal(spectra, frequencies, file, chiValue=1, threads=0):
    spectra = spectra.assign(Normal=0)
    frequencies = np.array(np.divide(list(frequencies.values()), np.sum(list(frequencies.values()))))
    if threads:
        pool = mp.Pool(processes=threads)
        results = [pool.apply(normalize, [row, frequencies, chiValue]) for idx, row in spectra.iterrows()]
        pool.close()
        pool.join()
        spectra['Normal'] = results
    else:
        for row in spectra.iterrows():
            spectra.iloc[row[0], 68] = normalize(row[1], frequencies, chiValue)

    resultsOut = spectra.copy()
    resultsOut.iloc[spectra['Normal'] == 1, 4:68] = 0
    del resultsOut['Normal']
    resultsOut.to_csv(f"normal_{file}", sep='\t', index=False)
    spectra.iloc[spectra['Normal'] == 0, 4:68] = 0
    del spectra['Normal']
    spectra.to_csv(f"outlier_{file}", sep='\t', index=False)


def reduceFrequencies(spectra, frequencies):
    width = spectra['End'][0] - spectra['Start'][0] - 1
    for mer in frequencies:
        count = width * frequencies[mer]
        # testing multiple modes of reduction/normalization:
        # method one - (Fobs - Fexp) if (Fobs - Fexp) > 0 else 0
        # spectra[mer] = spectra[mer].apply(lambda x: (round(x - count) if x - count > 0 else 0))
        # method two - |Fobs - Fexp|
        # spectra[mer] = spectra[mer].apply(lambda x: abs(round(x - count)))
        # method three - Fobs - Fexp
        # spectra[mer] = spectra[mer].apply(lambda x: round(x - count))
        # method four - (Fobs - Fexp) / Fexp
        spectra[mer] = spectra[mer].apply(lambda x: round((x - count)/count))

    return spectra

def execute(args):
    if not os.path.exists(args.input_tsv):
        print("Could not find input file '" + args.input_tsv + "'")
        exit()
    if args.nt:
        if args.nt > mp.cpu_count():
            threads = mp.cpu_count()
        else:
            threads = args.nt
    else:
        threads = 0

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
                    spectra = pd.concat([spectra, newSubset])

    frequencies = {}
    if args.weighted_filter or args.weighted_normalization or args.verbose:
        frequencies = getGlobalFrequencies(spectra)

    if args.verbose:
        print('Reported frequencies (JSON/pydict):')
        print(frequencies)

    if args.weighted_filter:
        filterNormal(spectra, frequencies, args.output, threads=threads)

    if args.weighted_normalization:
        spectra = reduceFrequencies(spectra, frequencies)

    spectra.to_csv(args.output, sep='\t', index=False)
