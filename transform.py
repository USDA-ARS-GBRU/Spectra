#!/usr/bin/env python3
import os
import numpy as np
import scipy.stats as sp
import pandas as pd
# from SpectraClass import Spectra, writeToFile

def groupByRows(group, rows):
    subGroups = group.groupby(group.index // rows)
    for subGroup in subGroups:
        subGroup[1][subGroup[1].keys()[4:]].agg('sum')
    return subGroups

def getGlobalFrequencies(spectra):
    counts = np.array(spectra.sum())
    return dict(zip(spectra.columns[4:], np.divide(counts[4:], counts[3] - counts[2] - len(spectra.Library))))

def filterNormal(spectra, frequencies, chiValue=1, normal=False):
    frequencies = np.array(np.divide(list(frequencies.values()), np.sum(list(frequencies.values()))))
    for row in spectra.iterrows():
        if row[1][4:].sum() == 0:
            spectra.drop(row[0], inplace=True)
        else:
            localFrequencies = np.array(np.divide(row[1][4:], row[1]['End']-row[1]['Start']-1))
            if np.sum(localFrequencies) < 1:
                localFrequencies = np.divide(localFrequencies, np.sum(localFrequencies))
            chiResult = sp.chisquare(localFrequencies, f_exp=frequencies)
            if normal:
                if chiResult[0] <= chiValue:
                    spectra.drop(row[0], inplace=True)
            else:
                if chiResult[0] > chiValue:
                    spectra.drop(row[0], inplace=True)
    return spectra

def reduceFrequencies(spectra, frequencies):
    for mer in frequencies:
        spectra[mer] = spectra[mer].apply(lambda x: round(x * (1-frequencies[mer])))
    return spectra

def execute(args):
    if not os.path.exists(args.input_tsv):
        print("Could not find input file '" + args.input_tsv + "'")
        exit()

    # spectra = Spectra(args.input_tsv, frequencies=args.frequencies)
    # if args.resize_window:
    #     spectra.resizeWindows(args.resize_window)
    #
    # if args.weighted_removal and args.weighted_normalization:
    #     print("Error: cannot process both commands, normalizing Spectra only")
    #     spectra.setGlobalFrequencies()
    #     spectra.reduceByGlobalFrequencies()
    # elif args.weighted_normalization:
    #     spectra.setGlobalFrequencies()
    #     spectra.reduceByGlobalFrequencies()
    # elif args.weighted_removal:
    #     spectra.setGlobalFrequencies()
    #     spectra.removeOutliers()
    #
    # writeToFile(spectra, args.output)

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
    if args.weighted_removal or args.weighted_normalization or args.weighted_keep or args.verbose:
        frequencies = getGlobalFrequencies(spectra)

    if args.verbose:
        print('Reported frequencies (JSON/pydict):')
        print(frequencies)

    if args.weighted_normalization:
        spectra = reduceFrequencies(spectra, frequencies)
    elif args.weighted_removal:
        spectra = filterNormal(spectra, frequencies, normal=False)
    elif args.weighted_keep:
        spectra = filterNormal(spectra, frequencies, normal=True)

    spectra.to_csv(args.output, sep='\t', index=False)
