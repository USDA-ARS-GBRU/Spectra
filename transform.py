#!/usr/bin/env python3

import os
from SpectraClass import Spectra, writeToFile

def execute(args):
    if not os.path.exists(args.input_tsv):
        print("Could not find input file '" + args.input_tsv + "'")
        exit()

    spectra = Spectra(args.input_tsv, frequencies=args.frequencies)
    if args.resize_window:
        spectra.resizeWindows(args.resize_window)

    if args.weighted_removal and args.weighted_normalization:
        print("Error: cannot process both commands, normalizing Spectra only")
        spectra.setGlobalFrequencies()
        spectra.reduceByGlobalFrequencies()
    elif args.weighted_normalization:
        spectra.setGlobalFrequencies()
        spectra.reduceByGlobalFrequencies()
    elif args.weighted_removal:
        spectra.setGlobalFrequencies()
        spectra.removeOutliers()

    writeToFile(spectra, args.output)
