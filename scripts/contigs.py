#!/usr/bin/env python3

# Project name: contigs
# Goal: spectra-directed contig-joining

# Input for this script will be a pre-ran spectra output.
# Program works through each pairing of sequence's spectra profiles in order to find potential contact points through high-identity matches.
# In order to do this, spectra profiles will be:
    # joined together and then ruptures will attempt to break them apart iteratively by starting with relaxed settings
    # and iteratively using more stringent criteria. Spectra pairings will be scored by their "stickiness"

import ruptures as rpt
import logging
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import spectral
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    if not os.path.exists(args.input_tsv):
        logging.error(f"Could not find input file '{args.input_tsv}'")
        exit()

    spectra = pd.read_csv(args.input_tsv, delimiter='\t')
    # barebones core structure for placeholder work:

    #if multi-library spectra is given, warn that only the first library will be processed.
    # then throw out all the unimportant data to simplify computation and memory
    if len(spectra['Library'].unique()) > 1:
        logger.warning(f'Spectra profile contains more than one library. Only the first found in the file will be used for contig-joining.')
        spectra = spectra.loc[spectra['Library'] == spectra['Library'][0]]

    spectraSequences = spectra['Sequence'].unique()
    for sequenceA in range(0, len(spectraSequences)-1):
        sequenceDataA = spectra.loc[spectra['Sequence'] == spectraSequences[sequenceA]]
        for sequenceB in range(sequenceA+1, len(spectraSequences)):
            logging.info(f"Comparing sequences {spectraSequences[sequenceA]} and {spectraSequences[sequenceB]}")
            sequenceDataB = spectra.loc[spectra['Sequence'] == spectraSequences[sequenceB]]
            sequenceDataJoined = pd.concat([sequenceDataA, sequenceDataB])

            indexLength = 4
            spectraDimensions = len(spectra.columns) - indexLength

            #do sequence rename and point location renumbering
            sequenceDataJoined['Start'] = list(range(len(sequenceDataJoined['Start'])))
            sequenceDataJoined['End'] = list(range(len(sequenceDataJoined['End'])))
            sequenceDataJoined = sequenceDataJoined.assign(Sequence=f'5-{spectraSequences[sequenceA]}_5-{spectraSequences[sequenceB]}')
            breakpoints = spectral.getBreakpoints(sequenceDataJoined, penalty=1000000, index=indexLength, dim=spectraDimensions)
            print(breakpoints)

            sequenceDataJoined = pd.concat([sequenceDataB, sequenceDataA])
            #do sequence rename and point location renumbering
            sequenceDataJoined['Start'] = list(range(len(sequenceDataJoined['Start'])))
            sequenceDataJoined['End'] = list(range(len(sequenceDataJoined['End'])))
            sequenceDataJoined = sequenceDataJoined.assign(Sequence=f'5-{spectraSequences[sequenceB]}_5-{spectraSequences[sequenceA]}')
            breakpoints = spectral.getBreakpoints(sequenceDataJoined, penalty=1000000, index=indexLength, dim=spectraDimensions)
            print(breakpoints)

            sequenceDataJoined = pd.concat([spectral.spectraRC(sequenceDataA), sequenceDataB])
            #do sequence rename and point location renumbering
            sequenceDataJoined['Start'] = list(range(len(sequenceDataJoined['Start'])))
            sequenceDataJoined['End'] = list(range(len(sequenceDataJoined['End'])))
            sequenceDataJoined = sequenceDataJoined.assign(Sequence=f'3-{spectraSequences[sequenceA]}_5-{spectraSequences[sequenceB]}')
            breakpoints = spectral.getBreakpoints(sequenceDataJoined, penalty=1000000, index=indexLength, dim=spectraDimensions)
            print(breakpoints)

            sequenceDataJoined = pd.concat([sequenceDataA, spectral.spectraRC(sequenceDataB)])
            #do sequence rename and point location renumbering
            sequenceDataJoined['Start'] = list(range(len(sequenceDataJoined['Start'])))
            sequenceDataJoined['End'] = list(range(len(sequenceDataJoined['End'])))
            sequenceDataJoined = sequenceDataJoined.assign(Sequence=f'5-{spectraSequences[sequenceA]}_3-{spectraSequences[sequenceB]}')
            breakpoints = spectral.getBreakpoints(sequenceDataJoined, penalty=1000000, index=indexLength, dim=spectraDimensions)
            print(breakpoints)

