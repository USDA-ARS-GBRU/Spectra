#!/usr/bin/env python

import os
from SpectraClass import Spectra, writeToFile

# ('-i', '--input', dest='input_tsv', type=str, help='Input spectra tsv', required=True)
# ('-o', '--output', '--output', dest='output', type=str, help='Output spectra tsv', default='transformed_output.tsv')
# ('-r', '--weighted-rm', dest='weighted_removal', action='store_true', help='Remove windows with spectra frequencies similar to the whole sequence', default=False)
# ('-n', '--weighted-norm', dest='weighted_normalization', action='store_true', help='Normalize spectra frequencies for each window by the frequencies for the whole sequence', default=False)
# ('-s', '--window-resize', dest='resize_window', type=int, help='Image resolution in DPI', default=300)

def execute(args):
    if not os.path.exists(args.input_tsv):
        print("Could not find input file '" + args.input_tsv + "'")
        exit()

    spectra = Spectra(args.input_tsv)

    #spectra.setGlobalFrequencies()
    spectra.removeOutliers()
    writeToFile(spectra, "modified.tsv")
    #print(len(spectra.getMerSubset(library=['VmandF.fna'], coords=(1, 6000))))

    #print(spectra.globalFrequencies)
