#!/usr/bin/env python3

# Project name: cluster
# Goal: spectra-directed contig-joining

# Input for this script will be a pre-ran spectra output.
# Program builds k-means clusters

import logging
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import spectral
from sklearn.cluster import KMeans
from sklearn.datasets import load_digits
from sklearn.decomposition import PCA
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)

    if not os.path.exists(args.input_tsv):
        logging.error(f"Could not find input file '{args.input_tsv}'")
        exit()

    spectra = pd.read_csv(args.input_tsv, delimiter='\t')
    bases = ["A", "C", "G", "T"]
    spectraLabels = [f"{a}{b}{c}" for c in bases for b in bases for a in bases]

    generations = 5
    current = 0
    outputData = pd.DataFrame(columns=spectra.columns)
    while current < generations:
        spectraData = spectra[spectra.columns.intersection(spectraLabels)]
        clusterSize = 2
        loop = True
        while loop:
            pca = PCA(2)
            df = pca.fit_transform(spectraData)
            kmeans = KMeans(n_clusters=clusterSize)
            label = kmeans.fit_predict(df)
            labelList = list(label)
            # shortest euclidian pairings will be the first pairs that appear as cluster total increases
            matches = [a for a in sorted(set(labelList)) if len([b for b in labelList if b == a]) == 2]
            if matches:
                outRows = []
                for match in matches:
                    matchingRows = [i for i in range(len(labelList)) if labelList[i] == match]
                    rows = spectra.iloc[matchingRows]
                    if list(rows['Sequence'])[0].split('_')[0] != list(rows['Sequence'])[1].split('_')[0]:
                        outRows = outRows + list(rows.index)
                if len(outRows) > 1:
                    loop = False
                else:
                    clusterSize += 1
            else:
                clusterSize += 1
        print(spectra.iloc[outRows])
        outputData = pd.merge(outputData, spectra.iloc[outRows])
        spectra = spectra.drop(outRows)
        current += 1
    # Writes this new pairing column to a new table.
    #spectra.insert(0, "Pairing", list(label), True)
    #if args.output_tsv:
    #    spectra.to_csv(args.output_tsv, sep='\t', index=False)
    #else:
    #    print(spectra)




