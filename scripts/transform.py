#!/usr/bin/env python3

import os
import numpy as np
import scipy.stats as sp
import pandas as pd
import logging
import spectral

logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')
logger = logging.getLogger()

def normalize(row, frequencies, chi_value=1, dim=64):
    if row.iloc[4:dim+4].sum() != 0:
        local_frequencies = np.array(np.divide(row.iloc[4:dim+4], row['End'] - row['Start'] - 1))
        if np.sum(local_frequencies) < 1:
            local_frequencies = np.divide(local_frequencies, np.sum(local_frequencies))
        chi_result = sp.chisquare(local_frequencies, f_exp=frequencies)
        return 1 if chi_result[0] > chi_value else 0
    else:
        return 1

def filter_normal(spectra, frequencies, file_path, chi_value=1, dim=64):
    spectra = spectra.assign(Normal=0)
    freq_values = np.array(list(frequencies.values()))
    freq_normalized = freq_values / np.sum(freq_values)

    for idx, row in spectra.iterrows():
        spectra.at[idx, 'Normal'] = normalize(row, freq_normalized, chi_value, dim)

    results_out = spectra.copy()
    results_out.loc[spectra['Normal'] == 1, spectra.columns[4:dim+4]] = 0
    results_out.drop(columns=['Normal']).to_csv(f"normal_{os.path.basename(file_path)}", sep='\t', index=False)

    spectra.loc[spectra['Normal'] == 0, spectra.columns[4:dim+4]] = 0
    spectra.drop(columns=['Normal']).to_csv(f"outlier_{os.path.basename(file_path)}", sep='\t', index=False)

def reduce_frequencies_weighted(spectra, frequencies):
    # This assumes all windows have the same width?
    width = spectra['End'].iloc[0] - spectra['Start'].iloc[0] - 1
    for mer, freq in frequencies.items():
        count = width * freq
        if count != 0:
            spectra[mer] = spectra[mer].apply(lambda x: round((x - count)/count))
    return spectra

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

    if args.resize_window:
        target_factor = int(args.resize_window)
        if target_factor > 1:
            spectra_groups = spectra.groupby(['Library', 'Sequence'])
            new_spectra_list = []
            for name, group in spectra_groups:
                for i in range(0, len(group), target_factor):
                    subset = group.iloc[i:i+target_factor]
                    new_row = subset.sum()
                    new_row['Library'] = name[0]
                    new_row['Sequence'] = name[1]
                    new_row['Start'] = subset['Start'].min()
                    new_row['End'] = subset['End'].max()
                    new_spectra_list.append(new_row)
            spectra = pd.DataFrame(new_spectra_list)
        spectra.reset_index(drop=True, inplace=True)

    frequencies = {}
    if args.weighted_filter or args.weighted_normalization or args.print:
        frequencies = spectral.get_global_frequencies(spectra)

    if args.print:
        logger.info('Reported frequencies:')
        logger.info(frequencies)

    if args.weighted_filter:
        filter_normal(spectra, frequencies, args.input_tsv)

    if args.weighted_normalization:
        spectra = reduce_frequencies_weighted(spectra, frequencies)

    if args.convert and args.frequencies:
        spectra = spectral.frequency_to_count(spectra)
    elif args.convert:
        spectra = spectral.count_to_frequency(spectra)

    if args.simplify:
        spectra = spectral.simplify(spectra)

    if args.output:
        try:
            spectra.to_csv(args.output, sep='\t', index=False)
        except Exception as e:
            logger.error(f"Error writing to {args.output}: {e}")
