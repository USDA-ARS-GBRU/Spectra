import ruptures as rpt
import pandas as pd
import numpy as np

# Counts the spectra of a sequence
def windowCount(seq):
    return seq[4] + [seq[2] + 1, seq[2] + len(seq[0]) if seq[2] + len(seq[0]) < seq[3] else seq[3]] + [
               seq[0].count_overlap(a) for a in seq[1]]

# Calculate breakpoints from spectra
# For 64 literal counts, a penalty of 1,000,000 is ideal, but for frequencies a penalty of 0.5 is ideal
def getBreakpoints(spectra, index=4, length=64, penalty=1000000, min_size=5):
    spectra = spectra.groupby(['Library', 'Sequence'])
    output = []
    # for each grouping, calculate rupture breakpoints, then convert them to starting window indices
    for group in spectra:
        data = group[1].iloc[0:len(group[1]), index:index + length].to_numpy()
        dataAlgo = rpt.KernelCPD(min_size=min_size).fit(data).predict(pen=penalty)
        output.append((group[0], group[1].iloc[[a-1 for a in dataAlgo], 2:3]['Start'].tolist()))
    return output

# Use breakpoints to append as bin column to spectra
def applyBreakpoints(spectra, breakpoints):
    breakCount = 0
    for bkp in breakpoints:
        breakName = f'{bkp[0][0]}_{bkp[0][1]}_{breakCount}'
        index = 1
        for bkpStart in bkp[1]:
            indices = list(spectra.loc[(spectra['Library'] == bkp[0][0]) & (spectra['Sequence'] == bkp[0][1]) & (
                    spectra['Start'] >= index) & (spectra['End'] < bkpStart)].index)
            spectra.loc[indices, 'Bin'] = breakName
            index = bkpStart
        indices = list(spectra.loc[(spectra['Library'] == bkp[0][0]) & (spectra['Sequence'] == bkp[0][1]) & (
                    spectra['Start'] >= index)].index)
        spectra.loc[indices, 'Bin'] = breakName
    return spectra

def getBreakpointFrequencies(spectra):
    spectra = spectra.groupby(['Bin', 'Library', 'Sequence'])
    outputs = pd.DataFrame()
    for group in spectra:
        frequencies = getGlobalFrequencies(group[1])
        outputs.loc[len(outputs.index), ['Library', 'Sequence', 'Bin', 'Length', ] + list(frequencies.keys())] = [group[0][1], group[0][2], str(group[0][0]), str(max(group[1]['End']) - min(group[1]['Start']) + 1)] + list(frequencies.values())
    return outputs

# Transform spectra counts to spectra frequencies
def countToFrequency():
    return

# Transform spectra counts to spectra frequencies
def frequencyToCount():
    return

# Calculate global frequencies across spectra
def getGlobalFrequencies(spectra):
    counts = np.array(spectra.sum())
    return dict(zip(spectra.columns[4:68], np.divide(counts[4:68], counts[3] - counts[2] - len(spectra.Library))))

# Reduce frequencies from a given set of global frequencies
def reduceFrequencies():
    return
