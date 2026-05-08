import ruptures as rpt
import pandas as pd
import numpy as np
from Bio import Seq
from collections import Counter
import itertools

# Fast reverse complement using translation table
RC_TRANS = str.maketrans("ACGTacgt", "TGCAtgca")

def set_mers(mer_size=3):
    bases = ["A", "C", "G", "T"]
    return ["".join(a) for a in list(itertools.product(bases, repeat=mer_size))]

# Shorthand for string reverse-complement
def rc(sequence):
    if isinstance(sequence, Seq.Seq):
        sequence = str(sequence)
    return sequence.translate(RC_TRANS)[::-1]

def canonical(kmer):
    return min(kmer, rc(kmer))

def map_canonical_mers(kmers):
    kmer_index = {k: i for i, k in enumerate(kmers)}
    seen = set()
    result = {}
    for kmer, i in kmer_index.items():
        rc_mer = rc(kmer)
        if kmer in seen or rc_mer in seen:
            continue
        j = kmer_index[rc_mer]
        # canonical = lexicographically smaller of the pair
        result[canonical(kmer)] = [i, j]
        seen.add(kmer)
        seen.add(rc_mer)
    return result

def collapse_rc(row, mers, index=4, dim=None):
    if dim is None:
        dim = len(row) - index
    meta = row[:index]
    counts = row[index:index+dim]
    collapsed_counts = []
    for canon, (i, j) in mers.items():
        collapsed_counts.append(counts[i] + counts[j])
    return meta + collapsed_counts

# Validate takes an input df, and returns a tuple of (0) approved spectra descriptor fields and (1) triplet or query names
def validate(df):
    columns = df.columns
    validated_descriptors = []
    for descriptor in ['Library', 'Sequence', 'Block', 'Length', 'Start', 'End']:
        if descriptor in columns:
            validated_descriptors.append(descriptor)
    validated_names = [a for a in columns if a[0] in ['A', 'C', 'G', 'T']]
    return validated_descriptors, validated_names

# Simplify bi-directional frequencies/counts to unidirectional
def simplify(spectra, index=4, dim=None):
    if dim is None:
        dim = len(spectra.columns) - index
    simple_queries = {}
    cols = list(spectra.columns)[index:index + dim]
    for query in cols:
        rc_query = rc(query)
        if query not in simple_queries.values() and rc_query not in simple_queries:
            simple_queries[query] = rc_query

    for mer, rc_mer in simple_queries.items():
        if rc_mer != mer and rc_mer in spectra.columns:
            spectra[mer] += spectra[rc_mer]

    drop_cols = [rc_mer for mer, rc_mer in simple_queries.items() if rc_mer != mer and rc_mer in spectra.columns]
    return spectra.drop(columns=drop_cols)

# Counts the spectra of a sequence
def window_count(seq_task):
    window_seq, queries, start, end, headers = seq_task
    window_seq = str(window_seq)
    mer_len = len(queries[0]) if queries else 3
    counts = Counter(window_seq[i:i + mer_len] for i in range(len(window_seq) - mer_len + 1))
    return headers + [start + 1, min(start + len(window_seq), end)] + [counts.get(q, 0) for q in queries]

def window_count_no_overlap(seq_task):
    window_seq, queries, start, end, headers = seq_task
    window_seq = str(window_seq)
    return headers + [start + 1, min(start + len(window_seq), end)] + [window_seq.count(q) for q in queries]

# Calculate breakpoints from spectra
# For literal counts, a high penalty is ideal, but for frequencies a lower penalty is ideal
def get_breakpoints(spectra, index=4, dim=None, penalty=1000000, min_size=5):
    if dim is None:
        dim = len(spectra.columns) - index
    spectra_grouped = spectra.groupby(['Library', 'Sequence'])
    output = []
    # for each grouping, calculate rupture breakpoints, then convert them to starting window indices
    for name, group in spectra_grouped:
        data = group.iloc[:, index:index + dim].to_numpy()
        if len(data) > min_size * 2:
            data_algo = rpt.KernelCPD(min_size=min_size).fit(data).predict(pen=penalty)
            # data_algo contains the indices of the end of each segment
            # we want the 'Start' value of the window at those indices
            bkps = group.iloc[[a-1 for a in data_algo if a < len(group)], :]['Start'].tolist()
            output.append((name, bkps))
        else:
            output.append((name, []))
    return output

# Use breakpoints to append as bin column to spectra
def apply_breakpoints(spectra, breakpoints):
    spectra['Block'] = None
    for (lib, seq), bkps in breakpoints:
        break_count = 0
        current_start = 1

        mask_base = (spectra['Library'] == lib) & (spectra['Sequence'] == seq)

        if bkps:
            for bkp_start in bkps:
                break_name = f'{lib}_{seq}_{break_count:02d}'
                mask = mask_base & (spectra['Start'] >= current_start) & (spectra['End'] < bkp_start)
                spectra.loc[mask, 'Block'] = break_name
                current_start = bkp_start
                break_count += 1

            # Last block
            break_name = f'{lib}_{seq}_{break_count:02d}'
            mask = mask_base & (spectra['Start'] >= current_start)
            spectra.loc[mask, 'Block'] = break_name
        else:
            break_name = f'{lib}_{seq}_00'
            spectra.loc[mask_base, 'Block'] = break_name
    return spectra

def get_breakpoint_frequencies(spectra, frequency, index=4, dim=None):
    if dim is None:
        dim = len(spectra.columns) - index
    spectra_grouped = spectra.groupby(['Block', 'Library', 'Sequence'])
    outputs = []
    for (block, lib, seq), group in spectra_grouped:
        freqs = get_global_frequencies(group, frequency, index=index, dim=dim)
        row = {
            'Library': lib,
            'Sequence': seq,
            'Bin': block,
            'Length': max(group['End']) - min(group['Start']) + 1,
            'Start': min(group['Start']),
            'End': max(group['End'])
        }
        row.update(freqs)
        outputs.append(row)
    return pd.DataFrame(outputs)

# Transform spectra counts to spectra frequencies
def count_to_frequency(spectra, index=4, dim=None, mer_len=None):
    if dim is None:
        dim = len(spectra.columns) - index
    if mer_len is None:
        mer_len = len(spectra.columns[index])

    # Vectorized calculation
    denominator = spectra['End'] - spectra['Start'] - (mer_len - 2)
    # Avoid division by zero
    denominator = denominator.replace(0, np.nan)

    cols = spectra.columns[index:index+dim]
    spectra[cols] = spectra[cols].astype(float).div(denominator, axis=0).fillna(0)

    return spectra

# Transform spectra counts to spectra frequencies (Placeholder)
def frequency_to_count(spectra, index=4, dim=None, mer_len=None):
    return spectra

# Calculate global frequencies across spectra
def get_global_frequencies(spectra, frequency=False, index=4, dim=None):
    if dim is None:
        dim = len(spectra.columns) - index

    cols = spectra.columns[index:index + dim]
    if frequency:
        # If it's already frequencies, we average them?
        # Original: np.array(spectra.iloc[0:len(spectra), index:index + dim].sum()/len(spectra))
        return spectra[cols].mean().to_dict()
    else:
        # Sum of counts / total possible k-mers
        total_counts = spectra[cols].sum()
        mer_len = len(cols[0])
        total_length = (spectra['End'] - spectra['Start'] + 1).sum()
        # This denominator logic in original was: counts[3] - counts[2] - (merLen - 2) * len(spectra)
        # where counts[3] is sum(End), counts[2] is sum(Start).
        # sum(End - Start) - (merLen - 2) * N
        # which is sum(End - Start - (merLen - 2))
        denominator = (spectra['End'] - spectra['Start'] - (mer_len - 2)).sum()
        if denominator <= 0:
            return {c: 0 for c in cols}
        return (total_counts / denominator).to_dict()

def reduce_frequencies():
    return

def spectra_rc(spectra, index=4, dim=None, mer_len=None):
    if dim is None:
        dim = len(spectra.columns) - index
    if mer_len is None:
        mer_len = len(spectra.columns[index])
    new_spectra = spectra.copy()
    mers = list(spectra.columns)[index:index + dim]
    for mer in mers:
        rc_mer = rc(mer)
        if rc_mer in spectra.columns:
            new_spectra[rc_mer] = spectra[mer]
    return new_spectra
