#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import os
import logging

# Logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger()

def calculate_jaccard(df, threshold):
    high_present = df['high'] > threshold
    low_present = df['low'] > threshold
    intersection = np.logical_and(high_present, low_present).sum()
    union = np.logical_or(high_present, low_present).sum()
    return intersection / union if union > 0 else 0.0

def plot_scatter(df, output_prefix):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(df['low'], df['high'], alpha=0.5, s=10)

    # Add diagonal line
    max_val = max(df['high'].max(), df['low'].max())
    ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.7)

    ax.set_xlabel('Low Extreme Kmer Count')
    ax.set_ylabel('High Extreme Kmer Count')
    ax.set_title('Scatter Plot: High vs Low Extreme Kmer Counts')

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_scatter.png", dpi=300)
    plt.close()

def plot_asymmetry_v_length(df_seq_stats, output_prefix):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(df_seq_stats['Length'], df_seq_stats['Asymmetry_Index'], alpha=0.7, color='green')
    ax.axhline(0, color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel('Asymmetry Index')
    ax.set_title('Sequence Asymmetry Index vs Length')
    ax.xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_asymmetry_v_length.png", dpi=300)
    plt.close()

def parse_ngaps(gff_path):
    # Simplified GFF parser for N-gaps
    gaps = {}
    if not gff_path or not os.path.exists(gff_path):
        return gaps

    try:
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                seqid = parts[0]
                start = int(parts[3])
                end = int(parts[4])

                if seqid not in gaps:
                    gaps[seqid] = []
                gaps[seqid].extend([start, end])

        # Sort for faster searching
        for seqid in gaps:
            gaps[seqid] = sorted(gaps[seqid])
    except Exception as e:
        logger.error(f"Error parsing N-gaps GFF: {e}")

    return gaps

def calculate_min_dist(row, seq_lengths, gaps):
    seq = row['Sequence']
    start = row['Start']
    end = row['End']
    length = seq_lengths.get(seq, end)

    # Initial distance to ends
    min_d = min(start - 1, length - end)

    # Check gaps
    if seq in gaps:
        # We can use binary search here for efficiency, but for now simple min
        for g_pos in gaps[seq]:
            # Distance from window to gap
            if g_pos < start:
                d = start - g_pos
            elif g_pos > end:
                d = g_pos - end
            else:
                d = 0 # Inside gap

            if d < min_d:
                min_d = d
                if min_d == 0:
                    break

    return min_d

def plot_density_v_distance(pivot_df, output_prefix, gaps):
    # Get sequence lengths
    seq_lengths = pivot_df.groupby('Sequence')['End'].max().to_dict()

    plot_data = pivot_df.copy()
    plot_data['DistFromEnd'] = plot_data.apply(calculate_min_dist, args=(seq_lengths, gaps), axis=1)

    # Density per kbp
    plot_data['WindowSize'] = plot_data['End'] - plot_data['Start'] + 1
    plot_data['High_per_kbp'] = (plot_data['high'] / plot_data['WindowSize']) * 1000
    plot_data['Low_per_kbp'] = (plot_data['low'] / plot_data['WindowSize']) * 1000

    fig, ax = plt.subplots(figsize=(12, 6))

    ax.scatter(plot_data['DistFromEnd'], plot_data['High_per_kbp'], alpha=0.4, color='red', s=5, label='High Bin')
    ax.scatter(plot_data['DistFromEnd'], plot_data['Low_per_kbp'], alpha=0.4, color='blue', s=5, label='Low Bin')

    ax.set_xlabel('Distance from Sequence End or N-gap (bp)')
    ax.set_ylabel('Extreme Kmer Density (counts per kbp)')
    ax.set_title('Extreme Kmer Density vs Distance from Nearest Feature (End or Gap)')
    ax.xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.legend()

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_density_v_distance.png", dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Kmer Mass Compare: Compare extreme kmer accumulations across genome windows")
    parser.add_argument('-i', '--input', required=True, help='Input TSV from mass-query.py')
    parser.add_argument('-o', '--output', default='mass_compare', help='Output prefix for plots and stats')
    parser.add_argument('--jaccard-minimum', type=int, default=0, help='Minimum count threshold for Jaccard coincidence [default 0]')
    parser.add_argument('--ngaps', help='GFF file of N-gap coordinates')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose mode')

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.INFO)

    if not os.path.exists(args.input):
        logger.error(f"Input file {args.input} not found.")
        return

    # Load data
    df = pd.read_csv(args.input, sep='\t')

    # Pivot to have one row per window
    pivot_df = df.pivot_table(index=['Sequence', 'Start', 'End'], columns='Bin', values='Count', fill_value=0).reset_index()

    # Ensure both bins exist
    if 'high' not in pivot_df.columns:
        pivot_df['high'] = 0
    if 'low' not in pivot_df.columns:
        pivot_df['low'] = 0

    # Global Statistics
    total_high = pivot_df['high'].sum()
    total_low = pivot_df['low'].sum()
    global_asymmetry = (total_high - total_low) / (total_high + total_low) if (total_high + total_low) > 0 else 0

    pearson_corr, _ = pearsonr(pivot_df['low'], pivot_df['high'])
    jaccard = calculate_jaccard(pivot_df, args.jaccard_minimum)

    # Sequence-level stats
    seq_stats = []
    for seq, group in pivot_df.groupby('Sequence'):
        s_high = group['high'].sum()
        s_low = group['low'].sum()
        s_asym = (s_high - s_low) / (s_high + s_low) if (s_high + s_low) > 0 else 0
        s_len = group['End'].max()

        seq_stats.append({
            'Sequence': seq,
            'Length': s_len,
            'Total_High': s_high,
            'Total_Low': s_low,
            'Asymmetry_Index': s_asym,
            'High_per_kbp': (s_high / s_len) * 1000 if s_len > 0 else 0,
            'Low_per_kbp': (s_low / s_len) * 1000 if s_len > 0 else 0
        })

    df_seq_stats = pd.DataFrame(seq_stats)

    # Write stats
    stats_file = f"{args.output}.stats"
    with open(stats_file, 'w') as f:
        f.write("# Global Statistics\n")
        f.write(f"Global_Total_High\t{total_high}\n")
        f.write(f"Global_Total_Low\t{total_low}\n")
        f.write(f"Global_Asymmetry_Index\t{global_asymmetry:.4f}\n")
        f.write(f"Pearson_Correlation\t{pearson_corr:.4f}\n")
        f.write(f"Jaccard_Coincidence_Index(>{args.jaccard_minimum})\t{jaccard:.4f}\n")
        if args.ngaps:
            f.write(f"NGaps_GFF\t{args.ngaps}\n")
        f.write("\n# Per-Sequence Statistics\n")

    df_seq_stats.to_csv(stats_file, sep='\t', index=False, mode='a')

    logger.info(f"Statistics written to {stats_file}")

    # Plotting
    plot_scatter(pivot_df, args.output)
    plot_asymmetry_v_length(df_seq_stats, args.output)

    gaps = parse_ngaps(args.ngaps) if args.ngaps else {}
    plot_density_v_distance(pivot_df, args.output, gaps)

    logger.info(f"Plots generated with prefix {args.output}")

if __name__ == "__main__":
    main()
