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

def calculate_near_coincidence(df, threshold, range_n):
    # df should be pivot_df
    total_high_active = 0
    total_low_active = 0
    high_with_low_near = 0
    low_with_high_near = 0

    results_per_seq = {}

    for seq, group in df.groupby('Sequence'):
        group = group.sort_values('Start')
        h = (group['high'] > threshold).values
        l = (group['low'] > threshold).values

        n = len(group)
        h_near_l = np.zeros(n, dtype=bool)
        l_near_h = np.zeros(n, dtype=bool)

        for i in range(n):
            start_idx = max(0, i - range_n)
            end_idx = min(n, i + range_n + 1)

            if h[i]:
                if np.any(l[start_idx:end_idx]):
                    h_near_l[i] = True
            if l[i]:
                if np.any(h[start_idx:end_idx]):
                    l_near_h[i] = True

        num_h = np.sum(h)
        num_l = np.sum(l)
        num_h_near_l = np.sum(h_near_l)
        num_l_near_h = np.sum(l_near_h)

        total_high_active += num_h
        total_low_active += num_l
        high_with_low_near += num_h_near_l
        low_with_high_near += num_l_near_h

        results_per_seq[seq] = {
            'high_count': num_h,
            'low_count': num_l,
            'high_with_low_near': num_h_near_l,
            'low_with_high_near': num_l_near_h,
            'h_near_l_frac': num_h_near_l / num_h if num_h > 0 else 0,
            'l_near_h_frac': num_l_near_h / num_l if num_l > 0 else 0
        }

    global_results = {
        'total_high_active': total_high_active,
        'total_low_active': total_low_active,
        'high_with_low_near': high_with_low_near,
        'low_with_high_near': low_with_high_near,
        'h_near_l_frac': high_with_low_near / total_high_active if total_high_active > 0 else 0,
        'l_near_h_frac': low_with_high_near / total_low_active if total_low_active > 0 else 0
    }

    return global_results, results_per_seq

def plot_scatter(df, output_prefix):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(df['low'], df['high'], alpha=0.5, s=10)

    # Add diagonal line
    max_val = max(df['high'].max(), df['low'].max())
    ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.7)
#    plt.xscale("log")
#    plt.yscale("log")
    ax.set_xlabel('Low Extreme Kmer Count')
    ax.set_ylabel('High Extreme Kmer Count')
    ax.set_title('Scatter Plot: High vs Low Extreme Kmer Counts')

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_scatter.png", dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Kmer Mass Compare: Compare extreme kmer accumulations across genome windows")
    parser.add_argument('-i', '--input', required=True, help='Input TSV from mass-query.py')
    parser.add_argument('-o', '--output', default='mass_compare', help='Output prefix for plots and stats')
    parser.add_argument('--jaccard-minimum', type=int, default=0, help='Minimum count threshold for Jaccard and coincidence [default 0]')
    parser.add_argument('--near-range', type=int, default=1, help='Range N for near-coincidence statistics [default 1]')
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

    # Near coincidence
    global_near, seq_near = calculate_near_coincidence(pivot_df, args.jaccard_minimum, args.near_range)

    # Sequence-level stats
    seq_stats = []
    for seq, group in pivot_df.groupby('Sequence'):
        s_high = group['high'].sum()
        s_low = group['low'].sum()
        s_asym = (s_high - s_low) / (s_high + s_low) if (s_high + s_low) > 0 else 0
        s_len = group['End'].max() # Approx

        n_res = seq_near[seq]

        seq_stats.append({
            'Sequence': seq,
            'Length': s_len,
            'Total_High': s_high,
            'Total_Low': s_low,
            'Asymmetry_Index': s_asym,
            'High_per_kbp': (s_high / s_len) * 1000 if s_len > 0 else 0,
            'Low_per_kbp': (s_low / s_len) * 1000 if s_len > 0 else 0,
            f'High_with_Low_near_{args.near_range}': n_res['h_near_l_frac'],
            f'Low_with_High_near_{args.near_range}': n_res['l_near_h_frac']
        })

    # Write stats
    stats_file = f"{args.output}.stats"
    with open(stats_file, 'w') as f:
        f.write("# Global Statistics\n")
        f.write(f"Global_Total_High\t{total_high}\n")
        f.write(f"Global_Total_Low\t{total_low}\n")
        f.write(f"Global_Asymmetry_Index\t{global_asymmetry:.4f}\n")
        f.write(f"Pearson_Correlation\t{pearson_corr:.4f}\n")
        f.write(f"Jaccard_Coincidence_Index(>{args.jaccard_minimum})\t{jaccard:.4f}\n")
        f.write(f"Global_High_with_Low_near_{args.near_range}\t{global_near['h_near_l_frac']:.4f}\n")
        f.write(f"Global_Low_with_High_near_{args.near_range}\t{global_near['l_near_h_frac']:.4f}\n")
        f.write("\n# Per-Sequence Statistics\n")

    df_seq_stats = pd.DataFrame(seq_stats)
    df_seq_stats.to_csv(stats_file, sep='\t', index=False, mode='a')

    logger.info(f"Statistics written to {stats_file}")

    # Plotting
    plot_scatter(pivot_df, args.output)

    logger.info(f"Plots generated with prefix {args.output}")

if __name__ == "__main__":
    main()