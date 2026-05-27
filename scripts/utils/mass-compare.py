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
    # GFF parser for N-gaps returning intervals
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
                gaps[seqid].append((start, end))

        # Sort intervals
        for seqid in gaps:
            gaps[seqid] = sorted(gaps[seqid])
    except Exception as e:
        logger.error(f"Error parsing N-gaps GFF: {e}")

    return gaps

def get_gap_overlap(w_start, w_end, seq_gaps):
    overlap = 0
    for g_start, g_end in seq_gaps:
        o_start = max(w_start, g_start)
        o_end = min(w_end, g_end)
        if o_start <= o_end:
            overlap += (o_end - o_start + 1)
    return overlap

def calculate_min_dist_v2(w_start, w_end, length, seq_gaps):
    # Initial distance to ends
    min_d = min(w_start - 1, length - w_end)

    # Check gaps
    for g_start, g_end in seq_gaps:
        # Distance from window to gap
        if g_end < w_start:
            d = w_start - g_end - 1
        elif g_start > w_end:
            d = g_start - w_end - 1
        else:
            d = 0 # Overlaps or touches gap

        if d < min_d:
            min_d = d
            if min_d == 0:
                break
    return max(0, min_d)

def plot_density_v_distance(plot_data, output_prefix):
    fig, ax = plt.subplots(figsize=(12, 6))

    ax.scatter(plot_data['MinDist'], plot_data['High_per_kbp'], alpha=0.4, color='red', s=5, label='High Bin')
    ax.scatter(plot_data['MinDist'], plot_data['Low_per_kbp'], alpha=0.4, color='blue', s=5, label='Low Bin')

    ax.set_xlabel('Distance from Sequence End or N-gap (bp)')
    ax.set_ylabel('Extreme Kmer Density (counts per kbp)')
    ax.set_title('Extreme Kmer Density vs Distance from Nearest Feature (End or Gap)')
    ax.xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.legend()

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_density_v_distance.png", dpi=300)
    plt.close()

def plot_end_comparison(stats, output_prefix):
    labels = ['High (Background)', 'High (End)', 'Low (Background)', 'Low (End)']
    values = [
        stats['High_Density_Background'], stats['High_Density_End'],
        stats['Low_Density_Background'], stats['Low_Density_End']
    ]
    colors = ['#ffcccc', 'red', '#ccccff', 'blue']

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(labels, values, color=colors)
    ax.set_ylabel('Mean Density (counts per kbp)')
    ax.set_title('Comparison of Extreme Kmer Densities: End vs Background')

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_end_comparison.png", dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Kmer Mass Compare: Compare extreme kmer accumulations across genome windows")
    parser.add_argument('-i', '--input', required=True, help='Input TSV from mass-query.py')
    parser.add_argument('-o', '--output', default='mass_compare', help='Output prefix for plots and stats')
    parser.add_argument('--jaccard-minimum', type=int, default=0, help='Minimum count threshold for Jaccard coincidence [default 0]')
    parser.add_argument('--ngaps', help='GFF file of N-gap coordinates')
    parser.add_argument('--end-threshold', type=int, default=0, help='Distance threshold for "end-associated" windows [default 0]')
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

    gaps = parse_ngaps(args.ngaps) if args.ngaps else {}
    seq_lengths = pivot_df.groupby('Sequence')['End'].max().to_dict()

    # Process all windows
    processed_windows = []
    for idx, row in pivot_df.iterrows():
        seq = row['Sequence']
        w_start = row['Start']
        w_end = row['End']

        seq_gaps = gaps.get(seq, [])
        gap_len = get_gap_overlap(w_start, w_end, seq_gaps)
        w_len = w_end - w_start + 1
        effective_len = w_len - gap_len

        if effective_len <= 0:
            continue # Ignore window entirely comprised of N

        min_dist = calculate_min_dist_v2(w_start, w_end, seq_lengths[seq], seq_gaps)

        processed_windows.append({
            'Sequence': seq,
            'Start': w_start,
            'End': w_end,
            'high': row['high'],
            'low': row['low'],
            'MinDist': min_dist,
            'EffectiveLen': effective_len,
            'High_per_kbp': (row['high'] / effective_len) * 1000,
            'Low_per_kbp': (row['low'] / effective_len) * 1000,
            'IsEnd': min_dist <= args.end_threshold
        })

    win_df = pd.DataFrame(processed_windows)

    if len(win_df) < 1:
        logger.error("No valid windows found after processing.")
        return

    # Global Statistics
    total_high = win_df['high'].sum()
    total_low = win_df['low'].sum()
    global_asymmetry = (total_high - total_low) / (total_high + total_low) if (total_high + total_low) > 0 else 0

    if len(win_df) > 1:
        pearson_corr, _ = pearsonr(win_df['low'], win_df['high'])
    else:
        pearson_corr = np.nan

    jaccard = calculate_jaccard(win_df, args.jaccard_minimum)

    # Global End-Associated Analysis
    end_wins = win_df[win_df['IsEnd']]
    bg_wins = win_df[~win_df['IsEnd']]

    global_end_stats = {
        'High_Density_End': end_wins['High_per_kbp'].mean() if not end_wins.empty else 0,
        'Low_Density_End': end_wins['Low_per_kbp'].mean() if not end_wins.empty else 0,
        'High_Density_Background': bg_wins['High_per_kbp'].mean() if not bg_wins.empty else 0,
        'Low_Density_Background': bg_wins['Low_per_kbp'].mean() if not bg_wins.empty else 0
    }

    # Sequence-level stats
    seq_stats = []
    for seq, group in win_df.groupby('Sequence'):
        s_high = group['high'].sum()
        s_low = group['low'].sum()
        s_asym = (s_high - s_low) / (s_high + s_low) if (s_high + s_low) > 0 else 0
        s_len = group['End'].max()

        s_end = group[group['IsEnd']]
        s_bg = group[~group['IsEnd']]

        s_high_end = s_end['High_per_kbp'].mean() if not s_end.empty else 0
        s_high_bg = s_bg['High_per_kbp'].mean() if not s_bg.empty else 0
        s_low_end = s_end['Low_per_kbp'].mean() if not s_end.empty else 0
        s_low_bg = s_bg['Low_per_kbp'].mean() if not s_bg.empty else 0

        seq_stats.append({
            'Sequence': seq,
            'Length': s_len,
            'Total_High': s_high,
            'Total_Low': s_low,
            'Asymmetry_Index': s_asym,
            'High_Density_End': s_high_end,
            'High_Density_Background': s_high_bg,
            'High_End_Ratio': s_high_end / s_high_bg if s_high_bg > 0 else np.nan,
            'High_End_Diff': s_high_end - s_high_bg,
            'Low_Density_End': s_low_end,
            'Low_Density_Background': s_low_bg,
            'Low_End_Ratio': s_low_end / s_low_bg if s_low_bg > 0 else np.nan,
            'Low_End_Diff': s_low_end - s_low_bg
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

        f.write(f"Global_Mean_High_Density_End(dist<={args.end_threshold})\t{global_end_stats['High_Density_End']:.4f}\n")
        f.write(f"Global_Mean_High_Density_Background\t{global_end_stats['High_Density_Background']:.4f}\n")
        f.write(f"Global_High_End_Background_Diff\t{global_end_stats['High_Density_End'] - global_end_stats['High_Density_Background']:.4f}\n")
        h_ratio = global_end_stats['High_Density_End'] / global_end_stats['High_Density_Background'] if global_end_stats['High_Density_Background'] > 0 else 0
        f.write(f"Global_High_End_Background_Ratio\t{h_ratio:.4f}\n")

        f.write(f"Global_Mean_Low_Density_End(dist<={args.end_threshold})\t{global_end_stats['Low_Density_End']:.4f}\n")
        f.write(f"Global_Mean_Low_Density_Background\t{global_end_stats['Low_Density_Background']:.4f}\n")
        f.write(f"Global_Low_End_Background_Diff\t{global_end_stats['Low_Density_End'] - global_end_stats['Low_Density_Background']:.4f}\n")
        l_ratio = global_end_stats['Low_Density_End'] / global_end_stats['Low_Density_Background'] if global_end_stats['Low_Density_Background'] > 0 else 0
        f.write(f"Global_Low_End_Background_Ratio\t{l_ratio:.4f}\n")

        if args.ngaps:
            f.write(f"NGaps_GFF\t{args.ngaps}\n")
        f.write("\n# Per-Sequence Statistics\n")

    df_seq_stats.to_csv(stats_file, sep='\t', index=False, mode='a')

    logger.info(f"Statistics written to {stats_file}")

    # Plotting
    plot_scatter(win_df, args.output)
    plot_asymmetry_v_length(df_seq_stats, args.output)
    plot_density_v_distance(win_df, args.output)
    plot_end_comparison(global_end_stats, args.output)

    logger.info(f"Plots generated with prefix {args.output}")

if __name__ == "__main__":
    main()
