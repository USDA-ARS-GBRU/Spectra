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

def calculate_nearest_feature(w_start, w_end, length, seq_gaps):
    # Features: Left End, Right End, or N-gap
    # min_dist, type, feature_id

    # Left End
    min_dist = w_start - 1
    f_type = "End-Adjacent"
    f_id = "Left-End"

    # Right End
    r_dist = length - w_end
    if r_dist < min_dist:
        min_dist = r_dist
        f_id = "Right-End"

    # Gaps
    for idx, (g_start, g_end) in enumerate(seq_gaps):
        # Distance from window to gap
        if g_end < w_start:
            d = w_start - g_end - 1
        elif g_start > w_end:
            d = g_start - w_end - 1
        else:
            d = 0 # Overlaps gap

        if d < min_dist:
            min_dist = d
            f_type = "N-Adjacent"
            f_id = f"Gap-{idx}"
            if min_dist == 0:
                break

    return max(0, min_dist), f_type, f_id

def plot_density_v_distance(plot_data, output_prefix):
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.scatter(plot_data['MinDist'], plot_data['High_per_kbp'], alpha=0.4, color='red', s=5, label='High Bin')
    ax.scatter(plot_data['MinDist'], plot_data['Low_per_kbp'], alpha=0.4, color='blue', s=5, label='Low Bin')
    ax.set_xlabel('Distance from Sequence End or N-gap (bp)')
    ax.set_ylabel('Extreme Kmer Density (counts per kbp)')
    ax.set_title('Extreme Kmer Density vs Distance from Nearest Feature')
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

def plot_ranked_features(feature_stats, bin_type, output_prefix):
    # feature_stats: list of {id, density, type}
    df = pd.DataFrame(feature_stats).sort_values('density', ascending=False).head(30)
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(12, 8))
    colors = ['red' if t == 'End-Adjacent' else 'orange' for t in df['type']]
    ax.bar(df['id'], df['density'], color=colors)
    ax.set_ylabel(f'Mean {bin_type} Density (counts per kbp)')
    ax.set_xlabel('Feature Instance')
    ax.set_title(f'Top 30 Ranked Features by {bin_type} Density')
    plt.xticks(rotation=90, fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_ranked_{bin_type.lower()}.png", dpi=300)
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

    df = pd.read_csv(args.input, sep='\t')
    pivot_df = df.pivot_table(index=['Sequence', 'Start', 'End'], columns='Bin', values='Count', fill_value=0).reset_index()
    if 'high' not in pivot_df.columns: pivot_df['high'] = 0
    if 'low' not in pivot_df.columns: pivot_df['low'] = 0

    gaps = parse_ngaps(args.ngaps) if args.ngaps else {}
    seq_lengths = pivot_df.groupby('Sequence')['End'].max().to_dict()

    processed_windows = []
    for idx, row in pivot_df.iterrows():
        seq = row['Sequence']
        w_start = row['Start']
        w_end = row['End']
        seq_gaps = gaps.get(seq, [])
        gap_len = get_gap_overlap(w_start, w_end, seq_gaps)
        w_len = w_end - w_start + 1
        effective_len = w_len - gap_len
        if effective_len <= 0: continue

        min_dist, f_type, f_id = calculate_nearest_feature(w_start, w_end, seq_lengths[seq], seq_gaps)

        processed_windows.append({
            'Sequence': seq,
            'Start': w_start,
            'End': w_end,
            'high': row['high'],
            'low': row['low'],
            'MinDist': min_dist,
            'FeatureType': f_type,
            'FeatureID': f"{seq}:{f_id}",
            'EffectiveLen': effective_len,
            'High_per_kbp': (row['high'] / effective_len) * 1000,
            'Low_per_kbp': (row['low'] / effective_len) * 1000,
            'IsEnd': min_dist <= args.end_threshold
        })

    win_df = pd.DataFrame(processed_windows)
    if len(win_df) < 1:
        logger.error("No valid windows found.")
        return

    # Statistics
    total_high = win_df['high'].sum()
    total_low = win_df['low'].sum()
    global_asymmetry = (total_high - total_low) / (total_high + total_low) if (total_high + total_low) > 0 else 0
    if len(win_df) > 1:
        pearson_corr, _ = pearsonr(win_df['low'], win_df['high'])
    else:
        pearson_corr = np.nan
    jaccard = calculate_jaccard(win_df, args.jaccard_minimum)

    # End analysis
    end_wins = win_df[win_df['IsEnd']]
    bg_wins = win_df[~win_df['IsEnd']]

    global_bg_high_mean = bg_wins['High_per_kbp'].mean() if not bg_wins.empty else win_df['High_per_kbp'].mean()
    global_bg_high_sd = bg_wins['High_per_kbp'].std() if not bg_wins.empty else win_df['High_per_kbp'].std()
    global_bg_low_mean = bg_wins['Low_per_kbp'].mean() if not bg_wins.empty else win_df['Low_per_kbp'].mean()
    global_bg_low_sd = bg_wins['Low_per_kbp'].std() if not bg_wins.empty else win_df['Low_per_kbp'].std()

    global_end_stats = {
        'High_Density_End': end_wins['High_per_kbp'].mean() if not end_wins.empty else 0,
        'Low_Density_End': end_wins['Low_per_kbp'].mean() if not end_wins.empty else 0,
        'High_Density_Background': global_bg_high_mean,
        'Low_Density_Background': global_bg_low_mean
    }

    # Outliers detection
    outliers = []

    # Sequence-level stats
    seq_stats_list = []
    for seq, group in win_df.groupby('Sequence'):
        s_len = seq_lengths[seq]
        s_high = group['high'].sum()
        s_low = group['low'].sum()
        s_asym = (s_high - s_low) / (s_high + s_low) if (s_high + s_low) > 0 else 0
        s_end = group[group['IsEnd']]
        s_bg = group[~group['IsEnd']]

        s_high_end = s_end['High_per_kbp'].mean() if not s_end.empty else 0
        s_high_bg_mean = s_bg['High_per_kbp'].mean() if not s_bg.empty else group['High_per_kbp'].mean()
        s_high_bg_sd = s_bg['High_per_kbp'].std() if not s_bg.empty else group['High_per_kbp'].std()

        s_low_end = s_end['Low_per_kbp'].mean() if not s_end.empty else 0
        s_low_bg_mean = s_bg['Low_per_kbp'].mean() if not s_bg.empty else group['Low_per_kbp'].mean()
        s_low_bg_sd = s_bg['Low_per_kbp'].std() if not s_bg.empty else group['Low_per_kbp'].std()

        seq_stats_list.append({
            'Sequence': seq, 'Length': s_len, 'Total_High': s_high, 'Total_Low': s_low,
            'Asymmetry_Index': s_asym, 'High_Density_End': s_high_end, 'High_Density_Background': s_high_bg_mean,
            'Low_Density_End': s_low_end, 'Low_Density_Background': s_low_bg_mean
        })

        # Detect outliers for this sequence
        for _, w in group.iterrows():
            is_global_high = w['High_per_kbp'] > (global_bg_high_mean + global_bg_high_sd)
            is_seq_high = w['High_per_kbp'] > (s_high_bg_mean + s_high_bg_sd)
            is_global_low = w['Low_per_kbp'] > (global_bg_low_mean + global_bg_low_sd)
            is_seq_low = w['Low_per_kbp'] > (s_low_bg_mean + s_low_bg_sd)

            if is_global_high or is_seq_high or is_global_low or is_seq_low:
                out_type = []
                if is_global_high: out_type.append("Global-High")
                if is_seq_high: out_type.append("Sequence-High")
                if is_global_low: out_type.append("Global-Low")
                if is_seq_low: out_type.append("Sequence-Low")

                outliers.append({
                    'Sequence': seq, 'Length': s_len, 'Window': f"{w['Start']}-{w['End']}",
                    'High_Density': w['High_per_kbp'], 'Low_Density': w['Low_per_kbp'],
                    'MinDist': w['MinDist'], 'FeatureType': w['FeatureType'],
                    'FeatureID': w['FeatureID'], 'OutlierType': ",".join(out_type)
                })

    df_seq_stats = pd.DataFrame(seq_stats_list)
    df_outliers = pd.DataFrame(outliers)

    # Feature Stats for Ranking
    high_feature_stats = []
    low_feature_stats = []
    for f_id, group in win_df.groupby('FeatureID'):
        high_feature_stats.append({'id': f_id, 'density': group['High_per_kbp'].mean(), 'type': group['FeatureType'].iloc[0]})
        low_feature_stats.append({'id': f_id, 'density': group['Low_per_kbp'].mean(), 'type': group['FeatureType'].iloc[0]})

    # Output
    stats_file = f"{args.output}.stats"
    with open(stats_file, 'w') as f:
        f.write("# Global Statistics\n")
        f.write(f"Global_Total_High\t{total_high}\n")
        f.write(f"Global_Total_Low\t{total_low}\n")
        f.write(f"Global_Asymmetry_Index\t{global_asymmetry:.4f}\n")
        f.write(f"Pearson_Correlation\t{pearson_corr:.4f}\n")
        f.write(f"Jaccard_Coincidence_Index(>{args.jaccard_minimum})\t{jaccard:.4f}\n")
        f.write(f"Global_Background_High_Mean\t{global_bg_high_mean:.4f}\n")
        f.write(f"Global_Background_High_SD\t{global_bg_high_sd:.4f}\n")
        f.write(f"Global_Background_Low_Mean\t{global_bg_low_mean:.4f}\n")
        f.write(f"Global_Background_Low_SD\t{global_bg_low_sd:.4f}\n")
        f.write("\n# Per-Sequence Statistics\n")
    df_seq_stats.to_csv(stats_file, sep='\t', index=False, mode='a')

    outliers_file = f"{args.output}.outliers.tsv"
    df_outliers.to_csv(outliers_file, sep='\t', index=False)

    logger.info(f"Statistics written to {stats_file}")
    logger.info(f"Outliers written to {outliers_file}")

    # Plotting
    plot_scatter(win_df, args.output)
    plot_asymmetry_v_length(df_seq_stats, args.output)
    plot_density_v_distance(win_df, args.output)
    plot_end_comparison(global_end_stats, args.output)
    plot_ranked_features(high_feature_stats, "High", args.output)
    plot_ranked_features(low_feature_stats, "Low", args.output)

    logger.info(f"Plots generated with prefix {args.output}")

if __name__ == "__main__":
    main()
