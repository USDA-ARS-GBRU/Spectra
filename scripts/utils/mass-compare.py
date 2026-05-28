#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import os
import logging
from matplotlib.colors import LinearSegmentedColormap

# Logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger()

def calculate_jaccard(df, threshold):
    high_present = df['high'] > threshold
    low_present = df['low'] > threshold
    intersection = np.logical_and(high_present, low_present).sum()
    union = np.logical_or(high_present, low_present).sum()
    return intersection / union if union > 0 else 0.0

def plot_scatter(win_df, output_prefix, end_threshold):
    fig, ax = plt.subplots(figsize=(10, 10))

    # Background windows (non-outliers)
    bg = win_df[~win_df['IsOutlier']]
    ax.scatter(bg['low'], bg['high'], alpha=0.3, s=10, color='blue', label='Background')

    # Outliers
    outliers = win_df[win_df['IsOutlier']]

    if not outliers.empty:
        if end_threshold > 0:
            dists = np.clip(outliers['MinDist'], 0, end_threshold)
            norm_dists = dists / end_threshold
            colors = [(1, 0, 0), (1, 0.75, 0.8)] # Red to Pink
            cm = LinearSegmentedColormap.from_list('outlier_cm', colors, N=100)
            sc = ax.scatter(outliers['low'], outliers['high'], c=norm_dists, cmap=cm, s=15, alpha=0.8, label='Outliers', vmin=0, vmax=1)
            cbar = plt.colorbar(sc, ax=ax)
            cbar.set_label(f'Distance from Feature (0 to {end_threshold})')
        else:
            ax.scatter(outliers['low'], outliers['high'], alpha=0.8, s=15, color='red', label='Outliers (Dist 0)')

    # Add diagonal line
    max_val = max(win_df['high'].max(), win_df['low'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)

    ax.set_xlabel('Low Extreme Kmer Count')
    ax.set_ylabel('High Extreme Kmer Count')
    ax.set_title('Scatter Plot: Extreme Kmer Counts (Outliers Highlighted)')
    ax.legend()

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_scatter.png", dpi=300)
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
    min_dist = w_start - 1
    f_type = "End-Adjacent"
    f_id = "Left-End"
    r_dist = length - w_end
    if r_dist < min_dist:
        min_dist = r_dist
        f_id = "Right-End"
    for idx, (g_start, g_end) in enumerate(seq_gaps):
        if g_end < w_start:
            d = w_start - g_end - 1
        elif g_start > w_end:
            d = g_start - w_end - 1
        else:
            d = 0
        if d < min_dist:
            min_dist = d
            f_type = "N-Adjacent"
            f_id = f"Gap-{idx}"
            if min_dist == 0:
                break
    return max(0, min_dist), f_type, f_id

def main():
    parser = argparse.ArgumentParser(description="Kmer Mass Compare: Compare extreme kmer accumulations across genome windows")
    parser.add_argument('-i', '--input', required=True, help='Input TSV from mass-query.py')
    parser.add_argument('-o', '--output', default='mass_compare', help='Output prefix for plots and stats')
    parser.add_argument('--jaccard-minimum', type=int, default=0, help='Minimum count threshold for Jaccard coincidence [default 0]')
    parser.add_argument('--ngaps', help='GFF file of N-gap coordinates')
    parser.add_argument('--end-threshold', type=int, default=0, help='Distance threshold for filtering outliers and coloring scatter [default 0]')
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
            'Sequence': seq, 'Start': w_start, 'End': w_end,
            'high': row['high'], 'low': row['low'],
            'MinDist': min_dist, 'FeatureType': f_type, 'FeatureID': f"{seq}:{f_id}",
            'EffectiveLen': effective_len,
            'High_per_kbp': (row['high'] / effective_len) * 1000,
            'Low_per_kbp': (row['low'] / effective_len) * 1000
        })

    win_df = pd.DataFrame(processed_windows)
    if len(win_df) < 1:
        logger.error("No valid windows found.")
        return

    # Background Statistics
    bg_high_mean = win_df['High_per_kbp'].mean()
    bg_high_sd = win_df['High_per_kbp'].std()
    bg_low_mean = win_df['Low_per_kbp'].mean()
    bg_low_sd = win_df['Low_per_kbp'].std()

    # Per-sequence background stats
    seq_backgrounds = {}
    for seq, group in win_df.groupby('Sequence'):
        seq_backgrounds[seq] = {
            'high_mean': group['High_per_kbp'].mean(),
            'high_sd': group['High_per_kbp'].std(),
            'low_mean': group['Low_per_kbp'].mean(),
            'low_sd': group['Low_per_kbp'].std()
        }

    # Outlier flagging
    outlier_rows = []
    win_df['IsOutlier'] = False

    for idx, row in win_df.iterrows():
        seq = row['Sequence']
        s_bg = seq_backgrounds[seq]

        is_global_high = row['High_per_kbp'] > (bg_high_mean + bg_high_sd)
        is_seq_high = row['High_per_kbp'] > (s_bg['high_mean'] + s_bg['high_sd'])
        is_global_low = row['Low_per_kbp'] > (bg_low_mean + bg_low_sd)
        is_seq_low = row['Low_per_kbp'] > (s_bg['low_mean'] + s_bg['low_sd'])

        if is_global_high or is_seq_high or is_global_low or is_seq_low:
            if row['MinDist'] <= args.end_threshold:
                win_df.at[idx, 'IsOutlier'] = True
                out_types = []
                if is_global_high: out_types.append("Global-High")
                if is_seq_high: out_types.append("Sequence-High")
                if is_global_low: out_types.append("Global-Low")
                if is_seq_low: out_types.append("Sequence-Low")

                outlier_rows.append({
                    'Sequence': seq, 'Length': seq_lengths[seq], 'Window': f"{row['Start']}-{row['End']}",
                    'High_Density': row['High_per_kbp'], 'Low_Density': row['Low_per_kbp'],
                    'MinDist': row['MinDist'], 'FeatureType': row['FeatureType'],
                    'FeatureID': row['FeatureID'], 'OutlierType': ",".join(out_types)
                })

    df_outliers = pd.DataFrame(outlier_rows)

    # Sequence-level Summary Stats
    seq_stats_list = []
    for seq, group in win_df.groupby('Sequence'):
        s_high = group['high'].sum()
        s_low = group['low'].sum()
        s_asym = (s_high - s_low) / (s_high + s_low) if (s_high + s_low) > 0 else 0

        s_end = group[group['MinDist'] == 0]
        s_bg = group[group['MinDist'] > 0]

        s_high_end = s_end['High_per_kbp'].mean() if not s_end.empty else 0
        s_high_bg = s_bg['High_per_kbp'].mean() if not s_bg.empty else group['High_per_kbp'].mean()
        s_low_end = s_end['Low_per_kbp'].mean() if not s_end.empty else 0
        s_low_bg = s_bg['Low_per_kbp'].mean() if not s_bg.empty else group['Low_per_kbp'].mean()

        seq_stats_list.append({
            'Sequence': seq, 'Length': seq_lengths[seq], 'Total_High': s_high, 'Total_Low': s_low,
            'Asymmetry_Index': s_asym, 'High_Density_End': s_high_end, 'High_Density_Background': s_high_bg,
            'Low_Density_End': s_low_end, 'Low_Density_Background': s_low_bg
        })

    df_seq_stats = pd.DataFrame(seq_stats_list)

    total_high = win_df['high'].sum()
    total_low = win_df['low'].sum()

    # Output files
    stats_file = f"{args.output}.stats"
    with open(stats_file, 'w') as f:
        f.write("# Global Statistics\n")
        f.write(f"Global_Asymmetry_Index\t{(total_high - total_low) / (total_high + total_low) if (total_high + total_low) > 0 else 0:.4f}\n")
        f.write(f"Pearson_Correlation\t{pearsonr(win_df['low'], win_df['high'])[0] if len(win_df)>1 else 0:.4f}\n")
        f.write(f"Jaccard_Coincidence_Index(>{args.jaccard_minimum})\t{calculate_jaccard(win_df, args.jaccard_minimum):.4f}\n")
        f.write(f"Global_Background_High_Mean\t{bg_high_mean:.4f}\n")
        f.write(f"Global_Background_High_SD\t{bg_high_sd:.4f}\n")
        f.write(f"Global_Background_Low_Mean\t{bg_low_mean:.4f}\n")
        f.write(f"Global_Background_Low_SD\t{bg_low_sd:.4f}\n")
        f.write("\n# Per-Sequence Statistics\n")
    df_seq_stats.to_csv(stats_file, sep='\t', index=False, mode='a')

    outliers_file = f"{args.output}.outliers.tsv"
    df_outliers.to_csv(outliers_file, sep='\t', index=False)

    logger.info(f"Statistics written to {stats_file}")
    logger.info(f"Outliers written to {outliers_file}")

    # Plotting
    plot_scatter(win_df, args.output, args.end_threshold)

    global_end_stats = {
        'High_Density_End': win_df[win_df['MinDist']==0]['High_per_kbp'].mean(),
        'Low_Density_End': win_df[win_df['MinDist']==0]['Low_per_kbp'].mean(),
        'High_Density_Background': win_df[win_df['MinDist']>0]['High_per_kbp'].mean(),
        'Low_Density_Background': win_df[win_df['MinDist']>0]['Low_per_kbp'].mean()
    }
    plot_end_comparison(global_end_stats, args.output)

    logger.info(f"Plots generated with prefix {args.output}")

if __name__ == "__main__":
    main()
