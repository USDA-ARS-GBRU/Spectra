#!/usr/bin/env python3
import argparse
import math
import random
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import logging
from scipy.stats import gaussian_kde

# CLI arguments
parser = argparse.ArgumentParser(description="Streamed comparison of raw vs assembly k-mer dump files")
parser.add_argument("-r", "--raw_dump", type=str, required=True, help="Raw k-mer dump file")
parser.add_argument("-a", "--asm_dump", type=str, required=True, help="Assembly k-mer dump file")
parser.add_argument("-k", "--kmer_size", type=int, default=20, help="K-mer size [default 20]")
parser.add_argument("-o", "--output_prefix", type=str, default="kmerComp_output", help="Output prefix")
parser.add_argument("-f", "--output_format", type=str, default="png", help="Output image format [default png]")
parser.add_argument("-s", "--plot_sample", type=int, default=1000000, help="Number of kmers to sample for plots")
parser.add_argument("--percentile-low", dest='percentile_low', type=float, default=1, help="Bottom N percent of kmers for extreme scatter plot [default 1]")
parser.add_argument("--percentile-high", dest='percentile_high', type=float, default=1, help="Top N percent of kmers for extreme scatter plot [default 1]")
parser.add_argument("--auto", action='store_true', help='Automatically determine low/high percentiles based on distribution')
parser.add_argument("-p", "--percentile", type=float, default=None, help="Deprecated: use --percentile-low and --percentile-high instead")
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)

args = parser.parse_args()

# Logging
logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')
logger = logging.getLogger()
if args.verbose:
    logger.setLevel(logging.INFO)


# Streaming merge
def stream_merge(raw_file, asm_file, sample_size):
    """
    Stream through two sorted dump files, yield matched kmers, and keep a reservoir sample.
    """
    sample = []
    total = 0

    try:
        with open(raw_file) as fr, open(asm_file) as fa:
            line_r = fr.readline().split()
            line_a = fa.readline().split()

            while line_r and line_a:
                if line_r[0] == line_a[0]:
                    kmer, count_r, count_a = line_r[0], int(line_r[1]), int(line_a[1])
                    total += 1

                    # Reservoir sampling
                    if len(sample) < sample_size:
                        sample.append((kmer, count_r, count_a))
                    else:
                        j = random.randint(0, total - 1)
                        if j < sample_size:
                            sample[j] = (kmer, count_r, count_a)

                    line_r = fr.readline().split()
                    line_a = fa.readline().split()

                elif line_r[0] < line_a[0]:
                    line_r = fr.readline().split()
                else:
                    line_a = fa.readline().split()
    except Exception as e:
        logger.error(f"Error during stream merge: {e}")
        return [], 0

    return sample, total

logger.info("Streaming merge and sampling...")
sample_kmers, total_pairs = stream_merge(args.raw_dump, args.asm_dump, args.plot_sample)

if not sample_kmers:
    logger.error("No matched kmers found or error occurred.")
    exit(1)

logger.info(f"Processed {total_pairs:,} matched kmers; kept {len(sample_kmers):,} in sample.")


# Transformations on sample
df = pd.DataFrame(sample_kmers, columns=["kmer", "RawCount", "AsmCount"])
df["logRaw"] = np.log10(df["RawCount"] + 1)
df["logAsm"] = np.log10(df["AsmCount"] + 1)
df["reduction"] = df["logAsm"] - df["logRaw"]
df["reductionRank"] = df["reduction"].rank(method="first")

if args.percentile is not None:
    args.percentile_low = args.percentile
    args.percentile_high = args.percentile

if args.auto:
    mu = df["reduction"].mean()
    sigma = df["reduction"].std()
    low_thresh = mu - 2 * sigma
    high_thresh = mu + 2 * sigma

    df_sorted = df.sort_values("reduction")
    low_cut_idx = df_sorted["reduction"].searchsorted(low_thresh, side='right')
    high_cut_idx = df_sorted["reduction"].searchsorted(high_thresh, side='left')

    args.percentile_low = (low_cut_idx / len(df)) * 100
    args.percentile_high = ((len(df) - high_cut_idx) / len(df)) * 100

    df_extreme = pd.concat([df_sorted.iloc[:low_cut_idx], df_sorted.iloc[high_cut_idx:]])
    logger.info(f"Auto-detected thresholds for sample: low={args.percentile_low:.2f}%, high={args.percentile_high:.2f}%")

    # Persist auto-calculated percentiles
    percentiles_file = f"{args.output_prefix}_percentiles.txt"
    try:
        with open(percentiles_file, 'w') as f:
            f.write(f"PERCENTILE_LOW={args.percentile_low:.4f}\n")
            f.write(f"PERCENTILE_HIGH={args.percentile_high:.4f}\n")
        logger.info(f"Saved auto-percentiles to {percentiles_file}")
    except Exception as e:
        logger.error(f"Failed to save percentiles file: {e}")
else:
    # Percentile filtering on sample
    n_rows = len(df)
    low_cut = int(math.ceil(n_rows * (args.percentile_low / 100.0)))
    high_cut = int(math.floor(n_rows * (1 - args.percentile_high / 100.0)))
    df_sorted = df.sort_values("reductionRank")
    df_extreme = pd.concat([df_sorted.iloc[:low_cut], df_sorted.iloc[high_cut:]])

xmin = df_sorted["RawCount"].min()
xmax = df_sorted["RawCount"].max()
ymin = df_sorted["AsmCount"].min()
ymax = df_sorted["AsmCount"].max()

# Plotting functions
sns.set(style="whitegrid")

# Scatter
plt.figure(figsize=(8, 8))

# Composite color mapping:
# Normalized log counts for color mapping
r_min, r_max = df["logRaw"].min(), df["logRaw"].max()
a_min, a_max = df["logAsm"].min(), df["logAsm"].max()


def normalize_color(vals, vmin, vmax):
    if vmin == vmax:
        return np.zeros_like(vals)
    return (vals - vmin) / (vmax - vmin)


norm_r = normalize_color(df["logRaw"].values, r_min, r_max)
norm_a = normalize_color(df["logAsm"].values, a_min, a_max)

# Color mapping: Bilinear interpolation between four corners
c00 = np.array([0.9, 0.8, 0.9])
c10 = np.array([0.5, 0.0, 0.0])
c01 = np.array([0.0, 0.0, 0.5])
c11 = np.array([0.25, 0.0, 0.25])

# Vectorized bilinear interpolation
colors = (np.outer((1 - norm_r) * (1 - norm_a), c00) +
          np.outer(norm_r * (1 - norm_a), c10) +
          np.outer((1 - norm_r) * norm_a, c01) +
          np.outer(norm_r * norm_a, c11))

plt.scatter(df["RawCount"], df["AsmCount"], s=1, alpha=0.3, c=colors, edgecolors='none')

plt.xscale("log")
plt.yscale("log")
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.xlabel("Kmers in raw data")
plt.ylabel("Kmers in assembly")
plt.title(f"K={args.kmer_size} coverage (scatter)")
plt.savefig(f"{args.output_prefix}_k{args.kmer_size}_scatter.{args.output_format}", dpi=200)
plt.close()

# Scatter (extremes only)
if not df_extreme.empty:
    plt.figure(figsize=(8, 8))
    plt.scatter(df_extreme["RawCount"], df_extreme["AsmCount"],
                c=df_extreme["reduction"], cmap="coolwarm", s=2, alpha=0.6)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel("Kmers in raw data")
    plt.ylabel("Kmers in assembly")
    plt.title(f"K={args.kmer_size} extreme kmers (Low {args.percentile_low}%, High {args.percentile_high}%)")
    plt.colorbar(label="Log-fold change")
    # For compatibility with pdfReport, we keep a generic name if they are equal, or use a new naming scheme
    suffix = f"{args.percentile_low}pct" if args.percentile_low == args.percentile_high else f"L{args.percentile_low}_H{args.percentile_high}pct"
    plt.savefig(f"{args.output_prefix}_k{args.kmer_size}_scatter_extreme_{suffix}.{args.output_format}", dpi=200)
    plt.close()

# ECDF
plt.figure(figsize=(8, 6))
ecdf_x = np.sort(df["reduction"])
ecdf_y = np.arange(1, len(ecdf_x) + 1) / len(ecdf_x)
plt.step(ecdf_x, ecdf_y, where="post", color="red")
plt.xlabel("Log-fold kmer change")
plt.ylabel("Cumulative probability")
plt.title(f"K={args.kmer_size} empirical cumulative distribution")
plt.savefig(f"{args.output_prefix}_k{args.kmer_size}_ecdf.{args.output_format}", dpi=200)
plt.close()

# Density
plt.figure(figsize=(8, 6))
sns.kdeplot(df["reduction"], fill=True, alpha=0.6, color="#399602")
plt.xlabel("Log-fold kmer change")
plt.ylabel("Density")
plt.title(f"K={args.kmer_size} coverage density")
plt.savefig(f"{args.output_prefix}_k{args.kmer_size}_density.{args.output_format}", dpi=200)
plt.close()

# Violin
plt.figure(figsize=(6, 6))
sns.violinplot(data=[df["logRaw"], df["logAsm"]], palette=["#1f77b4", "#ff7f0e"])
plt.xticks([0, 1], ["Raw", "Asm"])
plt.ylabel("Abundance (log10)")
plt.title(f"K={args.kmer_size} coverage abundance")
plt.savefig(f"{args.output_prefix}_k{args.kmer_size}_violin.{args.output_format}", dpi=200)
plt.close()

# Back-to-back density plot
plt.figure(figsize=(6, 6))

y_vals = np.linspace(
    min(df["logRaw"].min(), df["logAsm"].min()),
    max(df["logRaw"].max(), df["logAsm"].max()),
    500
)

raw_kde = gaussian_kde(df["logRaw"], bw_method=0.5)(y_vals)
asm_kde = gaussian_kde(df["logAsm"], bw_method=0.5)(y_vals)

# Mirror plot: Raw to the left (negative x), Asm to the right
plt.fill_betweenx(y_vals, -raw_kde, 0, color="#1f77b4", alpha=0.6, label="Raw")
plt.fill_betweenx(y_vals, 0, asm_kde, color="#ff7f0e", alpha=0.6, label="Asm")

plt.axvline(0, color="black", linewidth=0.8)
plt.xlabel("Density (mirrored)")
plt.ylabel("Abundance (log10)")
plt.title(f"K={args.kmer_size} back-to-back density")
plt.legend()
plt.savefig(f"{args.output_prefix}_k{args.kmer_size}_back2back_density.{args.output_format}", dpi=200)
plt.close()

logger.info("Done.")
