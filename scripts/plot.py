#!/usr/bin/env python3

import os
import pandas as pd
import logging
import spectral
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import matplotlib.ticker as ticker
import plotly.express as px

logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')
logger = logging.getLogger()

# Consistent plot parameters for alignment
PLOT_CONFIG = {
    'width': 10,
    'height_spectra': 2.5,
    'height_mass': 2.0,
    'height_gff': 0.4,
    'dpi': 300,
    'left_margin': 0.1,   # 10% of width
    'right_margin': 0.95, # 95% of width
    'bottom_margin': 0.2, # 20% of height
    'top_margin': 0.9     # 90% of height
}

def plotly_plot(df_melted, query_colors, queries):
    fig = px.area(df_melted, x="Start", y="Count", color="Triplet", color_discrete_map=query_colors,
                  line_group="Triplet", category_orders={"Triplet": list(reversed(queries))})
    fig.for_each_trace(lambda trace: trace.update(fillcolor=trace.line.color))
    fig.update_layout(margin_b=0, margin_t=0, margin_l=0, margin_r=0, xaxis=dict(rangeslider=dict(visible=True)))
    return fig

def palette_builder(triplet, palette="base", bases=None):
    if bases is None:
        bases = ["A", "C", "G", "T"]
    if palette == "base":
        colors = ["C6", "6C", "3C", "10"]
        try:
            return f"#{colors[bases.index(triplet[0])].lower()}{colors[bases.index(triplet[1])].lower()}{colors[bases.index(triplet[2])].lower()}"
        except (IndexError, ValueError):
            return "#000000"
    return "#000000"

def get_triplet_colors(queries, palette="base"):
    script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    palette_csv = os.path.join(script_dir, "includes", "paletteMatrix_base.csv")

    if os.path.exists(palette_csv):
        palette_order_df = pd.read_csv(palette_csv, header=None)
        palette_order = palette_order_df.values.flatten().tolist()
        ordered_queries = [q for q in palette_order if q in queries]
        # Add any queries not in the palette CSV at the end
        ordered_queries += [q for q in queries if q not in ordered_queries]
    else:
        ordered_queries = sorted(queries)

    colors = {q: palette_builder(q, palette=palette) for q in ordered_queries}
    return ordered_queries, colors

def setup_axes(ax, x_min, x_max, show_axes=True):
    ax.set_xlim(x_min, x_max)
    if not show_axes:
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_xlabel('')
        ax.spines['bottom'].set_visible(False)
    else:
        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.set_xlabel("Window Position (nucleotide)", fontsize=8)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def plot_spectra(df, output_path, sequence, x_min, x_max, show_axes=True, frequencies=False, ngaps_path=None, y_max=None):
    queries = [c for c in df.columns if len(c) == 3 and all(b in 'ACGT' for b in c)]
    ordered_queries, color_map = get_triplet_colors(queries)
    # Queries sorted to match the expected color gradient
    ordered_queries.sort(reverse=True)
    df = df.sort_values('Start')

    if not frequencies:
        window_sizes = df['End'] - df['Start'] + 1
        plot_df = df[ordered_queries].div(window_sizes, axis=0)
    else:
        plot_df = df[ordered_queries]

    x = (df['Start'] + df['End']) / 2
    y_values = plot_df[ordered_queries].values
    y_stack = np.cumsum(y_values, axis=1)

    fig, ax = plt.subplots(figsize=(PLOT_CONFIG['width'], PLOT_CONFIG['height_spectra']))
    plt.subplots_adjust(left=PLOT_CONFIG['left_margin'], right=PLOT_CONFIG['right_margin'],
                        bottom=PLOT_CONFIG['bottom_margin'], top=PLOT_CONFIG['top_margin'])

    y_prev = np.zeros(len(df))
    for i, q in enumerate(ordered_queries):
        # Looking to transition away from cum-sum plots for this aspect. Matplotlib in python ultimately is the answer for deprecating R, but the linewidth hackaround is not healthy
        ax.fill_between(x, y_prev, y_stack[:, i], color=color_map[q], step='mid', linewidth=.15, edgecolor=color_map[q])
        y_prev = y_stack[:, i]

    if ngaps_path and os.path.exists(ngaps_path):
        try:
            ngaps = pd.read_csv(ngaps_path, sep='\t', comment='#', header=None,
                                names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
            ngaps = ngaps[ngaps['seqid'] == sequence]
            for _, row in ngaps.iterrows():
                ax.axvspan(row['start'], row['end'], color='black', alpha=0.5, zorder=10)
        except Exception as e:
            logger.error(f"Error reading ngaps GFF {ngaps_path}: {e}")

    setup_axes(ax, x_min, x_max, show_axes)
    if y_max:
        ax.set_ylim(0, y_max)
    else:
        ax.set_ylim(0, 1.0)
    ax.set_ylabel('Proportion' if not frequencies else 'Frequency', fontsize=8)

    plt.savefig(output_path, dpi=PLOT_CONFIG['dpi'])
    plt.close()

def plot_mass(df, output_prefix, sequence, x_min, x_max, show_axes=True, y_max=None):
    df = df.sort_values(['Start', 'Bin'])
    bins = sorted(df['Bin'].unique())
    low_bins = [b for b in bins if int(b.replace('pct', '')) <= 50]
    high_bins = [b for b in bins if int(b.replace('pct', '')) > 50]

    all_starts = sorted(df['Start'].unique())
    midpoints = []
    for s in all_starts:
        midpoints.append(df[df['Start'] == s][['Start', 'End']].iloc[0].mean())
    x = np.array(midpoints)

    for bin_group, suffix in [(low_bins, 'low'), (high_bins, 'high')]:
        if not bin_group: continue

        pivot_df = df[df['Bin'].isin(bin_group)].pivot(index='Start', columns='Bin', values='Count').fillna(0)
        pivot_df = pivot_df.reindex(all_starts).fillna(0)

        y_stack = np.cumsum(pivot_df.values, axis=1)

        fig, ax = plt.subplots(figsize=(PLOT_CONFIG['width'], PLOT_CONFIG['height_mass']))
        plt.subplots_adjust(left=PLOT_CONFIG['left_margin'], right=PLOT_CONFIG['right_margin'],
                            bottom=PLOT_CONFIG['bottom_margin'], top=PLOT_CONFIG['top_margin'])

        y_prev = np.zeros(len(pivot_df))
        # Hot and cold color schemes added.
        cmap_hot = plt.get_cmap('YlOrRd', len(bin_group))
        cmap_cold = plt.get_cmap('GnBu', len(bin_group))
        for i in range(len(bin_group)):
            ax.fill_between(x, y_prev, y_stack[:, i], color=cmap_hot(i) if suffix=="high" else cmap_cold(i), step='mid', linewidth=0)
            y_prev = y_stack[:, i]

        setup_axes(ax, x_min, x_max, show_axes)
        if y_max:
            ax.set_ylim(0, y_max)
        else:
            current_max = y_stack.max()
            ax.set_ylim(0, max(current_max * 1.1, 1))

        ax.set_ylabel('Counts', fontsize=8)
        plt.savefig(f"{output_prefix}_{sequence}_{suffix}.png", dpi=PLOT_CONFIG['dpi'])
        plt.close()

def plot_gff(gff_path, output_path, sequence, x_min, x_max, tracks=None):
    if not os.path.exists(gff_path): return

    try:
        gff = pd.read_csv(gff_path, sep='\t', comment='#', header=None,
                         names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    except Exception as e:
        logger.error(f"Error reading GFF {gff_path}: {e}")
        return

    gff = gff[gff['seqid'] == sequence]
    if tracks:
        gff = gff[gff['type'].isin(tracks)]

    unique_types = sorted(gff['type'].unique())
    num_types = len(unique_types)
    if num_types == 0: return

    fig, axes = plt.subplots(nrows=num_types, figsize=(PLOT_CONFIG['width'], PLOT_CONFIG['height_gff'] * num_types), squeeze=False)
    plt.subplots_adjust(left=PLOT_CONFIG['left_margin'], right=PLOT_CONFIG['right_margin'],
                        bottom=0.1, top=0.9, hspace=0.1)

    for i, t in enumerate(unique_types):
        ax = axes[i, 0]
        type_gff = gff[gff['type'] == t]
        for _, row in type_gff.iterrows():
            color = 'blue' if row['strand'] == '+' else 'red' if row['strand'] == '-' else 'black'
            rect = Rectangle((row['start'], 0.1), row['end'] - row['start'], 0.8, color=color, linewidth=0)
            ax.add_patch(rect)

        ax.set_ylabel(t, rotation=0, ha='right', va='center', fontsize=7)
        setup_axes(ax, x_min, x_max, show_axes=(i == num_types - 1))
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)

    plt.savefig(output_path, dpi=PLOT_CONFIG['dpi'], transparent=True)
    plt.close()

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)
    if not os.path.exists(args.input_tsv):
        logger.error(f"Could not find input file '{args.input_tsv}'")
        return

    try:
        df = pd.read_csv(args.input_tsv, sep='\t')
    except Exception as e:
        logger.error(f"Error reading {args.input_tsv}: {e}")
        return

    if hasattr(args, 'html') and args.html:
        spectral.validate(df)
        if args.sequence:
            sequences = args.sequence.split(',')
            df = df.loc[df['Sequence'].isin(sequences)]

        if hasattr(args, 'zoom_width') and args.zoom_width:
            zoom = args.zoom_width.split(',')
            df = df.loc[df['Start'] >= int(zoom[0])]
            df = df.loc[df['End'] <= int(zoom[1])]

        spectra_groups = df.groupby(['Sequence'])
        for name, group in spectra_groups:
            libraries = group.Library.unique()
            for lib in libraries:
                lib_group = group[group['Library'] == lib]

                bases = ["A", "C", "G", "T"]
                queries = [f"{a}{b}{c}" for a in bases for b in bases for c in bases]
                present_queries = [q for q in queries if q in lib_group.columns]

                df_melted = pd.melt(
                    lib_group,
                    id_vars=['Library', 'Sequence', 'Start', 'End'],
                    value_vars=present_queries,
                    var_name="Triplet",
                    value_name="Count"
                )
                query_colors = {a: palette_builder(a) for a in present_queries}

                fig = plotly_plot(df_melted, query_colors, present_queries)
                output_name = f"{lib}_{name}_plot.html"
                fig.write_html(output_name)
                logger.info(f"Saved interactive plot to {output_name}")
        return

    is_mass = 'Bin' in df.columns and 'Count' in df.columns

    if args.sequence:
        sequences = args.sequence.split(',')
        df = df[df['Sequence'].isin(sequences)]
    else:
        sequences = df['Sequence'].unique()

    y_max_mass = None
    y_max_spectra = None
    if is_mass:
        # Calculate global y_max for consistent scaling across sequences
        bins = df['Bin'].unique()
        low_bins = [b for b in bins if int(str(b).replace('pct', '')) <= 50]
        high_bins = [b for b in bins if int(str(b).replace('pct', '')) > 50]

        max_low = 0
        if low_bins:
            max_low = df[df['Bin'].isin(low_bins)].groupby(['Sequence', 'Start'])['Count'].sum().max()
        max_high = 0
        if high_bins:
            max_high = df[df['Bin'].isin(high_bins)].groupby(['Sequence', 'Start'])['Count'].sum().max()

        y_max_mass = max(max_low, max_high) * 1.1 if not np.isnan(max(max_low, max_high)) else None
    else:
        if args.frequencies:
            # For frequencies, we might want consistent Y scale if they are not all 0-1
            queries = [c for c in df.columns if len(c) == 3 and all(b in 'ACGT' for b in c)]
            if queries:
                y_max_spectra = df[queries].sum(axis=1).max() * 1.1

    for seq in sequences:
        seq_df = df[df['Sequence'] == seq]
        if seq_df.empty: continue
        x_min = seq_df['Start'].min()
        x_max = seq_df['End'].max()
        if is_mass:
            plot_mass(seq_df, args.output, seq, x_min, x_max, show_axes=args.axes, y_max=y_max_mass)
        else:
            out_path = f"{args.output}_{seq}.png"
            plot_spectra(seq_df, out_path, seq, x_min, x_max, show_axes=args.axes, frequencies=args.frequencies, ngaps_path=args.ngaps if hasattr(args, 'ngaps') else None, y_max=y_max_spectra)
            if hasattr(args, 'gff_file') and args.gff_file:
                gff_tracks = args.gff_tracks.split(',') if args.gff_tracks else None
                plot_gff(args.gff_file, f"{args.output}_gff_{seq}.png", seq, x_min, x_max, tracks=gff_tracks)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input_tsv', required=True)
    parser.add_argument('-o', '--output', dest='output', required=True)
    parser.add_argument('-s', '--sequence', dest='sequence')
    parser.add_argument('-a', '--axes', action='store_true', default=False)
    parser.add_argument('-f', '--freq', dest='frequencies', action='store_true', default=False)
    parser.add_argument('--gff-file', dest='gff_file')
    parser.add_argument('--gff-tracks', dest='gff_tracks')
    parser.add_argument('--ngaps', dest='ngaps')
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()

    execute(args)
