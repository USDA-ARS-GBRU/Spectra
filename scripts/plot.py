#!/usr/bin/env python3

import os
import pandas as pd
import logging
import spectral
import plotly.express as px

logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')
logger = logging.getLogger()

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
        return f"#{''.join([colors[bases.index(a)] for a in triplet])}"
    if palette == "gc":
        gc_scale = [("D9", "3C"), ("C6", "5F"), ("B3", "76"), ("A0", "90")]
        gc_sum = 3-sum([1 for a in triplet if a in ["A", "T"]])
        ac_scale = ["33", "46", "59", "73"]
        ac_sum = 3-sum([1 for a in triplet if a in ["A", "C"]])
        return f"#{gc_scale[gc_sum][0]}{ac_scale[ac_sum]}{gc_scale[gc_sum][1]}"
    return "#000000"

def colorize_mer(mer):
    colors = {"A": "C6", "C": "6C", "G": "3C", "T": "10"}
    return "#" + colors[mer[2]] + colors[mer[1]] + colors[mer[0]]

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)
    if not os.path.exists(args.input_tsv):
        logger.error(f"Could not find input file '{args.input_tsv}'")
        return

    try:
        spectra_df = pd.read_csv(args.input_tsv, sep='\t')
    except Exception as e:
        logger.error(f"Error reading {args.input_tsv}: {e}")
        return

    # plot_colors = [colorize_mer(a) for a in spectra_df.columns if a[0] in ['A', 'C', 'G', 'T']]

    spectral.validate(spectra_df)
    if args.sequence:
        sequences = args.sequence.split(',')
        spectra_df = spectra_df.loc[spectra_df['Sequence'].isin(sequences)]

    if args.zoom_width:
        zoom = args.zoom_width.split(',')
        spectra_df = spectra_df.loc[spectra_df['Start'] >= int(zoom[0])]
        spectra_df = spectra_df.loc[spectra_df['End'] <= int(zoom[1])]

    spectra_groups = spectra_df.groupby(['Sequence'])
    for name, group in spectra_groups:
        libraries = group.Library.unique()
        for lib in libraries:
            lib_group = group[group['Library'] == lib]

            ### Plotly
            bases = ["A", "C", "G", "T"]
            # Hardcoded to 3-mer for now in original code
            queries = [f"{a}{b}{c}" for a in bases for b in bases for c in bases]

            # Filter queries to those present in the dataframe
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
