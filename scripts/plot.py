#!/usr/bin/env python3

import os
import pandas as pd
import logging
import spectral
import plotly.express as px
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()
def plotlyPlot(df_melted, query_colors, queries):
    fig = px.area(df_melted, x="Start", y="Count", color="Triplet", color_discrete_map=query_colors,
                  line_group="Triplet", category_orders={"Triplet": list(reversed(queries))})
    # fig = px.area(df_melted, x="Start", y="Count", color=query_colors, line_group="Triplet")
    fig.for_each_trace(lambda trace: trace.update(fillcolor=trace.line.color))
    # fig.update_layout(margin_b=0, margin_t=0, margin_l=0, margin_r=0, xaxis=dict(rangeslider=dict(visible=True)),yaxis=dict(autorange="reversed"))
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


def colorizeMer(mer):
    colors = {"A": "C6", "C": "6C", "G": "3C", "T": "10"}
    return "#" + colors[mer[2]] + colors[mer[1]] + colors[mer[0]]

def execute(args):
    if args.verbose:
        logger.setLevel(logging.INFO)
    if not os.path.exists(args.input_tsv):
        logging.error(f"Could not find input file '{args.input_tsv}'")
        exit()

    spectraPanda = pd.read_csv(args.input_tsv, sep='\t')
    plotColors = [colorizeMer(a) for a in spectraPanda.columns if a[0] in ['A', 'C', 'G', 'T']]

    validation = spectral.validate(spectraPanda)
    if args.sequence:
        sequences = args.sequence.split(',')
        spectraPanda = spectraPanda.loc[spectraPanda['Sequence'].isin(sequences)]

    if args.zoom_width:
        zoom = args.zoom_width.split(',')
        spectraPanda = spectraPanda.loc[spectraPanda['Start'] >= int(zoom[0])]
        spectraPanda = spectraPanda.loc[spectraPanda['End'] <= int(zoom[1])]

    spectraGroups = spectraPanda.groupby(['Sequence'])
    for group in spectraGroups:
        library = group[1].Library.unique()
        if len(library) > 1:
            pass
        else:
            library = library[0]
            sequence = group[1].Sequence.unique()[0]


            ### Plotly test
            bases = ["A", "C", "G", "T"]
            queries = [f"{a}{b}{c}" for a in bases for b in bases for c in bases]
            df_melted = pd.melt(
                group[1],
                id_vars=['Library', 'Sequence', 'Start', 'End'],
                value_vars=group[1].columns[4:],
                var_name="Triplet",
                value_name="Count"
            )
            query_colors = {a: palette_builder(a) for a in queries}

            fig = plotlyPlot(df_melted, query_colors, queries)
            fig.write_html(f"{library}_{sequence}_plot.html")

