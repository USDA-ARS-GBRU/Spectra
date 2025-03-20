#!/usr/bin/env python3

import os
import pandas as pd
#import plotnine
#from plotnine import ggplot, aes, geom_area, ggsave, theme, geom_segment, scale_y_continuous, scale_x_continuous, xlab, ylab, theme_bw, scale_fill_manual, themes
import logging
import spectral
import plotly.express as px

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

            ### Plotly test

### Plotnine ggplot port - potentially deprecating
#            sequencePanda = group[1].melt(id_vars=validation[0], value_vars=validation[1])
#            windowSize = sequencePanda['End'][0] - sequencePanda['Start'][0] - 1
#            subplotSize = (min(sequencePanda['Start']), max(sequencePanda['End']))

#            plot = ggplot(data=sequencePanda)
#            plot += aes(fill='variable', x=(sequencePanda['End'] + sequencePanda['Start']) / 2, y=sequencePanda['value'] / (sequencePanda['End']-sequencePanda['Start']+1))
#            plot += geom_area(mapping=aes(), stat="identity", position=plotnine.position_stack, outline_type='lower')
#            plot += scale_fill_manual(values=plotColors)
#            plot += scale_y_continuous(limits=(0, 1), expand=(0, 0))
#            plot += scale_x_continuous(limits=subplotSize, expand=(0, 0))
#            plot += xlab("Window Position (nucleotide)")
#            plot += ylab("Proportion")
#            plot += theme_bw()
#            plot += theme(axis_title_x=themes.element_text(size=6), axis_title_y=themes.element_text(size=6), legend_text=themes.element_text(size=5), legend_key_size=5, legend_spacing=-3, axis_text_x=themes.element_text(angle=0, vjust=0.5), line=themes.element_blank(), axis_ticks=themes.element_line(), legend_title=themes.element_blank())
#            ggsave(plot, filename=f"{library}_{sequence}.png", width=12, height=4, units="in", dpi=args.image_resolution, limitsize=False, verbose=False)
### Plotnine ggplot port