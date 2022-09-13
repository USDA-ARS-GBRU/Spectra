#!/usr/bin/env python3

import os
import pandas as pd
from plotnine import ggplot, aes, geom_area, ggsave, theme, geom_segment, scale_y_continuous, scale_x_continuous, xlab, ylab, theme_bw, scale_fill_manual, themes
import logging

logging.basicConfig(level=logging.INFO)

def colorizeMer(mer):
    colors = {"A": "DD", "G": "9F", "C": "60", "T": "11"}
    return "#" + colors[mer[0]] + colors[mer[1]] + colors[mer[2]]

#('-i', '--input', dest='input_tsv', type=str, help='Input spectra tsv', required=True)
#('-o', '--output', dest='output', type=str, help='Output spectra plot', default='spectra_plot.png')
#('-z', '--zoom', dest='zoom_width', type=str, help='Plot only a portion of the windows from between X,Y')
#('-s', '--sequence', dest='sequence', type=str, help='Plot only sequences matching Name1,Name2,Name3')
#('-g', '--gff-file', dest='gff_file', type=str, help='Plot annotations from gff')
#('-t', '--gff-tracks', dest='gff_tracks', type=str, help='Plot only annotations matching category Name1,Name2,Name3')
#('-l', '--legend', dest='show_legend', action='store_true', help='Draw legend', default=False)
#('-r', '--dpi', dest='image_resolution', type=int, help='Image resolution in DPI', default=300)

def execute(args):
    if not os.path.exists(args.input_tsv):
        logging.error(f"Could not find input file '{args.input_tsv}'")
        exit()
    spectraPanda = pd.read_csv(args.input_tsv, sep='\t')
    plotColors = [colorizeMer(a) for a in spectraPanda.columns[-64:]]
    library = spectraPanda.Library.unique()[0]
    libraryPanda = spectraPanda[spectraPanda['Library'] == library]
    sequence = libraryPanda.Sequence.unique()[1]
    sequencePanda = libraryPanda[libraryPanda['Sequence'] == sequence]

    sequencePanda = sequencePanda.melt(id_vars=['Library', 'Sequence', 'Start', 'End'], value_vars=sequencePanda.columns[-64:])

    # TODO: add boolean check if these are counts or frequencies
    windowSize = sequencePanda['End'][0] - sequencePanda['Start'][0] - 1
    subplotSize = (min(sequencePanda['Start']), max(sequencePanda['End']))

    plot = ggplot(data=sequencePanda)
    plot += aes(fill='variable', x=(sequencePanda['End'] + sequencePanda['Start']) / 2, y=sequencePanda['value'] / windowSize)
    plot += geom_area(stat="identity", position="stack")
    plot += scale_fill_manual(values=plotColors)
    plot += scale_y_continuous(limits=(0, 1), expand=(0, 0))
    plot += scale_x_continuous(limits=subplotSize, expand=(0, 0))
    plot += xlab("Window Position (nucleotide)")
    plot += ylab("Proportion")
    plot += theme_bw()
    plot += theme(axis_title_x=themes.element_text(size=6), axis_title_y=themes.element_text(size=6), legend_text=themes.element_text(size=5), legend_key_size=5, legend_spacing=-3, axis_text_x=themes.element_text(angle=0, vjust=0.5), line=themes.element_blank(), axis_ticks=themes.element_line(), legend_title=themes.element_blank())
    ggsave(plot, filename=library+"_"+sequence+".png", width=12, height=4, units="in", dpi=args.image_resolution, limitsize=False, verbose=False)

    # TODO: Plot window size needs to be wider
    # TODO: Fix DPI
    # TODO: Color values by colorizeMer()
    # TODO: Fix limits
    # TODO: Change values to frequencies instead of counts
    # TODO: Legend needs to be optional and moved
    # TODO: Fix graph padding/margin



    #seqPlot = sequencePanda.plot.area(x='Start', y=list(sequencePanda.columns)[-64:], stacked=True, figsize=(args.image_resolution / 15, args.image_resolution / 40))
    #outPlot = seqPlot.get_figure()
    #outPlot.savefig(library+"_"+sequence+".png")
    #plt.close('all')
    # for library in spectraPanda.Library.unique():
    #     libraryPanda = spectraPanda[spectraPanda['Library'] == library]
    #     for sequence in libraryPanda.Sequence.unique():
    #         sequencePanda = libraryPanda[libraryPanda['Sequence'] == sequence]
    #         sequencePanda = sequencePanda.melt(id_vars=['Library', 'Sequence', 'Start', 'End'], value_vars=sequencePanda.head()[-64:])
    #         seqPlot = sequencePanda.plot.area(x='Start', y='value', stacked=True)
    #         outPlot = seqPlot.get_figure()
    #         outPlot.savefig(library+"_"+sequence+".png")
    #         outPlot.close()

