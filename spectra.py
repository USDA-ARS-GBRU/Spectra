#!/usr/bin/env python3

# this is the script to house all components
# Core imports: csv, argparse, and os will be needed for all three arms
import csv
import argparse
import importlib
import os

# TODO: modify argparse to be a universal set of inputs
parser = argparse.ArgumentParser(description='Spectra genetic profiling. Counting, processing, and visualization of 3-mers')
subparser = parser.add_subparsers(title='spectra', dest='subparser', required=True)

parserCount = subparser.add_parser("count", description="Generate tsv file of spectra counts")
parserCount.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
parserCount.add_argument('-f', '--format', dest='sequence_format', type=str, help='Input file type', default='fasta')
parserCount.add_argument('-w', '--width', dest='width', type=int, help='Window width', default='3000')
parserCount.add_argument('-s', '--spacing', dest='spacing', type=int, help='Window spacing', default='3000')
parserCount.add_argument('-o', '--output', dest='output', type=str, help='Output tsv file', default='spectra_report.tsv')
parserCount.add_argument('-l', '--libraries', dest='libraries', action='store_true', help='Sequence names include multiple libraries, prefixed by LIBRARY_', default=False)
parserCount.add_argument('-p', '--proportions', dest='proportions', action='store_true', help='Return Spectra 3-mer proportions instead of raw counts', default=False)

parserCollate = subparser.add_parser('collate', description='Collate multiple spectra output tsv into a multi-library tsv')
parserCollate.add_argument('-i', '--input', dest='input_tsvs', help='Input spectra tsvs, separated by spaces', nargs='*', required=True)
parserCollate.add_argument('-o', '--output', dest='output', help='Output spectra tsv', default='collated_spectra.tsv')

parserTransform = subparser.add_parser('transform', description='Transform spectra data for additional insight')
parserTransform.add_argument('-i', '--input', dest='input_tsv', type=str, help='Input spectra tsv', required=True)
parserTransform.add_argument('-o', '--output', '--output', dest='output', type=str, help='Output spectra tsv', default='transformed_output.tsv')
parserTransform.add_argument('-r', '--weighted-rm', dest='weighted_removal', action='store_true', help='Remove windows with spectra frequencies similar to the whole sequence', default=False)
parserTransform.add_argument('-n', '--weighted-norm', dest='weighted_normalization', action='store_true', help='Normalize spectra frequencies for each window by the frequencies for the whole sequence', default=False)
parserTransform.add_argument('-f', '--freq', dest='frequencies', action='store_true', help='Mark this is Spectra data is already in frequencies', default=False)
parserTransform.add_argument('-s', '--window-resize', dest='resize_window', type=int, help='Resize windows to summarize N for every 1 window')

parserPlot = subparser.add_parser('plot', description='Plot spectra profiles')
parserPlot.add_argument('-i', '--input', dest='input_tsv', type=str, help='Input spectra tsv', required=True)
parserPlot.add_argument('-o', '--output', dest='output', type=str, help='Output spectra plot', default='spectra_plot.png')
parserPlot.add_argument('-z', '--zoom', dest='zoom_width', type=str, help='Plot only a portion of the windows from between X,Y')
parserPlot.add_argument('-s', '--sequence', dest='sequence', type=str, help='Plot only sequences matching Name1,Name2,Name3')
parserPlot.add_argument('-g', '--gff-file', dest='gff_file', type=str, help='Plot annotations from gff')
parserPlot.add_argument('-t', '--gff-tracks', dest='gff_tracks', type=str, help='Plot only annotations matching category Name1,Name2,Name3')
parserPlot.add_argument('-l', '--legend', dest='show_legend', action='store_true', help='Draw legend', default=False)
parserPlot.add_argument('-r', '--dpi', dest='image_resolution', type=int, help='Image resolution in DPI', default=300)

args = parser.parse_args()

scriptDirectory = os.path.dirname(__file__)
subparserScript = scriptDirectory + "/" + args.subparser + ".py"

if os.path.exists(subparserScript):
    script = importlib.import_module(args.subparser)
    script.execute(args)
# Set a switch here that determines which subprogram to run

#TODO: create Spectra-count functions

#TODO: create Spectra-collate functions

#TODO: create Spectra-plot functions

#TODO: create Spectra-mutate functions