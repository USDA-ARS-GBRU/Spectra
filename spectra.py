#!/usr/bin/env python3

# this is the script to house all components
# Core imports: argparse, importlib,  and os will be needed for all three arms
import argparse
import importlib
# import os
# import logging
# import configparser


def moduleFromPath(path):
    return path[1:].replace(".py", "").replace("/", ".")


# TODO: add configure.ini for some hard-coded values like chi-square p-value
# TODO: modify argparse to be a universal set of inputs
parser = argparse.ArgumentParser(description='Spectra genetic profiling. Counting, processing, and visualization of 3-mers')
subparsers = parser.add_subparsers(title='spectra', dest='subparsers', required=True)

parserCount = subparsers.add_parser("count", description="Generate tsv file of spectra counts")
parserCount.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
parserCount.add_argument('-f', '--format', dest='sequence_format', type=str, help='Input file type', default='fasta')
parserCount.add_argument('-w', '--width', dest='width', type=int, help='Window width', default='3000')
parserCount.add_argument('-s', '--spacing', dest='spacing', type=int, help='Window spacing', default='3000')
parserCount.add_argument('-o', '--output', dest='output', type=str, help='Output tsv file', default='spectra_report.tsv')
parserCount.add_argument('-c', '--complement', dest='complement', action='store_true', help='Complement sequence file name. If set, calculates spectra for sequence complement (not reversed-complemented)', default=False)
parserCount.add_argument('-l', '--libraries', dest='libraries', action='store_true', help='Sequence names include multiple libraries, prefixed by LIBRARY_', default=False)
parserCount.add_argument('-p', '--proportions', dest='proportions', action='store_true', help='Return Spectra 3-mer proportions instead of raw counts', default=False)
parserCount.add_argument('-m', '--memory', dest='memory', action='store_true', help='Use memory-conservation mode', default=False)
parserCount.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
parserCount.add_argument('-n', '--no-overlap', dest='overlap', action='store_false', help='Verbose mode', default=True)

parserQuery = subparsers.add_parser("query", description="Generate tsv file of spectra counts")
parserQuery.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
parserQuery.add_argument('-f', '--format', dest='sequence_format', type=str, help='Input file type', default='fasta')
parserQuery.add_argument('-q', '--query', dest='query', type=str, help='Query sequences, separated by commas', required=True)
parserQuery.add_argument('-w', '--width', dest='width', type=int, help='Window width', default='3000')
parserQuery.add_argument('-s', '--spacing', dest='spacing', type=int, help='Window spacing', default='3000')
parserQuery.add_argument('-o', '--output', dest='output', type=str, help='Output tsv file', default='spectra_report.tsv')
parserQuery.add_argument('-l', '--libraries', dest='libraries', action='store_true', help='Sequence names include multiple libraries, prefixed by LIBRARY_', default=False)
parserQuery.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
parserQuery.add_argument('-m', '--memory', dest='memory', action='store_true', help='Use memory-conservation mode', default=False)
parserQuery.add_argument('-c', '--complement', dest='complement', action='store_true', help='Complement sequence file name. If set, calculates spectra for sequence complement (not reversed-complemented)', default=False)
parserQuery.add_argument('-n', '--no-overlap', dest='overlap', action='store_false', help='Verbose mode', default=True)


parserCollate = subparsers.add_parser('collate', description='Collate multiple spectra output tsv into a multi-library tsv')
parserCollate.add_argument('-i', '--input', dest='input_tsvs', help='Input spectra tsvs, separated by spaces', nargs='*', required=True)
parserCollate.add_argument('-o', '--output', dest='output', help='Output spectra tsv', default='collated_spectra.tsv')
parserCollate.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)


parserTransform = subparsers.add_parser('transform', description='Transform spectra data for additional insight')
parserTransform.add_argument('-i', '--input', dest='input_tsv', type=str, help='Input spectra tsv', required=True)
parserTransform.add_argument('-o', '--output', '--output', dest='output', type=str, help='Output spectra tsv', default=False)
parserTransform.add_argument('-r', '--weighted-filter', dest='weighted_filter', action='store_true', help='Produce two additional outputs that have outlier windows and normal windows', default=False)
parserTransform.add_argument('-n', '--weighted-norm', dest='weighted_normalization', action='store_true', help='Normalize spectra frequencies for each window by the frequencies for the whole sequence', default=False)
parserTransform.add_argument('-f', '--freq', dest='frequencies', action='store_true', help='Mark this is Spectra data is already in frequencies', default=False)
parserTransform.add_argument('-c', '--convert', dest='convert', action='store_true', help='Convert between counts and frequencies', default=False)
parserTransform.add_argument('-s', '--window-resize', dest='resize_window', type=int, help='Resize windows to summarize N for every 1 window')
parserTransform.add_argument('-p', '--print', dest='print', action='store_true', help='Print global frequencies', default=False)
parserTransform.add_argument('-y', '--simplify', dest='simplify', action='store_true', help='Simplify forward and reverse-complement counts per window', default=False)
parserTransform.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='verbose mode', default=False)


parserPlot = subparsers.add_parser('plot', description='Plot spectra profiles')
parserPlot.add_argument('-i', '--input', dest='input_tsv', type=str, help='Input spectra tsv', required=True)
parserPlot.add_argument('-o', '--output', dest='output', type=str, help='Output spectra plot', default='spectra_plot.png')
parserPlot.add_argument('-z', '--zoom', dest='zoom_width', type=str, help='Plot only a portion of the windows from between X,Y')
parserPlot.add_argument('-s', '--sequence', dest='sequence', type=str, help='Plot only sequences matching Name1,Name2,Name3')
parserPlot.add_argument('-g', '--gff-file', dest='gff_file', type=str, help='Plot annotations from gff')
parserPlot.add_argument('-t', '--gff-tracks', dest='gff_tracks', type=str, help='Plot only annotations matching category Name1,Name2,Name3')
parserPlot.add_argument('-l', '--legend', dest='show_legend', action='store_true', help='Draw legend', default=False)
parserPlot.add_argument('-r', '--dpi', dest='image_resolution', type=int, help='Image resolution in DPI', default=300)
parserPlot.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)


parserAnalyze = subparsers.add_parser('analyze', description='Analyze spectra profiles')
parserAnalyze.add_argument('-i', '--input', dest='input_tsv', type=str, help='Input spectra tsv', required=True)
parserAnalyze.add_argument('-o', '--output', dest='output_tsv', type=str, help='Output spectra tsv', default=False)
parserAnalyze.add_argument('-p', '--penalty', dest='penalty', type=float, help='Ruptures breakpoint penalty criterion', default=1000000)
parserAnalyze.add_argument('-a', '--aligned', dest='is_aligned', action='store_true', help='Check for if input tsv comes from alignment or from sequence data', default=False)
parserAnalyze.add_argument('-s', '--size', dest='size', type=int, help='Minimum windows to be considered a novel segment', default=5)
parserAnalyze.add_argument('-f', '--frequencies', dest='frequency', action='store_true', help='Process breaks by frequencies instead of raw counts', default=False)
parserAnalyze.add_argument('-b', '--blocked', dest='is_blocked', action='store_true', help='If data already has breakpoints, ', default=False)
parserAnalyze.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)

parserContigs = subparsers.add_parser('contigs', description='Find congruent contig joins')
parserContigs.add_argument('-i', '--input', dest='input_tsv', type=str, help='Input spectra tsv', required=True)
parserContigs.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)

parserCluster = subparsers.add_parser('cluster', description='Find congruent contig joins through k-means pairing')
parserCluster.add_argument('-i', '--input', dest='input_tsv', type=str, help='Input spectra tsv', required=True)
parserCluster.add_argument('-o', '--output', dest='output_tsv', type=str, help='Output spectra tsv')
parserCluster.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)

args = parser.parse_args()

subparsersScript = f"/scripts/{args.subparsers}.py"
script = importlib.import_module(moduleFromPath(subparsersScript))
script.execute(args)
