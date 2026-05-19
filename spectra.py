#!/usr/bin/env python3

import argparse
import importlib
import sys
import os

def main():
    parser = argparse.ArgumentParser(description='Spectra genetic profiling. Counting, processing, and visualization of 3-mers')
    subparsers = parser.add_subparsers(title='commands', dest='command', required=True)

    # Count
    parser_count = subparsers.add_parser("count", description="Generate tsv file of spectra counts")
    parser_count.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
    parser_count.add_argument('-f', '--format', dest='sequence_format', type=str, help='Input file type', default='fasta')
    parser_count.add_argument('-w', '--width', dest='width', type=int, help='Window width', default=10000)
    parser_count.add_argument('-s', '--spacing', dest='spacing', type=int, help='Window spacing', default=10000)
    parser_count.add_argument('-o', '--output', dest='output', type=str, help='Output tsv file', default='spectra_report.tsv')
    parser_count.add_argument('-c', '--complement', dest='complement', action='store_true', help='Complement sequence file name. If set, calculates spectra for sequence complement (not reversed-complemented)', default=False)
    parser_count.add_argument('-l', '--libraries', dest='libraries', action='store_true', help='Sequence names include multiple libraries, prefixed by LIBRARY_', default=False)
    parser_count.add_argument('-p', '--proportions', dest='proportions', action='store_true', help='Return Spectra 3-mer proportions instead of raw counts', default=False)
    parser_count.add_argument('-k', '--kmer-size', dest='mer_size', type=int, help='kmer size to tabulate.', default=3)
    parser_count.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
    parser_count.add_argument('-n', '--no-overlap', dest='overlap', action='store_false', help='Count base pairs in repetitive runs of nucleotides only once.', default=True)
    parser_count.add_argument('-z', '--chunk-size', dest='chunk_size', type=int, help='Max chunk size to work on', default=30000000)
    parser_count.add_argument('-t', '--threads', dest='threads', type=int, help='Number of threads for parallel processing', default=1)
    parser_count.add_argument('--minimum-size', dest='minimum_size', type=int, help='Minimum sequence size to include.', default=15000)

    # Query
    parser_query = subparsers.add_parser("query", description="Generate tsv file of spectra counts for specific motifs")
    parser_query.add_argument('-i', '--input', dest='input_sequence', type=str, help='Input sequence file', required=True)
    parser_query.add_argument('-f', '--format', dest='sequence_format', type=str, help='Input file type', default='fasta')
    parser_query.add_argument('-q', '--query', dest='query', type=str, help='Query sequences, separated by commas', required=True)
    parser_query.add_argument('-w', '--width', dest='width', type=int, help='Window width', default=3000)
    parser_query.add_argument('-s', '--spacing', dest='spacing', type=int, help='Window spacing', default=3000)
    parser_query.add_argument('-o', '--output', dest='output', type=str, help='Output tsv file', default='spectra_report.tsv')
    parser_query.add_argument('-l', '--libraries', dest='libraries', action='store_true', help='Sequence names include multiple libraries, prefixed by LIBRARY_', default=False)
    parser_query.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
    parser_query.add_argument('-m', '--memory', dest='memory', action='store_true', help='Use memory-conservation mode', default=False)
    parser_query.add_argument('-d', '--consolidate', dest='consolidate', action='store_true', help='Consolidate all records into single column', default=False)
    parser_query.add_argument('-c', '--complement', dest='complement', action='store_true', help='Complement sequence file name. If set, calculates spectra for sequence complement (not reversed-complemented)', default=False)
    parser_query.add_argument('-n', '--no-overlap', dest='overlap', action='store_false', help='Count base pairs in repetitive runs of nucleotides only once.', default=True)
    parser_query.add_argument('-z', '--chunk-size', dest='chunk_size', type=int, help='Max chunk size to work on', default=30000000)
    parser_query.add_argument('-t', '--threads', dest='threads', type=int, help='Number of threads for parallel processing', default=1)

    # Collate
    parser_collate = subparsers.add_parser('collate', description='Collate multiple spectra output tsv into a multi-library tsv')
    parser_collate.add_argument('-i', '--input', dest='input_tsvs', help='Input spectra tsvs, separated by spaces', nargs='*', required=True)
    parser_collate.add_argument('-o', '--output', dest='output', help='Output spectra tsv', default='collated_spectra.tsv')
    parser_collate.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)

    # Transform
    parser_transform = subparsers.add_parser('transform', description='Transform spectra data for additional insight')
    parser_transform.add_argument('-i', '--input', dest='input_tsv', type=str, help='Input spectra tsv', required=True)
    parser_transform.add_argument('-o', '--output', dest='output', type=str, help='Output spectra tsv', required=True)
    parser_transform.add_argument('-r', '--weighted-filter', dest='weighted_filter', action='store_true', help='Produce two additional outputs that have outlier windows and normal windows', default=False)
    parser_transform.add_argument('-n', '--weighted-norm', dest='weighted_normalization', action='store_true', help='Normalize spectra frequencies for each window by the frequencies for the whole sequence', default=False)
    parser_transform.add_argument('-f', '--freq', dest='frequencies', action='store_true', help='Mark this if Spectra data is already in frequencies', default=False)
    parser_transform.add_argument('-c', '--convert', dest='convert', action='store_true', help='Convert between counts and frequencies', default=False)
    parser_transform.add_argument('-s', '--window-resize', dest='resize_window', type=int, help='Resize windows to summarize N for every 1 window')
    parser_transform.add_argument('-p', '--print', dest='print', action='store_true', help='Print global frequencies', default=False)
    parser_transform.add_argument('-y', '--simplify', dest='simplify', action='store_true', help='Simplify forward and reverse-complement counts per window', default=False)
    parser_transform.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='verbose mode', default=False)

    # Plot
    parser_plot = subparsers.add_parser('plot', description='Plot spectra profiles')
    parser_plot.add_argument('-i', '--input', dest='input_tsv', type=str, help='Input spectra or mass-query tsv', required=True)
    parser_plot.add_argument('-o', '--output', dest='output', type=str, help='Output prefix or filename', required=True)
    parser_plot.add_argument('-z', '--zoom', dest='zoom_width', type=str, help='Plot only a portion of the windows from between X,Y')
    parser_plot.add_argument('-s', '--sequence', dest='sequence', type=str, help='Plot only sequences matching Name1,Name2,Name3')
    parser_plot.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
    parser_plot.add_argument('-a', '--axes', dest='axes', action='store_true', help='Display axes', default=False)
    parser_plot.add_argument('-f', '--freq', dest='frequencies', action='store_true', help='Data is frequencies', default=False)
    parser_plot.add_argument('--gff-file', dest='gff_file', type=str, help='GFF file for annotations')
    parser_plot.add_argument('--gff-tracks', dest='gff_tracks', type=str, help='GFF tracks to include')
    parser_plot.add_argument('--ngaps', dest='ngaps', type=str, help='GFF of N-gap coordinates')
    parser_plot.add_argument('--html', dest='html', action='store_true', help='Generate interactive HTML plot (Plotly)', default=False)

    # Analyze
    parser_analyze = subparsers.add_parser('analyze', description='Analyze spectra profiles to detect breakpoints')
    parser_analyze.add_argument('-i', '--input', dest='input_tsv', type=str, help='Input spectra tsv', required=True)
    parser_analyze.add_argument('-o', '--output', dest='output_prefix', type=str, help='Output prefix', default=False)
    parser_analyze.add_argument('-p', '--penalty', dest='penalty', type=float, help='Ruptures breakpoint penalty criterion.', default=1000000)
    parser_analyze.add_argument('-a', '--aligned', dest='is_aligned', action='store_true', help='Check for if input tsv comes from alignment or from sequence data', default=False)
    parser_analyze.add_argument('-s', '--size', dest='size', type=int, help='Minimum windows to be considered a novel segment', default=5)
    parser_analyze.add_argument('-f', '--frequencies', dest='frequency', action='store_true', help='Process breaks by frequencies instead of raw counts', default=False)
    parser_analyze.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)

    args = parser.parse_args()

    try:
        module_name = f"scripts.{args.command}"
        script = importlib.import_module(module_name)
        script.execute(args)
    except ImportError as e:
        print(f"Error: Could not find script for command '{args.command}': {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error executing command '{args.command}': {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
