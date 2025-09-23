#!/usr/bin/env python3
import argparse
import os
import time
from Bio import SeqIO, Seq
import csv
import sys
import logging

def find_N_regions(fasta_file, output_tsv):
    with open(output_tsv, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_name = record.id
            seq_str = str(record.seq).upper()

            i = 0
            while i < len(seq_str):
                if seq_str[i] == "N":
                    start = i + 1  # 1-based index
                    while i < len(seq_str) and seq_str[i] == "N":
                        i += 1
                    end = i  # still 1-based
                    out.write(f"{seq_name}\tn-counter\tngap\t{start}\t{end}\t.\t+\t.\tNGAP:{end-start+1}bp\n")
                else:
                    i += 1


parser = argparse.ArgumentParser(description="Kmer Mass Query: localize percentile kmers in genomic sequences (multi-pass)")
parser.add_argument('-i', '--input', dest='input', required=True, help='Input sequence file (FASTA/FASTQ)')
parser.add_argument('-f', '--format', default='fasta', help='Input file type [default fasta]')
parser.add_argument('-o', '--output', default='n-positions.gff', help='Output gff file')
parser.add_argument('-v', '--verbose', action='store_true', help='Verbose mode')
args = parser.parse_args()

# Logging
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()
if args.verbose:
    logger.setLevel(logging.INFO)
startTime = time.time()

if not os.path.exists(args.input):
    logger.error(f"Couldn't find input sequence file '{args.input}'")
    exit()


find_N_regions(args.input, args.output)
logger.info(f"Execution time in seconds: {time.time() - startTime:.2f}")