#!/usr/bin/env python3
import argparse
import os
import time
from Bio import SeqIO
import logging

def find_N_regions(fasta_file, output_gff):
    try:
        with open(output_gff, "w") as out:
            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_id = record.id
                seq_str = str(record.seq).upper()

                i = 0
                seq_len = len(seq_str)
                while i < seq_len:
                    if seq_str[i] == "N":
                        start = i + 1  # 1-based index
                        while i < seq_len and seq_str[i] == "N":
                            i += 1
                        end = i  # still 1-based
                        out.write(f"{seq_id}\tn-counter\tngap\t{start}\t{end}\t.\t+\t.\tNGAP:{end-start+1}bp\n")
                    else:
                        i += 1
    except Exception as e:
        logging.error(f"Error finding N regions: {e}")

def main():
    parser = argparse.ArgumentParser(description="n-counter: localize N-gaps in genomic sequences")
    parser.add_argument('-i', '--input', dest='input', required=True, help='Input sequence file (FASTA/FASTQ)')
    parser.add_argument('-f', '--format', default='fasta', help='Input file type [default fasta]')
    parser.add_argument('-o', '--output', default='n-positions.gff', help='Output gff file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose mode')
    args = parser.parse_args()

    # Logging
    logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')
    logger = logging.getLogger()
    if args.verbose:
        logger.setLevel(logging.INFO)

    start_time = time.time()

    if not os.path.exists(args.input):
        logger.error(f"Couldn't find input sequence file '{args.input}'")
        return

    find_N_regions(args.input, args.output)
    logger.info(f"Execution time in seconds: {time.time() - start_time:.2f}")

if __name__ == "__main__":
    main()
