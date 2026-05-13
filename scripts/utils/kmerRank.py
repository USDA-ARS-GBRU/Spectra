#!/usr/bin/env python3
import argparse
import math
import tempfile
import heapq
import os
import logging

def stream_compute_reductions(raw_file, asm_file, out_file, chunk_size=5000000):
    """
    Stream merge raw and asm files, compute log-fold reductions,
    write them to disk in sorted chunks, then external-merge them into ranked output.
    """
    tmp_files = []

    # Pass 1: stream merge + chunked sort
    try:
        with open(raw_file) as fr, open(asm_file) as fa:
            line_r = fr.readline().split()
            line_a = fa.readline().split()
            buffer = []

            while line_r and line_a:
                if line_r[0] == line_a[0]:
                    kmer, rc, ac = line_r[0], int(line_r[1]), int(line_a[1])
                    reduction = math.log10(ac + 1) - math.log10(rc + 1)
                    buffer.append((reduction, kmer, rc, ac))

                    if len(buffer) >= chunk_size:
                        buffer.sort(key=lambda x: x[0])
                        tf = tempfile.NamedTemporaryFile(delete=False, mode="w")
                        for red, k, rcount, acount in buffer:
                            tf.write(f"{red}\t{k}\t{rcount}\t{acount}\n")
                        tf.close()
                        tmp_files.append(tf.name)
                        buffer.clear()

                    line_r = fr.readline().split()
                    line_a = fa.readline().split()

                elif line_r[0] < line_a[0]:
                    line_r = fr.readline().split()
                else:
                    line_a = fa.readline().split()

            # flush last buffer
            if buffer:
                buffer.sort(key=lambda x: x[0])
                tf = tempfile.NamedTemporaryFile(delete=False, mode="w")
                for red, k, rcount, acount in buffer:
                    tf.write(f"{red}\t{k}\t{rcount}\t{acount}\n")
                tf.close()
                tmp_files.append(tf.name)
                buffer.clear()
    except Exception as e:
        logging.error(f"Error during reduction computation: {e}")
        # Cleanup any temporary files created before returning
        for f in tmp_files:
            if os.path.exists(f): os.remove(f)
        return

    # Pass 2: external merge of sorted chunks
    def file_iter(fname):
        with open(fname) as f:
            for line in f:
                red, k, rc, ac = line.strip().split("\t")
                yield float(red), k, int(rc), int(ac)

    iterators = [file_iter(f) for f in tmp_files]
    merged = heapq.merge(*iterators, key=lambda x: x[0])

    # Pass 3: write ranked table
    try:
        with open(out_file, "w") as out:
            out.write("kmer\tRaw\tAsm\treduction\treductionRank\n")
            rank = 1
            for red, k, rc, ac in merged:
                out.write(f"{k}\t{rc}\t{ac}\t{red:.6f}\t{rank}\n")
                rank += 1
    except Exception as e:
        logging.error(f"Error writing ranked table to {out_file}: {e}")

    # cleanup
    for f in tmp_files:
        if os.path.exists(f):
            os.remove(f)


def main():
    # CLI arguments
    parser = argparse.ArgumentParser(description="Rank kmers by log-fold change between raw and assembly")
    parser.add_argument("-r", "--raw_dump", required=True, help="Raw jellyfish dump file")
    parser.add_argument("-a", "--asm_dump", required=True, help="Assembly jellyfish dump file")
    parser.add_argument("-o", "--output", required=True, help="Output ranked table (TSV)")
    parser.add_argument("-c", "--chunk_size", type=int, default=5000000,
                        help="Number of kmers to hold in memory before spilling to disk [default 5,000,000]")
    parser.add_argument("-e", "--extreme", type=float, default=None,
                        help="Keep only the top/bottom PERCENT of reductions (e.g. 5 = keep 5%% lowest and 5%% highest)")
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)
    args = parser.parse_args()

    # Logging
    logging.basicConfig(level=logging.ERROR, format='%(levelname)s: %(message)s')
    logger = logging.getLogger()
    if args.verbose:
        logger.setLevel(logging.INFO)

    stream_compute_reductions(args.raw_dump, args.asm_dump, args.output, args.chunk_size)
    logger.info(f"Ranked table written to {args.output}")

if __name__ == "__main__":
    main()
