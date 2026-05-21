#!/usr/bin/env python3
import argparse
import math
import tempfile
import heapq
import os
import logging

def stream_compute_reductions(raw_file, asm_file, out_file, chunk_size=5000000, extreme=None):
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
                if not line_r or not line_a: break
                if line_r[0] == line_a[0]:
                    kmer, r_count, a_count = line_r[0], int(line_r[1]), int(line_a[1])
                    # Methodology: normalize by total counts if needed?
                    # For now, stick to log-fold change but ensure matched kmers are handled correctly.
                    reduction = math.log10(a_count + 1) - math.log10(r_count + 1)
                    buffer.append((reduction, kmer, r_count, a_count))

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
        if extreme is not None:
            # We need to know the total count to apply extreme filtering
            # Since merged is an iterator, we might need a list or a two-pass if it's large.
            # However, external merge sort usually implies it's too large for memory.
            # Let's write everything to a temp file first to get the count if extreme is set.
            temp_ranked = tempfile.NamedTemporaryFile(delete=False, mode="w")
            count = 0
            for red, k, rc, ac in merged:
                temp_ranked.write(f"{red}\t{k}\t{rc}\t{ac}\n")
                count += 1
            temp_ranked.close()

            low_bound = count * (extreme / 100.0)
            high_bound = count * (1 - extreme / 100.0)

            with open(out_file, "w") as out:
                out.write("kmer\tRaw\tAsm\treduction\treductionRank\n")
                with open(temp_ranked.name) as tr:
                    rank = 1
                    for line in tr:
                        if rank <= low_bound or rank > high_bound:
                            red, k, rc, ac = line.strip().split("\t")
                            out.write(f"{k}\t{rc}\t{ac}\t{float(red):.6f}\t{rank}\n")
                        rank += 1
            os.remove(temp_ranked.name)
        else:
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
    parser.add_argument("-r", "--raw_dump", required=True, help="Raw k-mer dump file")
    parser.add_argument("-a", "--asm_dump", required=True, help="Assembly k-mer dump file")
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

    stream_compute_reductions(args.raw_dump, args.asm_dump, args.output, args.chunk_size, extreme=args.extreme)
    logger.info(f"Ranked table written to {args.output}")

if __name__ == "__main__":
    main()
