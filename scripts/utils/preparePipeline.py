#!/usr/bin/env python3
import argparse
import os

# CLI arguments
parser = argparse.ArgumentParser(description="Prepare inputs for specta pipeline, and generate a bash script for running.")
parser.add_argument('-r', '--raw', dest='raw', required=True, nargs='+',help='Input raw fasta/fastq read file(s). These can be gzipped, but must end in ".gz". If multiple, separate with spaces')
parser.add_argument('-a', '--assembled', dest='assembled', required=True, help='Input fasta/gzipped fasta sequence assembly file')
parser.add_argument('-o', '--output-script', dest='output', default='spectra-pipeline.sh', help='Output bash file')
parser.add_argument('-p', '--output-prefix', dest='prefix', default='spectra_pipeline', help='Output files prefix. A directory will be created with this name for storing images')
parser.add_argument('-t', '--threads', dest='threads', type=int, help='processing threads [default 1]', default=1)
parser.add_argument('-k', '--kmer-size', dest='mer_size', type=int, help='kmer size in query [default 20]', default=20)
parser.add_argument('--jellyfish-bloom', dest='jf_bloom', type=str, default='100M', help='Jellyfish2 count bloomfilter initial size [default 100M]')
parser.add_argument('--jellyfish-path', dest='jf_path', type=str, default='jellyfish', help='Jellyfish2 path. Default assumes it is in your env [default jellyfish]')
parser.add_argument('--python-callable', dest='python', type=str, default='python', help='python3 path. Default assumes it is in your env [default python]')
parser.add_argument('--spectra-callable', dest='spectra', type=str, default=None, help='Spectra path. If not set, automatically detected from this script')
parser.add_argument('--rscript-callable', dest='rscript', type=str, default='Rscript', help='Rscript path. Default assumes it is in your env [default Rscript]')
parser.add_argument('-R', '--raw-min', dest='raw_min', type=int, default=100, help='Jellyfish2 raw kmer minimum count to retain [default 100]')
parser.add_argument('-A', '--asm-min', dest='asm_min', type=int, default=2, help='Jellyfish2 assembly kmer minimum count to retain [default 2]')
parser.add_argument('-m', '--mq-window', dest='mq_window', type=int, default=200000, help='Window and spacing width for kmer mass-query.py localization [default 200000]')
parser.add_argument('-w', '--spectra-window', dest='spectra_window', type=int, default=10000, help='Window and spacing width for spectra.py K=3 localization [default 10000]')
parser.add_argument('-c', '--clean-workspace', dest='clean', action='store_false', help='Clean workspace as files are processed. Jellyfish kmer counts are very large. By default, these files are removed after processing.', default=True)
args = parser.parse_args()

### Hard values. These will be modified to arguments when they have been tested further
kmer_sample_size = 5000000
chunk_size = 5000000
percentile_keep = 5

spectra_path = args.spectra if args.spectra else os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
if spectra_path.endswith('/'):
    spectra_path = spectra_path[:-1]

stop = False
if not os.path.exists(args.assembled):
    print(f"Error: {args.assembled} not found.")
    stop = True
for rawIn in args.raw:
    if not os.path.exists(rawIn):
        print(f"Error: {rawIn} not found.")
        stop = True
if stop:
    print("Missing input files. Please check these files.")
    exit()

with open(args.output, 'w') as f:
    f.write(
        "#!/bin/bash\n"
        "\n"
        "###### This code requires your path to have:\n"
        "### jellyfish2\n"
        "### python3\n"
        "### R4\n"
        "### Spectra + dependencies\n"
        "######\n\n"
        "###### Run raw jellyfish calculations, then dump and sort kmers above minimum.\n"
    )

    if len(args.raw)>1:
        for rawIn in args.raw:
            if os.path.exists(rawIn):
                f.write(f"{args.jf_path} count -t {args.threads} -s {args.jf_bloom} -m {args.mer_size} -o {args.prefix}_rcp_{os.path.basename(rawIn)}.jfc -C " + (f"<(zcat {rawIn})\n" if rawIn.endswith(".gz") else f"{rawIn}\n"))
                f.write(f"{args.jf_path} stats {args.prefix}_rcp_{os.path.basename(rawIn)}.jfc > {args.prefix}_rcp_{os.path.basename(rawIn)}.jstats")
            else:
                print(f"Error: file {rawIn} not found. Excluded from script.")
        f.write(f"{args.jf_path} merge -o {args.prefix}_raw_count.jfc {args.prefix}_rcp_*.jfc\n")
    else:
        f.write(f"{args.jf_path} count -t {args.threads} -s {args.jf_bloom} -m {args.mer_size} -o {args.prefix}_raw_count.jfc -C "+ (f"<(zcat {args.raw[0]})\n" if args.raw[0].endswith(".gz") else f"{args.raw[0]}\n"))
    f.write(f"{args.jf_path} stats {args.prefix}_raw_count.jfc > {args.prefix}_raw_count.jstats")
    f.write(f"{args.jf_path} dump -L {args.raw_min} -c {args.prefix}_raw_count.jfc |sort > {args.prefix}_raw.jdump\n")
    if args.clean:
        f.write(f"rm {args.prefix}_r*.jfc\n\n")

    f.write("###### Run assembly jellyfish calculations, then dump and sort kmers above minimum.\n")
    f.write(f"{args.jf_path} count -t {args.threads} -s {args.jf_bloom} -m {args.mer_size} -o {args.prefix}_asm_count.jfc -C {args.assembled}\n")
    f.write(f"{args.jf_path} stats {args.prefix}_asm_count.jfc > {args.prefix}_asm_count.jstats")
    f.write(f"{args.jf_path} dump -L {args.asm_min} -c {args.prefix}_asm_count.jfc |sort > {args.prefix}_asm.jdump\n")
    if args.clean:
        f.write(f"rm {args.prefix}_asm_count.jfc\n\n")

    f.write(f"mkdir {args.prefix}\n\n")

    f.write(f"###### Generate kmer report files and figures\n")
    f.write(f"{args.python} {spectra_path}/scripts/utils/kmerComp.py -r {args.prefix}_raw.jdump -a {args.prefix}_asm.jdump -k {args.mer_size} -o {args.prefix}/{args.prefix}_kmer_comp -s {kmer_sample_size} -p {percentile_keep} -v\n")
    f.write(f"{args.python} {spectra_path}/scripts/utils/kmerRank.py -r {args.prefix}_raw.jdump -a {args.prefix}_asm.jdump -o {args.prefix}_kmer_rank.tsv -c {chunk_size} -e {percentile_keep} -v\n\n")

    if args.clean:
        f.write(f"rm {args.prefix}_asm.jdump {args.prefix}_raw.jdump\n\n")

    f.write(f"###### Generate and plot localization of extreme kmers\n")
    f.write(f"{args.python} {spectra_path}/scripts/utils/mass-query.py -i {args.assembled} -q {args.prefix}_kmer_rank.tsv -m {args.mer_size} -o {args.prefix}_mass_query.tsv -c -w {args.mq_window} -s {args.mq_window} -v\n")
    f.write(f"{args.rscript} {spectra_path}/scripts/utils/mass-query-plot.r -i {args.prefix}_mass_query.tsv -o {args.prefix}/{args.prefix}_mass\n\n")

    f.write("###### Generate Spectra\n")
    f.write(f"{args.python} {spectra_path}/spectra.py count -w {args.spectra_window} -s {args.spectra_window} -i {args.assembled} -o {args.prefix}_spectra.tsv -v\n")
    f.write(f"{args.rscript} {spectra_path}/spectra-plot.r -i {args.prefix}_spectra.tsv -o {args.prefix}/{args.prefix}_spectra\n\n")

    f.write(f"###### Collate information into PDF report\n")
    f.write(f"{args.python} {spectra_path}/scripts/utils/pdfReport.py -i {args.prefix} -o {args.prefix}_report.pdf -m {args.mer_size} -p {args.prefix}\n")
