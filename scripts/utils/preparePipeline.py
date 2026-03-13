#!/usr/bin/env python3
import argparse
import os
import sys
# CLI arguments
parser = argparse.ArgumentParser(description="Prepare inputs for specta pipeline, and generate a bash script for running.")
parser.add_argument('-r', '--raw', dest='raw', required=True, nargs='+',help='Input raw fasta/fastq read file(s). These can be gzipped, but must end in ".gz". If multiple, separate with spaces')
parser.add_argument('-a', '--assembled', dest='assembled', required=True, help='Input fasta/bgzipped fasta sequence assembly file')
parser.add_argument('-o', '--output-script', dest='output', default='spectra-pipeline.sh', help='Output bash file')
parser.add_argument('-p', '--output-prefix', dest='prefix', default='spectra_pipeline', help='Output files prefix. A directory will be created with this name for storing images')
parser.add_argument('-t', '--threads', dest='threads', type=int, help='Processing threads for Jellyfish kmer counting', required=True)
parser.add_argument('-k', '--kmer-size', dest='mer_size', type=int, help='kmer size in query [default 20]', default=20)
parser.add_argument('-m', '--minimum-sequence-size', dest='minimum_size', type=int, help='Minimum sequence size to include in reports [100,000 bp]', default=100000)
parser.add_argument('--n-gaps', dest='ngaps', action='store_true', help='Label gaps in the assembly in the final report', default=False)
parser.add_argument('--bin-identify', dest='bins', action='store_true', help='Label bin regions in the genome assembly', default=False)
parser.add_argument('--jellyfish-bloom', dest='jf_bloom', type=str, default='100M', help='Jellyfish2 count bloomfilter initial size [default 100M]')
parser.add_argument('--jellyfish-path', dest='jf_path', type=str, default='jellyfish', help='Jellyfish2 path. Default assumes it is in your env [default jellyfish]')
parser.add_argument('--jellyfish-disk', dest='jf_disk', action='store_true', default=False, help='Use Jellyfish2 count disk parameter for large raw data files [default False]')
parser.add_argument('--python-callable', dest='python', type=str, default='python', help='python3 path. Default assumes it is in your env [default python]')
parser.add_argument('--spectra-callable', dest='spectra', type=str, default=None, help='Spectra path. If not set, automatically detected from this script')
parser.add_argument('--rscript-callable', dest='rscript', type=str, default='Rscript', help='Rscript path. Default assumes it is in your env [default Rscript]')
parser.add_argument('--time', dest='time', action='store_true', default=False, help='Write timestamps for program progress [default False]')
parser.add_argument('--sample-size', dest='sample_size', type=int, default=5000000, help='Number of randomly sampled k-mers to show in comparison plots.[default 5,000,000]')
parser.add_argument('--chunk-size', dest='chunk_size', type=int, default=5000000, help='Maximum size of sequences to process on. Larger sequences will be segmented before processing.[default 5,000,000 bp]')
parser.add_argument('--percentile', dest='percentile', type=int, default=5, help='Percentiles of highest/lowest to plot.[default 5]')
parser.add_argument('--raw-min', dest='raw_min', type=int, default=100, help='Jellyfish2 raw kmer minimum count to retain [default 100]')
parser.add_argument('--asm-min', dest='asm_min', type=int, default=2, help='Jellyfish2 assembly kmer minimum count to retain [default 2]')
parser.add_argument('--mq-window', dest='mq_window', type=int, default=200000, help='Window and spacing width for kmer mass-query.py localization [default 200,000 bp]')
parser.add_argument('--spectra-window', dest='spectra_window', type=int, default=10000, help='Window and spacing width for spectra.py K=3 localization [default 10,000 bp]')
parser.add_argument('--keep', dest='keep', action='store_false', help='Clean workspace as files are processed. Jellyfish kmer counts are very large. By default, these files are removed after processing.', default=True)
parser.add_argument('--variable-paths', dest='variable', action='store_true', help='Code will use variables for naming of analysis files. Default is hard paths.', default=False)
args = parser.parse_args()

spectra_path = args.spectra if args.spectra else os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
if spectra_path.endswith('/'):
    spectra_path = spectra_path[:-1]

# Check if input files exist, then terminate if any do not.
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

if args.assembled.endswith("gz") or args.assembled.endswith("gzip"):
    print(f"WARNING: Assembly file ending in 'gz/gzip' detected. Non-BGZF compression will cause crashing, see https://biopython.org/docs/latest/Tutorial/chapter_seqio.html#sec-seqio-index-bgzf for details.")

# Begin writing the script file
with open(args.output, 'w') as f:
    f.write(
        "#!/bin/bash\n\n"
        "###### This code requires your path to have:\n"
        "### jellyfish2\n"
        "### python3\n"
        "### R4\n"
        "### Spectra + dependencies\n"
        "######\n\n"
    )

    f.write("##### Code generated using the command:\n")
    f.write(f"# {' '.join(sys.argv)}\n")
    f.write("#####\n\n")

    # If variables required, define variables from argument parser. Then define whether argument or variable will be called.
    # This is redundant for arguments that are not changed, but done for consistency during development. These are taken directly from the CLI.
    # Rationale for redefining is that this is more flexible to keep track of downstream. Variables are in the dict variables, boolean settings are in args.
    variable_names = ["output", "prefix", "threads", "mer_size", "minimum_size", "jf_bloom", "jf_path", "python", "rscript", "sample_size", "chunk_size", "percentile", "raw_min", "asm_min", "mq_window", "spectra_window", "assembled"]

    if args.variable:
        variables = {name: "${" + f"{name}" + "}" for name in variable_names}
        f.write("##### Naming variables to be used in analysis.\n")
        for name in variable_names:
            f.write(f'{name}="{args.__dict__[name]}"\n')
        for i in range(len(args.raw)):
            f.write(f'raw_{i}="{args.raw[i]}"\n')
        variables['raw'] = [("${raw_" + f"{i}" + "}", any([args.raw[i].endswith("gz"), args.raw[i].endswith("gzip")])) for i in range(len(args.raw))]
        f.write("#####\n\n")
    else:
        variables = {name: args.__dict__[name] for name in variable_names}
        variables['raw'] = [(i, any([i.endswith("gz"), i.endswith("gzip")])) for i in args.raw]
    print(variables['raw'])
    f.write("##### Image output directory.\n")
    f.write(f"mkdir {variables['prefix']}\n\n")

    # Begin writing raw jellyfish code
    f.write("###### Run raw jellyfish calculations, then dump and sort kmers above minimum.\n")
    if args.time:
        f.write(f"echo 'Starting {variables['mer_size']}-mer processing on raw data at:'\ndate\n")
    if len(variables['raw'])>1:
        for i in range(len(variables['raw'])):
            f.write(f"{variables['jf_path']} count {'--disk ' if args.jf_disk else ''}-t {variables['threads']} -s {variables['jf_bloom']} -m {variables['mer_size']} -o {variables['prefix']}_rcp_{os.path.basename(variables['raw'][i][0])}.jfc -C " + (f"<(zcat {variables['raw'][i][0]})\n" if variables['raw'][i][1] else f"{variables['raw'][i][0]}\n"))
            f.write(f"{variables['jf_path']} stats {variables['prefix']}_rcp_{os.path.basename(variables['raw'][i][0])}.jfc > {variables['prefix']}_rcp_{os.path.basename(variables['raw'][i][0])}.jstats\n")
        f.write(f"{variables['jf_path']} merge -o {variables['prefix']}_raw_count.jfc {variables['prefix']}_rcp_*.jfc\n")
    else:
        f.write(f"{variables['jf_path']} count {'--disk ' if args.jf_disk else ''}-t {variables['threads']} -s {variables['jf_bloom']} -m {variables['mer_size']} -o {variables['prefix']}_raw_count.jfc -C "+ (f"<(zcat {variables['raw'][0][0]})\n" if variables['raw'][0][1] else f"{variables['raw'][0][0]}\n"))
    f.write(f"{variables['jf_path']} stats {variables['prefix']}_raw_count.jfc > {variables['prefix']}_raw_count.jstats\n")
    f.write(f"{variables['jf_path']} histo {variables['prefix']}_raw_count.jfc > {variables['prefix']}_raw_count.jhisto\n")
    f.write(f"{variables['jf_path']} dump -L {variables['raw_min']} -c {variables['prefix']}_raw_count.jfc |sort > {variables['prefix']}_raw.jdump\n")
    if args.time:
        f.write(f"echo 'Ending {variables['mer_size']}-mer processing on raw data at:'\ndate\n\n")

    if args.keep:
        f.write(f"rm {variables['prefix']}_r*.jfc\n\n")
    else:
        f.write('\n')

    # Begin writing assembly jellyfish code
    f.write("###### Run assembly jellyfish calculations, then dump and sort kmers above minimum.\n")
    if args.time:
        f.write(f"echo 'Starting {variables['mer_size']}-mer processing on assembly data at:'\ndate\n")
    f.write(f"{variables['jf_path']} count -t {variables['threads']} -s {variables['jf_bloom']} -m {variables['mer_size']} -o {variables['prefix']}_asm_count.jfc -C {variables['assembled']}\n")
    f.write(f"{variables['jf_path']} stats {variables['prefix']}_asm_count.jfc > {variables['prefix']}_asm_count.jstats\n")
    f.write(f"{variables['jf_path']} histo {variables['prefix']}_asm_count.jfc > {variables['prefix']}_asm_count.jhisto\n")
    f.write(f"{variables['jf_path']} dump -L {variables['asm_min']} -c {variables['prefix']}_asm_count.jfc |sort > {variables['prefix']}_asm.jdump\n")
    if args.time:
        f.write(f"echo 'Ending {variables['mer_size']}-mer processing on assembly data at:'\ndate\n\n")

    if args.keep:
        f.write(f"rm {variables['prefix']}_asm_count.jfc\n\n")
    else:
        f.write('\n')

    # Begin writing kmer comparison code
    f.write(f"###### Generate kmer comparison\n")
    if args.time:
        f.write(f"echo 'Starting k-mer comparison and ranking at:'\ndate\n")
    f.write(f"{variables['python']} {spectra_path}/scripts/utils/kmerComp.py -r {variables['prefix']}_raw.jdump -a {variables['prefix']}_asm.jdump -k {variables['mer_size']} -o {variables['prefix']}/{variables['prefix']}_kmer_comp -s {variables['sample_size']} -p {variables['percentile']} -v\n")
    f.write(f"{variables['python']} {spectra_path}/scripts/utils/kmerRank.py -r {variables['prefix']}_raw.jdump -a {variables['prefix']}_asm.jdump -o {variables['prefix']}_kmer_rank.tsv -c {variables['chunk_size']} -e {variables['percentile']} -v\n")
    if args.time:
        f.write(f"echo 'Ending k-mer comparison and ranking at:'\ndate\n\n")

    if args.keep:
        f.write(f"rm {variables['prefix']}_asm.jdump {variables['prefix']}_raw.jdump\n\n")
    else:
        f.write('\n')

    # Begin writing localization code
    f.write(f"###### Generate and plot localization of extreme kmers\n")
    if args.time:
        f.write(f"echo 'Starting {variables['mer_size']}-mer localization at:'\ndate\n")
    f.write(f"{variables['python']} {spectra_path}/scripts/utils/mass-query.py -i {variables['assembled']} -q {variables['prefix']}_kmer_rank.tsv -m {variables['mer_size']} -o {variables['prefix']}_mass_query.tsv -c -w {variables['mq_window']} -t {variables['threads']} -s {variables['mq_window']} --minimum-size {variables['minimum_size']} -v\n")
    f.write(f"{variables['rscript']} {spectra_path}/scripts/utils/mass-query-plot.r -i {variables['prefix']}_mass_query.tsv -o {variables['prefix']}/{variables['prefix']}_mass\n")
    if args.time:
        f.write(f"echo 'Ending {variables['mer_size']}-mer localization at:'\ndate\n\n")
    else:
        f.write(f"\n")

    # Begin writing spectra 3-mer code
    f.write("###### Generate Spectra\n")
    if args.time:
        f.write(f"echo 'Starting 3-mer localization at:'\ndate\n")
    f.write(f"{variables['python']} {spectra_path}/spectra.py count -w {variables['spectra_window']} -s {variables['spectra_window']} -i {variables['assembled']} -o {variables['prefix']}_spectra.tsv --minimum-size {variables['minimum_size']} -v\n")
    f.write(f"{variables['rscript']} {spectra_path}/spectra-plot.r -i {variables['prefix']}_spectra.tsv -o {variables['prefix']}/{variables['prefix']}_circular -c -a\n")
    spectraString=f"{variables['rscript']} {spectra_path}/spectra-plot.r -i {variables['prefix']}_spectra.tsv -o {variables['prefix']}/{variables['prefix']}_spectra"
    if args.bins:
        f.write(f"{variables['python']} {spectra_path}/spectra.py analyze -i {variables['prefix']}_spectra.tsv -o {variables['prefix']}_spectra -v\n")
        spectraString += f" -g {variables['prefix']}_spectra_bins.gff -t bin-region"
    if args.ngaps:
        f.write(f"{variables['python']} {spectra_path}/scripts/utils/n-counter.py -i {variables['assembled']} -o {variables['prefix']}_ngaps.gff -v\n")
        spectraString += f" -j {variables['prefix']}_ngaps.gff"
    f.write(spectraString+"\n")
    if args.time:
        f.write(f"echo 'Ending 3-mer localization at:'\ndate\n\n")
    else:
        f.write(f"\n")

    # Begin writing PDF report code
    f.write(f"###### Collate information into PDF report\n")
    if args.time:
        f.write(f"echo 'Starting PDF report generation at:'\ndate\n")
    f.write(f"{variables['python']} {spectra_path}/scripts/utils/pdfReport.py -i {variables['prefix']} -o {variables['prefix']}_report.pdf -m {variables['mer_size']} -p {variables['prefix']}{' -b' if args.bins else ''}\n")
    if args.time:
        f.write(f"echo 'Ending PDF report generation at:'\ndate\n")
