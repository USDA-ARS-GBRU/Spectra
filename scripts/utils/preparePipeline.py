#!/usr/bin/env python3
import argparse
import os
import sys
import shlex
import subprocess

def get_memory_limit():
    """Attempt to detect total system memory in GB."""
    try:
        # Linux/Unix using os.sysconf
        return (os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')) // (1024**3)
    except (ValueError, AttributeError, OSError):
        try:
            # Fallback for macOS if sysconf fails
            return int(subprocess.check_output(['sysctl', '-n', 'hw.memsize']).strip()) // (1024**3)
        except (subprocess.CalledProcessError, FileNotFoundError, ValueError):
            return 4  # Default to 4GB if detection fails

# CLI arguments
def main():
    parser = argparse.ArgumentParser(description="Prepare inputs for specta pipeline, and generate a bash script for running.")
    parser.add_argument('-r', '--raw', dest='raw', required=True, nargs='+',help='Input raw fasta/fastq read file(s). These can be gzipped, but must end in ".gz". If multiple, separate with spaces')
    parser.add_argument('-a', '--assembled', dest='assembled', required=True, help='Input fasta/bgzipped fasta sequence assembly file')
    parser.add_argument('-o', '--output-script', dest='output', default='spectra-pipeline.sh', help='Output bash file')
    parser.add_argument('-p', '--output-prefix', dest='prefix', default='spectra_pipeline', help='Output files prefix. A directory will be created with this name for storing images')
    parser.add_argument('-t', '--threads', dest='threads', type=int, help='Processing threads for Jellyfish kmer counting', required=True)
    parser.add_argument('-k', '--kmer-size', dest='mer_size', type=int, help='kmer size in query [default 20]', default=20)
    parser.add_argument('-m', '--minimum-sequence-size', dest='minimum_size', type=int, help='Minimum sequence size to include in reports [100,000 bp]', default=100000)
    parser.add_argument('-c', '--counter', dest='counter', choices=['jellyfish', 'meryl'], default='jellyfish', help='K-mer counter to use [default jellyfish]')
    parser.add_argument('--canonical', dest='canonical', action='store_true', help='Generate canonical spectra in addition to non-canonical', default=False)
    parser.add_argument('--n-gaps', dest='ngaps', action='store_true', help='Label gaps in the assembly in the final report', default=False)
    parser.add_argument('--bin-identify', dest='bins', action='store_true', help='Label bin regions in the genome assembly', default=False)
    parser.add_argument('--bin-penalty', dest='bin_penalty', type=int, default=1000000, help='Penalty for bin identification [default 1,000,000]')
    parser.add_argument('--bin-size', dest='bin_size', type=int, default=5, help='Minimum size for bin identification [default 5]')
    parser.add_argument('--jellyfish-bloom', dest='jf_bloom', type=str, default='100M', help='Jellyfish2 count bloomfilter initial size [default 100M]')
    parser.add_argument('--jellyfish-path', dest='jf_path', type=str, default='jellyfish', help='Jellyfish2 path. Default assumes it is in your env [default jellyfish]')
    parser.add_argument('--jellyfish-disk', dest='jf_disk', action='store_true', default=False, help='Use Jellyfish2 count disk parameter for large raw data files [default False]')
    parser.add_argument('--jellyfish-count-sep', dest='jf_sep', action='store_true', default=False, help='Process multiple raw inputs separately before joining together instead of as one count [default False]')
    parser.add_argument('--meryl-memory', dest='meryl_memory', type=str, default=None, help='Meryl memory parameter (e.g. 16G). If not set, it will be automatically detected')
    parser.add_argument('--meryl-path', dest='meryl_path', type=str, default='meryl', help='Meryl path. Default assumes it is in your env [default meryl]')
    parser.add_argument('--hard-links', dest='hard_links', action='store_true', default=False, help='Hard-link input files into the current directory for processing [default False]')
    parser.add_argument('--python-callable', dest='python', type=str, default='python', help='python3 path. Default assumes it is in your env [default python]')
    parser.add_argument('--spectra-callable', dest='spectra', type=str, default=None, help='Spectra path. If not set, automatically detected from this script')
    parser.add_argument('--rscript-callable', dest='rscript', type=str, default='Rscript', help='Rscript path. Default assumes it is in your env [default Rscript]')
    parser.add_argument('--time', dest='time', action='store_true', default=False, help='Write timestamps for program progress [default False]')
    parser.add_argument('--sample-size', dest='sample_size', type=int, default=5000000, help='Number of randomly sampled k-mers to show in comparison plots.[default 5,000,000]')
    parser.add_argument('--chunk-size', dest='chunk_size', type=int, default=5000000, help='Maximum size of sequences to process on. Larger sequences will be segmented before processing.[default 5,000,000 bp]')
    parser.add_argument('--percentile-low', dest='percentile_low', type=float, default=5, help='Bottom N percent of kmers to keep [default 5]')
    parser.add_argument('--percentile-high', dest='percentile_high', type=float, default=5, help='Top N percent of kmers to keep [default 5]')
    parser.add_argument('--auto-percentile', dest='auto_percentile', action='store_true', default=False, help='Automatically determine low/high percentiles based on distribution')
    parser.add_argument('--raw-min', dest='raw_min', type=int, default=100, help='Jellyfish2 raw kmer minimum count to retain [default 100]')
    parser.add_argument('--asm-min', dest='asm_min', type=int, default=2, help='Jellyfish2 assembly kmer minimum count to retain [default 2]')
    parser.add_argument('--mq-window', dest='mq_window', type=int, default=200000, help='Window and spacing width for kmer mass-query.py localization [default 200,000 bp]')
    parser.add_argument('--spectra-window', dest='spectra_window', type=int, default=10000, help='Window and spacing width for spectra.py K=3 localization [default 10,000 bp]')
    parser.add_argument('--max-output', dest='max_output', type=int, default=50, help='Maximum number of contigs to include in final report to limit filesize [default 50]')
    parser.add_argument('--keep', dest='keep', action='store_false', help='Clean workspace as files are processed. Jellyfish kmer counts are very large. By default, these files are removed after processing.', default=True)
    parser.add_argument('--variable-paths', dest='variable', action='store_true', help='Code will use variables for naming of analysis files. Default is hard paths.', default=False)
    args = parser.parse_args()

    if args.counter == 'meryl' and args.meryl_memory is None:
        args.meryl_memory = f"{get_memory_limit()}G"

    spectra_path = args.spectra if args.spectra else os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
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
        sys.exit(1)

    if args.assembled.lower().endswith(("gz","gzip")):
        print(f"WARNING: Assembly file ending in 'gz/gzip' detected. Non-BGZF compression will cause crashing, see https://biopython.org/docs/latest/Tutorial/chapter_seqio.html#sec-seqio-index-bgzf for details.")

    # Begin writing the script file
    with open(args.output, 'w') as f:
        f.write(
            "#!/bin/bash\n"
            "set -e\n"
            "set -o pipefail\n\n"
            "###### This code requires your path to have:\n"
            f"### {args.counter} (or counter of choice)\n"
            "### python3\n"
            "### R4\n"
            "### Spectra + dependencies\n"
            "######\n\n"
        )

        f.write("##### Code generated using the command:\n")
        f.write(f"# {' '.join([shlex.quote(arg) for arg in sys.argv])}\n")
        f.write("#####\n\n")

        # If variables required, define variables from argument parser.
        variable_names = ["output", "prefix", "threads", "mer_size", "minimum_size", "jf_bloom", "jf_path", "python", "rscript", "sample_size", "chunk_size", "raw_min", "asm_min", "mq_window", "spectra_window", "assembled", "meryl_path", "meryl_memory", "counter", "max_output", "bin_penalty", "bin_size"]

        if args.variable:
            variables = {name: f'"${{{name}}}"' for name in variable_names}
            f.write("##### Naming variables to be used in analysis.\n")
            for name in variable_names:
                f.write(f'{name}={shlex.quote(str(args.__dict__[name]))}\n')
            for i in range(len(args.raw)):
                f.write(f'raw_{i}={shlex.quote(args.raw[i])}\n')
            variables['raw'] = [(f'"${{raw_{i}}}"', args.raw[i].lower().endswith((".gz", ".gzip"))) for i in range(len(args.raw))]
            f.write("#####\n\n")
        else:
            variables = {name: shlex.quote(str(args.__dict__[name])) for name in variable_names}
            variables['raw'] = [(shlex.quote(i), i.lower().endswith(("gz", "gzip"))) for i in args.raw]

        f.write("##### Image output directory.\n")
        f.write(f"mkdir -p {variables['prefix']}\n\n")

        if args.hard_links:
            f.write("##### Hard-linking input files to current directory for processing.\n")
            new_raw = []
            for i, (raw_path, is_gz) in enumerate(variables['raw']):
                local_raw = shlex.quote(f"local_raw_{i}.fasta" + (".gz" if is_gz else ""))
                f.write(f"ln -f {raw_path} {local_raw}\n")
                new_raw.append((local_raw, is_gz))
            variables['raw'] = new_raw

            local_asm = shlex.quote("local_assembled.fasta")
            f.write(f"ln -f {variables['assembled']} {local_asm}\n")
            variables['assembled'] = local_asm
            f.write("\n")

        if args.counter == 'jellyfish':
            # Begin writing raw jellyfish code
            f.write("###### Run raw jellyfish calculations, then dump and sort kmers above minimum.\n")
            if args.time:
                f.write(f'echo "Starting {variables["mer_size"]}-mer processing on raw data at:"\ndate\n')

            if len(variables['raw']) > 1:
                if args.jf_sep:
                    for i in range(len(variables['raw'])):
                        prefix_rcp = f"{variables['prefix']}_rcp_{shlex.quote(os.path.basename(str(args.raw[i])))}"
                        f.write(f"{variables['jf_path']} count {'--disk ' if args.jf_disk else ''}-t {variables['threads']} -s {variables['jf_bloom']} -m {variables['mer_size']} -o {prefix_rcp}.kcount -C " + (f"<(zcat {variables['raw'][i][0]})\n" if variables['raw'][i][1] else f"{variables['raw'][i][0]}\n"))
                        f.write(f"{variables['jf_path']} stats {prefix_rcp}.kcount > {prefix_rcp}.kstats\n")
                    f.write(f"{variables['jf_path']} merge -o {variables['prefix']}_raw_count.kcount {variables['prefix']}_rcp_*.kcount\n")
                else:
                    f.write(f"{variables['jf_path']} count {'--disk ' if args.jf_disk else ''}-t {variables['threads']} -s {variables['jf_bloom']} -m {variables['mer_size']} -o {variables['prefix']}_raw_count.kcount -C " + ' '.join([f"<({'z' if i[1] else ''}cat {i[0]})" for i in variables['raw']]) + '\n')
            else:
                f.write(f"{variables['jf_path']} count {'--disk ' if args.jf_disk else ''}-t {variables['threads']} -s {variables['jf_bloom']} -m {variables['mer_size']} -o {variables['prefix']}_raw_count.kcount -C "+ (f"<(zcat {variables['raw'][0][0]})\n" if variables['raw'][0][1] else f"{variables['raw'][0][0]}\n"))

            f.write(f"{variables['jf_path']} stats {variables['prefix']}_raw_count.kcount > {variables['prefix']}_raw_count.kstats\n")
            f.write(f"{variables['jf_path']} histo {variables['prefix']}_raw_count.kcount > {variables['prefix']}_raw_count.khisto\n")
            f.write(f"{variables['jf_path']} dump -L {variables['raw_min']} -c {variables['prefix']}_raw_count.kcount |sort > {variables['prefix']}_raw.kdump\n")

            if args.time:
                f.write(f'echo "Ending {variables["mer_size"]}-mer processing on raw data at:"\ndate\n\n')

            if args.keep:
                f.write(f"rm {variables['prefix']}_r*.kcount\n\n")
            else:
                f.write('\n')

            # Begin writing assembly jellyfish code
            f.write("###### Run assembly jellyfish calculations, then dump and sort kmers above minimum.\n")
            if args.time:
                f.write(f'echo "Starting {variables["mer_size"]}-mer processing on assembly data at:"\ndate\n')
            f.write(f"{variables['jf_path']} count -t {variables['threads']} -s {variables['jf_bloom']} -m {variables['mer_size']} -o {variables['prefix']}_asm_count.kcount -C {variables['assembled']}\n")
            f.write(f"{variables['jf_path']} stats {variables['prefix']}_asm_count.kcount > {variables['prefix']}_asm_count.kstats\n")
            f.write(f"{variables['jf_path']} histo {variables['prefix']}_asm_count.kcount > {variables['prefix']}_asm_count.khisto\n")
            f.write(f"{variables['jf_path']} dump -L {variables['asm_min']} -c {variables['prefix']}_asm_count.kcount |sort > {variables['prefix']}_asm.kdump\n")
            if args.time:
                f.write(f'echo "Ending {variables["mer_size"]}-mer processing on assembly data at:"\ndate\n\n')

            if args.keep:
                f.write(f"rm {variables['prefix']}_asm_count.kcount\n\n")
            else:
                f.write('\n')

        elif args.counter == 'meryl':
            # Begin writing raw meryl code
            f.write("###### Run raw meryl calculations, then dump and sort kmers above minimum.\n")
            if args.time:
                f.write(f'echo "Starting {variables["mer_size"]}-mer processing on raw data at:"\ndate\n')

            f.write(f"{variables['meryl_path']} count k={variables['mer_size']} memory={variables['meryl_memory']} threads={variables['threads']} output {variables['prefix']}_raw_count " + ' '.join([i[0] for i in variables['raw']]) + "\n")
            f.write(f"{variables['meryl_path']} statistics {variables['prefix']}_raw_count > {variables['prefix']}_raw_count.kstats\n")
            f.write(f"{variables['meryl_path']} histogram {variables['prefix']}_raw_count > {variables['prefix']}_raw_count.khisto\n")
            f.write(f"{variables['meryl_path']} print at-least {variables['raw_min']} {variables['prefix']}_raw_count | sort > {variables['prefix']}_raw.kdump\n")

            if args.time:
                f.write(f'echo "Ending {variables["mer_size"]}-mer processing on raw data at:"\ndate\n\n')

            if args.keep:
                f.write(f"rm -rf {variables['prefix']}_raw_count\n\n")
            else:
                f.write('\n')

            # Begin writing assembly meryl code
            f.write("###### Run assembly meryl calculations, then dump and sort kmers above minimum.\n")
            if args.time:
                f.write(f'echo "Starting {variables["mer_size"]}-mer processing on assembly data at:"\ndate\n')
            f.write(f"{variables['meryl_path']} count k={variables['mer_size']} memory={variables['meryl_memory']} threads={variables['threads']} output {variables['prefix']}_asm_count {variables['assembled']}\n")
            f.write(f"{variables['meryl_path']} statistics {variables['prefix']}_asm_count > {variables['prefix']}_asm_count.kstats\n")
            f.write(f"{variables['meryl_path']} histogram {variables['prefix']}_asm_count > {variables['prefix']}_asm_count.khisto\n")
            f.write(f"{variables['meryl_path']} print at-least {variables['asm_min']} {variables['prefix']}_asm_count | sort > {variables['prefix']}_asm.kdump\n")
            if args.time:
                f.write(f'echo "Ending {variables["mer_size"]}-mer processing on assembly data at:"\ndate\n\n')

            if args.keep:
                f.write(f"rm -rf {variables['prefix']}_asm_count\n\n")
            else:
                f.write('\n')

        # Begin writing kmer comparison code
        f.write(f"###### Generate kmer comparison\n")
        if args.time:
            f.write(f"echo 'Starting k-mer comparison and ranking at:'\ndate\n")

        kmer_comp_cmd = f"{variables['python']} {shlex.quote(spectra_path + '/scripts/utils/kmerComp.py')} -r {variables['prefix']}_raw.kdump -a {variables['prefix']}_asm.kdump -k {variables['mer_size']} -o {variables['prefix']}/{variables['prefix']}_kmer_comp -s {variables['sample_size']} --percentile-low {args.percentile_low} --percentile-high {args.percentile_high}"
        if args.auto_percentile:
            kmer_comp_cmd += " --auto"
        kmer_comp_cmd += " -v"
        f.write(kmer_comp_cmd + "\n")

        if args.auto_percentile:
            f.write(f"if [ -f {variables['prefix']}/{variables['prefix']}_kmer_comp_percentiles.txt ]; then\n")
            f.write(f"    source {variables['prefix']}/{variables['prefix']}_kmer_comp_percentiles.txt\n")
            f.write(f"else\n")
            f.write(f"    PERCENTILE_LOW={args.percentile_low}\n")
            f.write(f"    PERCENTILE_HIGH={args.percentile_high}\n")
            f.write(f"fi\n")
        else:
            f.write(f"PERCENTILE_LOW={args.percentile_low}\n")
            f.write(f"PERCENTILE_HIGH={args.percentile_high}\n")

        f.write(f"{variables['python']} {shlex.quote(spectra_path + '/scripts/utils/kmerRank.py')} -r {variables['prefix']}_raw.kdump -a {variables['prefix']}_asm.kdump -o {variables['prefix']}_kmer_rank.tsv -c {variables['chunk_size']} -v\n")
        if args.time:
            f.write(f"echo 'Ending k-mer comparison and ranking at:'\ndate\n\n")

        if args.keep:
            f.write(f"rm {variables['prefix']}_asm.kdump {variables['prefix']}_raw.kdump\n\n")
        else:
            f.write('\n')

        # Begin writing localization code
        f.write(f"###### Generate and plot localization of extreme kmers\n")
        if args.time:
            f.write(f'echo "Starting {variables["mer_size"]}-mer localization at:"\ndate\n')

        mass_query_cmd = f"{variables['python']} {shlex.quote(spectra_path + '/scripts/utils/mass-query.py')} -i {variables['assembled']} -q {variables['prefix']}_kmer_rank.tsv -m {variables['mer_size']} -o {variables['prefix']}_mass_query.tsv -c -w {variables['mq_window']} -t {variables['threads']} -s {variables['mq_window']} --minimum-size {variables['minimum_size']} --percentile-low $PERCENTILE_LOW --percentile-high $PERCENTILE_HIGH"

        f.write(mass_query_cmd + "\n")
        f.write(f"{variables['python']} {shlex.quote(spectra_path + '/spectra.py')} plot -i {variables['prefix']}_mass_query.tsv -o {variables['prefix']}/{variables['prefix']}_mass -a\n")
        if args.time:
            f.write(f'echo "Ending {variables["mer_size"]}-mer localization at:"\ndate\n\n')
        else:
            f.write(f"\n")

        # Begin writing spectra 3-mer code
        f.write("###### Generate Spectra\n")
        if args.time:
            f.write(f"echo 'Starting 3-mer localization at:'\ndate\n")

        # Non-canonical spectra
        f.write(f"{variables['python']} {shlex.quote(spectra_path + '/spectra.py')} count -w {variables['spectra_window']} -s {variables['spectra_window']} -i {variables['assembled']} -o {variables['prefix']}_spectra_standard.tsv --minimum-size {variables['minimum_size']} -t {variables['threads']} -v\n")
        # Circular plot still uses R as it's not implemented in Python yet
        f.write(f"{variables['rscript']} {shlex.quote(spectra_path + '/spectra-plot.r')} -i {variables['prefix']}_spectra_standard.tsv -o {variables['prefix']}/{variables['prefix']}_standard_circular -c -a\n")

        spectra_string = f"{variables['python']} {shlex.quote(spectra_path + '/spectra.py')} plot -i {variables['prefix']}_spectra_standard.tsv -o {variables['prefix']}/{variables['prefix']}_spectra_standard"
        if args.bins:
            f.write(f"{variables['python']} {shlex.quote(spectra_path + '/spectra.py')} analyze -i {variables['prefix']}_spectra_standard.tsv -o {variables['prefix']}_spectra_standard -p {variables['bin_penalty']} -s {variables['bin_size']} -v\n")
            spectra_string += f" --gff-file={variables['prefix']}_spectra_standard_bins.gff --gff-tracks=bin-region"
        if args.ngaps:
            f.write(f"{variables['python']} {shlex.quote(spectra_path + '/scripts/utils/n-counter.py')} -i {variables['assembled']} -o {variables['prefix']}_ngaps.gff -v\n")
            spectra_string += f" --ngaps={variables['prefix']}_ngaps.gff"
        f.write(spectra_string + "\n")

        # Canonical spectra
        if args.canonical:
            f.write(f"{variables['python']} {shlex.quote(spectra_path + '/spectra.py')} count -c -w {variables['spectra_window']} -s {variables['spectra_window']} -i {variables['assembled']} -o {variables['prefix']}_spectra_canonical.tsv --minimum-size {variables['minimum_size']} -t {variables['threads']} -v\n")
            spectra_string_canon = f"{variables['python']} {shlex.quote(spectra_path + '/spectra.py')} plot -i {variables['prefix']}_spectra_canonical.tsv -o {variables['prefix']}/{variables['prefix']}_spectra_canonical"
            if args.bins:
                f.write(f"{variables['python']} {shlex.quote(spectra_path + '/spectra.py')} analyze -i {variables['prefix']}_spectra_canonical.tsv -o {variables['prefix']}_spectra_canonical -p {variables['bin_penalty']} -s {variables['bin_size']} -v\n")
                spectra_string_canon += f" --gff-file={variables['prefix']}_spectra_canonical_bins.gff --gff-tracks=bin-region"
            if args.ngaps:
                spectra_string_canon += f" --ngaps={variables['prefix']}_ngaps.gff"
            f.write(spectra_string_canon + "\n")

        if args.time:
            f.write(f"echo 'Ending 3-mer localization at:'\ndate\n\n")
        else:
            f.write(f"\n")

        # Begin writing PDF report code
        f.write(f"###### Collate information into PDF report\n")
        if args.time:
            f.write(f"echo 'Starting PDF report generation at:'\ndate\n")

        raw_files_str = ' '.join([i[0] for i in variables['raw']])
        pdf_report_cmd = (f"{variables['python']} {shlex.quote(spectra_path + '/scripts/utils/pdfReport.py')} "
                          f"-i {variables['prefix']} -o {variables['prefix']}_report.pdf -m {variables['mer_size']} "
                          f"-p {variables['prefix']} --max-output {variables['max_output']}{' -b' if args.bins else ''} "
                          f"{'--canonical ' if args.canonical else ''}"
                          f" --percentile-low $PERCENTILE_LOW --percentile-high $PERCENTILE_HIGH")

        pdf_report_cmd += f" -r {raw_files_str} -a {variables['assembled']} -c {variables['counter']} -s {shlex.quote(spectra_path)}"

        f.write(pdf_report_cmd + "\n")

        if args.time:
            f.write(f"echo 'Ending PDF report generation at:'\ndate\n")

if __name__ == "__main__":
    main()
