# Spectra quick start pipeline
This pipeline is designed for quick use of Spectra 3-mer plotting and comparison of pre/post-assembly N-mer composition.

To run this pipeline, you will need the following:
* Computing needs:
  * Linux is preferred, but Unix systems should be able to run this pipeline. Systems were tested in Debian WSL2 for
  Windows.
  * High-throughput computing access preferable, but local machines can run this entire pipeline. Storing and running
  k-mer analysis on large raw read files is ill-advised. All successive steps can be run easily.
* [Jellyfish2](http://academic.oup.com/bioinformatics/article/27/6/764/234905)
* [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/getting-started.html) or a comparable conda
packager
* Raw HiFi or comparable long read data, formatted as FASTQ or FASTA. No quality-control is run on this read data,
so the exact filtered data that was used in the assembly is preferred.
* Assembled genome, formatted as a FASTA. Chromosome-length contig assemblies are preferable, but any assembly will
work. For sequence visualization and the report, highly fragmented genomes will generate a lot of image data and final
report will be truncated to keep PDF size small.
* Follow installation instructions for [README.md](README.md).

# Preparing the pipeline
Run [preparePipeline.py](scripts/utils/preparePipeline.py)
* This will generate and annotate the scripts necessary for each step of Spectra.
* Requires raw data (`-r raw.fastq` or `-r raw1.fastq raw2.fastq ...`) and assembled data (`-a assembled.fasta`)
* If only raw and assembled data are supplied, scripts will use these default settings:
  * Comparative K-mer size of 20 (`-k 20`). Larger K-mer sizes will generate more K-mer information, but exponentially increase runtime.
  * Sliding window sizes 10,000bp for 3-mer composition (`--spectra-window 10000`), 200,000bp for extreme k-mer localization (`--mq-window 200000`).
  * Raw data k-mers with fewer than 100 hits (`--raw-min 100`) and assembled data k-mers with fewer than 2 hits (`--asm-min 2`) removed prior to analyses.
  * Intermediate K-mer files are removed during analyses.
  * *(OPTIONAL)* Sequence pattern shifts can additionally be plotted (`--bin-identify`). N-gaps on assembled 3-mer profiles can be highlighted (`--n-gaps`).
  * *(OPTIONAL)* The k-mer counter can be switched between Jellyfish2 and Meryl (`-c meryl`). Meryl memory can be specified with `--meryl-memory`.
  * *(OPTIONAL)* Input files can be hard-linked to the current directory to avoid path issues (`--hard-links`).
  * WARNING: When working with large raw data files, Jellyfish2 may require more memory than the resources provided. A flag (`--jellyfish-disk`) will require `jellyfish count` to use the `--disk` parameter.

Example call (Jellyfish2):
`python scripts/utils/preparePipeline.py -r raw1.fq.gz raw2.fq.gz -a scaffolded.fa -o spectra_script.sh -p spectra_run -t 20`

Example call (Meryl):
`python scripts/utils/preparePipeline.py -c meryl -r raw1.fq.gz raw2.fq.gz -a scaffolded.fa -o spectra_script.sh -p spectra_run -t 20 --meryl-memory 32G`

This will prepare the Jellyfish2, Python, and R scripts necessary for comparing the raw sequence data in `raw1.fq.gz` and `raw2.fq.gz` to the scaffolded genome assembly `scaffolded.fa`. The bash script with commands will be named `spectra_script.sh`, and all files produced will have the prefix `spectra_run`. Jellyfish2 will have 20 computing threads for k-mer counting.

# Running the pipeline
The generated `spectra-script.sh` script can be run directly with `bash spectra-script.sh`.
If running using shared computing resources, modification for submitting the job to a manager such as SLURM is required.
The script attempts to run all programs (jellyfish, python, Rscript) from your PATH.

# Outputted data
Spectra will supply files for further analyses, a folder of individual images, and a PDF report comprised of those images.

**Output files are:**
* `spectra_run_raw_count.jstats` - N-mer count stats of raw data
* `spectra_run_asm_count.jstats` - N-mer count stats of assembled data
* `spectra_run_spectra.tsv` - 3-mer composition data
* `spectra_run/` - image folder
* `spectra_run_report.pdf` - final output report

**Optional outputted files**
* `spectra_run_spectra_bins.tsv` - predicted region shifts from 3-mer composition data in assembly
* `spectra_run_spectra_bins.gff` - GFF track of predicted regions in assembly
* `spectra_run_ngaps.gff` - GFF track of Ns in assembly
