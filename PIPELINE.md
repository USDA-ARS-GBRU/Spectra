# Spectra quick start pipeline
This pipeline is designed for quick use of Spectra 3-mer plotting and comparison of pre/post-assembly N-mer composition.

To run this pipeline, you will need the following:
* Computing needs:
  * Linux is preferred, but Unix systems should be able to run this pipeline. Systems were tested in Debian WSL2 for
  Windows
  * High-throughput computing access preferable,but local machines can run this entire pipeline. Storing and running
  k-mer analysis on large raw read files is ill-advised. All successive steps can be run easily.
* [Jellyfish2](http://academic.oup.com/bioinformatics/article/27/6/764/234905)
* [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/getting-started.html) or a comparable conda
packager
* Raw HiFi or comparable long read data, formatted as FASTQ or FASTA. No quality-control is run on this read data,
so the exact filtered data that was used in the assembly is preferred
* Assembled genome, formatted as a FASTA. Chromosome-length contig assemblies are preferable, but any assembly will
work. For sequence visualization and the report, highly fragmented genomes will generate a lot of image data and final
report will be truncated to keep PDF size small.
* Follow installation instructions for [README.md](README.md) as either a conda environment (and make sure the
environment is activated before running any scripts) or as a standalone. This will install the necessary Python 3.X and
R 4.X packages necessary for analyses.

# Preparing the pipeline
* Run [preparePipeline.py](scripts/utils/preparePipeline.py)
  * This will generate and annotate the scripts necessary for each step of Spectra.
  * Requires raw data (-r raw.fastq or -r raw1.fastq raw2.fastq ...) and assembled data (-a assembled.fasta)
  * If only raw and assembled data are supplied, scripts will use these default settings:
    * Comparative K-mer size of 20 (-k 20). Larger K-mer sizes will generate more K-mer information, but exponentially increase runtime.
    * Sliding window sizes 10000bp for 3-mer compoisiton (-w 10000), 200000bp for extreme k-mer localization (-m 200000). These windows are ideal for 200MBp-2GBp genomes. Larger genome sizes can use these settings, but localization information can be traded for final output sizes and runtimes with using larger windows.
    * Raw data k-mers with fewer than 100 hits (-R 100) and assembled data k-mers with fewer than 2 hits (-A 2) removed prior to analyses. Setting the raw k-mer threshold below the X-coverage of the genome is advised (Total Raw Data / Expected Genome Size = NX coverage. E.G, 50 GBp Raw / 1 Gbp Genome size = 50X coverage)
    * Intermediate K-mer files are removed during analyses. Files can optionally be kept (-c) for use in other k-mer profiling.
    * *(OPTIONAL)* Sequence pattern shifts can additionally be plotted (-b). N-gaps on assembled 3-mer profiles can be highlighted (-n).

# Running the pipeline
The `spectra-pipeline.sh` script can be run numerous ways. It can be run directly with `bash spectra-pipeline.sh`.
If running using shared computing resources, modification for submitting the job to a manager such as slurm is required.
Each script in `spectra-pipeline.sh` can be run individually. The script attempts to run all programs (jellyfish2, python, Rscript) from your PATH, and attempts to run all Spectra scripts from where the `preparePipeline.py` was ran.

# Outputted data
Spectra will supply files for further analyses, a folder of individual images, and a pdf report comprised of those images.

**Output files are:**
* spectra_pipeline_raw_count.jstats - N-mer count stats (Unique, Distinct, Total, Max) of raw data
* spectra_pipeline_asm_count.jstats - N-mer count stats (Unique, Distinct, Total, Max) of assembled data
* spectra_pipeline_spectra.tsv - 3-mer composition data
* spectra_pipeline/ - image folder
* spectra_pipeline_report.pdf - final output report

**Optional outputted files**
* spectra_pipeline_spectra_bins.tsv - predicted region shifts from 3-mer composition data in assembly
* spectra_pipeline_spectra_bins.gff - gff track of predicted regions in assembly
* spectra_pipeline_ngaps.gff - gff track of Ns in assembly