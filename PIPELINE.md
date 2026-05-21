# Spectra Quick Start Pipeline

## Overview
The Spectra pipeline is designed to compare the genetic composition of raw sequencing reads (the "raw data") with a genome assembly (the "assembled data"). By analyzing short sequences of DNA, it helps assess assembly quality, identify shifts in sequence identity, and localize specific genetic elements.

### Key Concepts
* **K-mer**: A short sequence of DNA of length *k*. For example, a 20-mer is a sequence of 20 nucleotides.
* **3-mer (Tri-nucleotide)**: A k-mer of length 3 (e.g., AAA, ACG). Spectra uses 3-mer distributions to create a "fingerprint" of different genomic regions.
* **Abundance/Count**: The number of times a specific k-mer appears in a dataset.
* **Log-fold Change**: A mathematical way to compare abundance. A high positive value means a k-mer is much more common in the assembly than the raw reads, while a negative value means it is under-represented in the assembly.

## Prerequisites
To run this pipeline, you will need the following:
* Computing needs:
  * Linux is preferred, but Unix systems should be able to run this pipeline. Systems were tested in Debian WSL2 for
  Windows.
  * High-throughput computing access preferable, but local machines can run this entire pipeline. Storing and running
  k-mer analysis on large raw read files is ill-advised. All successive steps can be run easily.
* [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/getting-started.html) or a comparable conda
* [Jellyfish2](http://academic.oup.com/bioinformatics/article/27/6/764/234905) or [meryl](https://link.springer.com/article/10.1186/s13059-020-02134-9) if not using Conda
packager
* Raw HiFi or comparable long read data, formatted as FASTQ or FASTA. No quality-control is run on this read data,
so the exact filtered data that was used in the assembly is preferred.
* Assembled genome, formatted as a FASTA. Chromosome-length contig assemblies are preferable, but any assembly will
work. For sequence visualization and the report, highly fragmented genomes will generate a lot of image data and final
report will be truncated to keep PDF size small.
* Follow installation instructions for [README.md](README.md).

## Step 1: Preparing the pipeline
The pipeline is managed by a helper script, `preparePipeline.py`, which generates a customized bash script to run the various analysis steps.

1.  **Run `preparePipeline.py`**: Provide your raw reads, assembled genome, and the number of processing threads.

    *   **Required inputs**:
        *   `-r raw.fastq`: One or more raw read files (FASTA/FASTQ, can be gzipped).
        *   `-a assembled.fasta`: The assembled genome file.
        *   `-t 20`: Number of threads to use for k-mer counting.

    *   **Commonly used options**:
        *   `-k 20`: K-mer size (default: 20).
        *   `-o script.sh`: Name of the output bash script.
        *   `-p prefix`: Prefix for all output files and the image directory.
        *   `--bin-identify`: Identify and plot sequence pattern shifts.
        *   `--n-gaps`: Highlight gaps (N's) in the assembly report.

    *   **Example (using Jellyfish2)**:
        ```bash
        python scripts/utils/preparePipeline.py -r raw1.fq.gz raw2.fq.gz -a scaffolded.fa -o spectra_script.sh -p spectra_run -t 20
        ```

    *   **Example (using Meryl)**:
        ```bash
        python scripts/utils/preparePipeline.py -c meryl -r raw1.fq.gz raw2.fq.gz -a scaffolded.fa -o spectra_script.sh -p spectra_run -t 20 --meryl-memory 32G
        ```

2.  **Understand Default Settings**: If not specified, the pipeline uses:
    *   Sliding window sizes: 10,000bp for 3-mers, 200,000bp for larger k-mers.
    *   Minimum counts: Raw k-mers must appear at least 100 times, assembly k-mers at least 2 times.
    *   Intermediate files are automatically cleaned up to save space.

## Step 2: Running the pipeline
Once you have generated the bash script (e.g., `spectra_script.sh`), you can execute the entire analysis.

1.  **Run the script**:
    ```bash
    bash spectra_script.sh
    ```
2.  **Resource Management**: For large datasets, it is recommended to run this on a high-performance computing (HPC) cluster. You may need to modify the generated script to include scheduler headers (like SLURM `#SBATCH` directives).
3.  **Environment**: Ensure `jellyfish` (or `meryl`), `python`, and `Rscript` are in your system's PATH.

# Detailed Output Explanation

The pipeline generates several files and a comprehensive PDF report. Below is a detailed breakdown of what these outputs represent.

## Summary Files
*   `{prefix}_raw_count.jstats`: Summary statistics for k-mers found in the raw reads (e.g., total k-mers, unique k-mers).
*   `{prefix}_asm_count.jstats`: Summary statistics for k-mers found in the assembly.
*   `{prefix}_spectra.tsv`: A large table containing the counts of all 64 possible 3-mers across sliding windows of the assembly.
*   `{prefix}_kmer_rank.tsv`: A list of all k-mers (size *K*) ranked by their log-fold change between raw reads and the assembly.

## The PDF Report (`{prefix}_report.pdf`)
The report is the primary way to visualize the results. It is divided into several sections:

### 1. K-mer Distribution Plots
These plots help you see how well the k-mers in your assembly match the k-mers in your raw reads.
*   **K-mer Prevalence (Scatter Plot)**: Compares k-mer counts in raw data (x-axis) vs. assembly (y-axis). Ideally, most points should fall along a diagonal line. Points far from the diagonal indicate k-mers that are over- or under-represented in the assembly.
*   **Abundance Shift (Density Plot)**: Shows the distribution of log-fold changes. A peak at zero means most k-mers have similar representation in both datasets.
*   **Violin/Back-to-Back Density Plots**: Compare the overall "shape" of k-mer abundances in raw vs. assembled data.

### 2. Empirical Cumulative Distribution (ECDF)
This plot shows the cumulative probability of log-fold changes. It's a technical way to see what proportion of k-mers fall below a certain change threshold.

### 3. Sequence-Specific Breakdowns
For each major sequence (contig/chromosome) in your assembly, the report provides:
*   **High-Abundance K-mers (Top Plot)**: Localizes k-mers that are much more common in the assembly than expected (potential collapses or repetitive elements).
*   **Spectra (Middle Plot)**: Visualizes the 3-mer "fingerprint" along the sequence. Different colors represent different 3-mers. Sudden shifts in the color patterns can indicate boundaries between different types of genetic material (e.g., transitions into centromeres or telomeres).
*   **Low-Abundance K-mers (Bottom Plot)**: Localizes k-mers that are missing or under-represented in the assembly compared to the raw reads (potential mis-assemblies or missing sequence).

## Optional Outputs
*   `{prefix}_spectra_bins.gff`: If `--bin-identify` is used, this file contains the coordinates of regions where the 3-mer composition shifts significantly.
*   `{prefix}_ngaps.gff`: If `--n-gaps` is used, this file marks the locations of gaps (represented by 'N's) in the assembly.
