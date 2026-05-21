![Spectra logo](includes/Spectra-Logo.png)

Spectra: Software for analysis and visualization of 3-mers in genetic sequence data
==================================================================================================

## Introduction

Spectra is a software toolkit for the analysis of 3-mer (tri-nucleotide) distributions. Spectra was developed as a way
to assess regions of DNA that are comprised of unique sets of k-mers from one another, and to break those k-mers down
to their basal elements. This offers the ability to measure and localize shifts in genetic composition and identify the
relation between various tandem repetitive elements.

Spectra is run as a script through Python 3.x and R 4.x with a low overhead of required packages. It's primarily
developed and tested in WSL2 Debian.

## Installing with conda
For Linux:
```shell
cd Spectra
conda env create -n spectra -f environment.yml
conda activate spectra
```

## Installing standalone
For Linux:
```shell
pip install biopython numpy pandas plotly ruptures scipy matplotlib seaborn reportlab
Rscript -e "install.packages(c('ape','egg','ggplot2','dplyr','optparse','svglite','tidyr', 'readr', 'scales'), dependencies=TRUE)"
```

## Test data - [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/) human genome
These steps can be used to test the functionality of Spectra for counting and visualizing kmer frequencies
and pairing this with gene annotation data. Requires about 4gb of storage space and 4gb of RAM. This will produce the
use-case for larger genome sizes, so smaller genome sizes or less-contiguous assemblies will run faster. Run these
commands in a folder containing the `Spectra` folder.
```shell
# Download and unzip human genome data from genbank
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
gunzip GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
# Run 3-mer counting
python Spectra/spectra.py count -i GCF_009914755.1_T2T-CHM13v2.0_genomic.fna -o Hsap.tsv
# Run visualization
Rscript Spectra/spectra-plot.r -i Hsap.tsv
```

### Chromosome 1 (NC_060925.1)
![Spectra result](includes/example_data/Hsap.tsv_NC_060925.1.png)

## Comparative k-mer analysis pipeline
A k-mer profiling of raw/assembled data can be generated following the instructions on [PIPELINE.md](PIPELINE.md). This will run Spectra 3-mer counting and visualization as well as aggregate information about k-mer composition to a PDF report.

## Advanced usage

All Python functionalities are accessible through the `spectra.py` wrapper.

---
### Count
Basic usage: `python spectra.py count -i INPUT_SEQUENCE -o OUTPUT_TSV`

Generate absolute trinucleotide counts or canonical counts (`-c`) from input nucleotide sequence data
(`-i INPUT_SEQUENCE`) and output a tabular-separated values file (`-o OUTPUT_TSV`). Window size set with `-w WIDTH` and
spaced with `-s SPACING`. For non-overlapping windows, width and spacing must be the same. For data with multiple data
sources, this can be run together with `-l` and data sources will be defined by the first underscore in sequence
headers. Frequencies can optionally be reported instead of raw counts (`-p`). K-mer size can be changed with `-k` (default 3).

### Plot (R)
Basic usage: `Rscript spectra-plot.r -i INPUT_TSV`

Generate a plot for each sequence as a png. Output file prefix can be supplied with `-o output_prefix` and image format with `-f png/svg/tiff`. DPI resolution set with `-r 300`. Plot subset of libraries with `-n LIBRARY1,LIBRARY2,...`,
subset of sequences with `-s SEQUENCE1,SEQUENCE2,...`, or partial sequence range with `-w START,END`. Regular
expressions to filter libraries and sequences supported with `-e`. Legend shown with `-l`. Axis labels omitted with
`-a`. Plotting scale set with `-k` to keep absolute scale and `-x N` to set scale to N million base pairs every 1 inch. If data
already supplied as frequencies, use `-q`. Datasets with multiple libraries defaults to a faceted plot for each sequence
header shared between libraries.

External data sources can be supplied to plot in addition to spectra profile. GFF data can be plotted using
`--gff-file GFF_FILE --gff-tracks TRACK1,TRACK2,...`. Outputs from Tandem Repeat Finder can be supplied by converting the output
with `scripts/utils/trfWindows.py` and visualized with `--trf-file TRFFILE.tsv`. N-gaps can be highlighted with `--ngaps GFF_FILE`.

A circular plot can be produced with `-c` and the plotting length can be set with `--graphlength LENGTH` to scale alongside other
datasets.

### Interactive Plot (Python)
Basic usage: `python spectra.py plot -i INPUT_TSV -o OUTPUT_HTML`

Generates interactive [Plotly](https://plotly.com/python/) webpages for navigating across the sequence and selectively
viewing a subset of trinucleotides. It is not recommended to use this with small window sizes or dense spacing as these
are computationally costly to display in browser.

### Analyze
Basic usage: `python spectra.py analyze -i INPUT_TSV -o OUTPUT_PREFIX`

Detect breakpoints along sequences where spectra identity shifts using
[Ruptures](https://centre-borelli.github.io/ruptures-docs/). Reports the boundaries of each segment, and generates
a modified Spectra tsv and frequency profile of segments. Breakpoint penalty can be set with
`-p PENALTY`, which defaults to 1,000,000, and the minimum number of Spectra windows in a segment can be set with
`-s WINDOWS`, which defaults to 5. Breakpoints can be inferred from Spectra frequencies with `-f`.

### Transform
Basic usage: `python spectra.py transform -i INPUT_TSV -o OUTPUT_TSV`

Transform Spectra data tsv in a number of ways. Convert between counts and frequencies with `-c`. Windows can be increased by
whole number factors with `-s N` to summarize N windows into 1. Outlier frequencies from the genome-wide means can be
identified with `-r`. Global frequencies can be printed with `-p`. Simplify forward and reverse-complement counts with `-y`.

### Collate
Basic usage: `python spectra.py collate -i INPUT_TSV1 INPUT_TSV2 ... -o OUTPUT_TSV`

Multiple Spectra data tsvs can be collated together, such as individual runs of multiple genomes.

### Query
Basic usage: `python spectra.py query -i INPUT_SEQUENCE -q QUERY1,QUERY2,... -o OUTPUT_TSV`

Using the Spectra-count algorithm, query motifs can be supplied for counting and visualization. Parameterization
follows Spectra-count.
