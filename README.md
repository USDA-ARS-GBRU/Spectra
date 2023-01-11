# Spectra
A software toolkit for the analysis of trinucleotide distributions.

## Introduction
Spectra was developed as a way to measure and visualize the location of highly repetitious elements along chromosomal sequences, as evidenced through trinucleotide distributions. The included tools allow for counting trinucleotide elements, transforming count results, and visualizing with Spectra graphs.

## Features
- Fast counting of trinucleotide elements across sequence windows broken down to any preferred width and any amount of spacing or overlap.
- Many supported input formats. If Biopython.SeqIO supports the format, it should be able to take the input.
- Spectra data can be transformed through normalizations from a global frequency, removal of conforming frequencies, or consolidation of multiple windows without needing to rerun the primary analysis.
- Data can be analyzed to detect changes in 3-mer frequencies across sequences.
- Spectra graphs can be generated in numerous formats supported by R-`ggplot2` or python3-`plotnine`
- Can be run entirely with python from count through visualization or with visualization in R.
