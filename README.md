# Variant Calling and Annotation

This script performs variant calling and annotation from GISAID fasta files.

## Usage

The script takes the following command-line arguments:

```shell
python variant_calling.py [-h] [-d DATADIR] [-r REFGENOMEDIR] [-c CPU] [-v]
```
`-d` or `--datadir`: Path for GISAID fasta files. The default directory is ./Data, which is the same directory where the script is located or running from.

`-r` or `--refGenomeDir`: Path for the reference genome directory. Use genomic fna sequence. The default directory is ./referenceGenome.

`-c` or `--cpu`: Maximum number of parallel CPU workers to use for multithreading.

`-v` or `--version`: Show the version of the script.


## Prerequisites
Before running the script, make sure you have the following dependencies installed:

- python3
- samtools
- minimap2
- freebayes
- snpEff
