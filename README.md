# Overview

This pipeline downloads and processes raw FASTQ files into TPM values from SRA using kallisto and tximport.

This pipeline is currently processing data from *Distinct Immune Cell Populations Define Response to Anti-PD-1 Monotherapy and Anti-PD-1/Anti-CTLA-4 Combined Therapy*, *Genomic correlates of response to CTLA-4 blockade in metastatic melanoma* and *Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma*.

## Usage

To submit the jobs in Snakemake on biowulf, use the following

`sbatch --time=10:00:00 --mem=4g --partition=norm,ccr scripts/run-kallisto-biowulf.sh`
