# Overview

This pipeline downloads and processes raw FASTQ files into TPM values from SRA using kallisto and tximport.

This pipeline is currently processing data from *Distinct Immune Cell Populations Define Response to Anti-PD-1 Monotherapy and Anti-PD-1/Anti-CTLA-4 Combined Therapy* (Gide2019), *Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma.* (Hugo2016)  and *Tumor and Microenvironment Evolution during Immunotherapy with Nivolumab.* (Riaz2017).

## Usage

To submit the jobs in Snakemake on biowulf, use the following

`sbatch --time=10:00:00 --mem=4g --partition=norm,ccr scripts/run-kallisto-biowulf.sh`
