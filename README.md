# Overview

This pipeline downloads and processes the RNA-seq data from "Distinct Immune Cell Populations Define Response to Anti-PD-1 Monotherapy and Anti-PD-1/Anti-CTLA-4 Combined Therapy"

## Usage

To submit the jobs in Snakemake on biowulf, use the following

`sbatch --time=10:00:00 --mem=4g --partition=norm,ccr scripts/run-kallisto-biowulf.sh`
