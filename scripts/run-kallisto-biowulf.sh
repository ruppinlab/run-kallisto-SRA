#!/bin/bash
module load snakemake || exit 1

snakemake --use-conda --cluster "sbatch --partition=norm,ccr --time=08:00:00 --mem=16g --cpus-per-task=4" --jobs 100 --latency-wait 60 --keep-going --local-cores 4 all
