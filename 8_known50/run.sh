#!/bin/bash -e
# Run benchmark 8 (12 genomes, half novel in r214, half Streptomyces) on the aqua
# queue via the snakemake aqua profile. Each rule is submitted as its own PBS job
# (snakemake_mqsub).
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking
PYTHONPATH=.. pixi run snakemake \
  --directory 8_known50 \
  -s 8_known50/Snakefile \
  --profile aqua \
  --configfile 8_known50/config-8threads.yaml \
  --keep-going "$@"
