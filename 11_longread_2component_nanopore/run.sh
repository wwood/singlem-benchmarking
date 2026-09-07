#!/bin/bash -e
# Run benchmark 11 (two known species at 10x each, Oxford Nanopore R10.4.1 reads
# simulated with Badread) on the aqua queue via the snakemake aqua profile. Each rule
# is submitted as its own PBS job (snakemake_mqsub).
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking
PYTHONPATH=.. pixi run snakemake \
  --directory 11_longread_2component_nanopore \
  -s 11_longread_2component_nanopore/Snakefile \
  --profile aqua \
  --configfile 11_longread_2component_nanopore/config-8threads.yaml \
  --keep-going "$@"
