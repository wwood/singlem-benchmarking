#!/bin/bash -e
# Run benchmark 10 (four genomes: two genera, each with a known r207 species at 50x
# plus an r214-novel congener, at 100x in one genus and 10x in the other) on the
# aqua queue via the snakemake aqua profile. Each rule is submitted as its own PBS
# job (snakemake_mqsub).
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking
PYTHONPATH=.. pixi run snakemake \
  --directory 10_congeneric_novelty \
  -s 10_congeneric_novelty/Snakefile \
  --profile aqua \
  --configfile 10_congeneric_novelty/config-8threads.yaml \
  --keep-going "$@"
