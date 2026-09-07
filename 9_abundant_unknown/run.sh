#!/bin/bash -e
# Run benchmark 9 (10 genomes, 9 known in r207 plus one abundant genome new in
# r214) on the aqua queue via the snakemake aqua profile. Each rule is submitted
# as its own PBS job (snakemake_mqsub).
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking
PYTHONPATH=.. pixi run snakemake \
  --directory 9_abundant_unknown \
  -s 9_abundant_unknown/Snakefile \
  --profile aqua \
  --configfile 9_abundant_unknown/config-8threads.yaml \
  --keep-going "$@"
