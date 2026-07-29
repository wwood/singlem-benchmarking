#!/bin/bash -e
# Run benchmark 5 (5 known species @ 0.2x) on the aqua queue via the snakemake
# aqua profile. Each rule is submitted as its own PBS job (snakemake_mqsub).
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking
PYTHONPATH=.. pixi run snakemake \
  --directory 5_known_species_low_coverage \
  -s 5_known_species_low_coverage/Snakefile \
  --profile aqua \
  --configfile 5_known_species_low_coverage/config-8threads.yaml \
  --keep-going "$@"
