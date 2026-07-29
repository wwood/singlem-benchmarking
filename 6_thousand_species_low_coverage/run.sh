#!/bin/bash -e
# Run benchmark 6 (1000 known species, all below 1x) on the aqua queue via the snakemake
# aqua profile. Each rule is submitted as its own PBS job (snakemake_mqsub).
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking
PYTHONPATH=.. pixi run snakemake \
  --directory 6_thousand_species_low_coverage \
  -s 6_thousand_species_low_coverage/Snakefile \
  --profile aqua \
  --configfile 6_thousand_species_low_coverage/config-8threads.yaml \
  --keep-going "$@"
