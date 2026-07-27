#!/bin/bash -e
# Run benchmark 7 (r232 metapackage, benchmark 6's community) on the aqua queue via
# the snakemake aqua profile. Each rule is submitted as its own PBS job
# (snakemake_mqsub). Only singlem, singlem-regime3 and sylph are run.
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking
PYTHONPATH=.. pixi run snakemake \
  --directory 7_r232_metapackage \
  -s 7_r232_metapackage/Snakefile \
  --profile aqua \
  --configfile 7_r232_metapackage/config-8threads.yaml \
  --keep-going "$@"
