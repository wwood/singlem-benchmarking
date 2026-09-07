#!/bin/bash -e
# Run benchmark 12 (benchmark 11's two known species at 10x each, but as PacBio HiFi
# reads simulated with Badread) on the aqua queue via the snakemake aqua profile. Each
# rule is submitted as its own PBS job (snakemake_mqsub).
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking
PYTHONPATH=.. pixi run snakemake \
  --directory 12_longread_2component_hifi \
  -s 12_longread_2component_hifi/Snakefile \
  --profile aqua \
  --configfile 12_longread_2component_hifi/config-8threads.yaml \
  --keep-going "$@"
