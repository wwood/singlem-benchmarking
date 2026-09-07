#!/bin/bash -e
# Run benchmark 16 (benchmark 2's 120 novel-lineage + known-species pairs at 10x each,
# as PacBio HiFi reads simulated with Badread) on the aqua queue via the snakemake aqua
# profile. Each rule is submitted as its own PBS job (snakemake_mqsub).
#
# 120 samples x (read simulation + 3 tools) is a lot of jobs; --keep-going means one
# failed sample does not stop the rest.
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking
PYTHONPATH=.. pixi run snakemake \
  --directory 16_phylogenetic_novelty_hifi \
  -s 16_phylogenetic_novelty_hifi/Snakefile \
  --profile aqua \
  --configfile 16_phylogenetic_novelty_hifi/config-8threads.yaml \
  --keep-going "$@"
