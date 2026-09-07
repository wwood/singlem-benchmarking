#!/bin/bash -e
# Run benchmark 15 (benchmark 2's 120 novel-lineage + known-species pairs at 10x each,
# as Oxford Nanopore R10.4.1 reads simulated with Badread) on the aqua queue via the
# snakemake aqua profile. Each rule is submitted as its own PBS job (snakemake_mqsub).
#
# 120 samples x (read simulation + 3 tools) is a lot of jobs; --keep-going means one
# failed sample does not stop the rest.
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking
PYTHONPATH=.. pixi run snakemake \
  --directory 15_phylogenetic_novelty_nanopore \
  -s 15_phylogenetic_novelty_nanopore/Snakefile \
  --profile aqua \
  --configfile 15_phylogenetic_novelty_nanopore/config-8threads.yaml \
  --keep-going "$@"
