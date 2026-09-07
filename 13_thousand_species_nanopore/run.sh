#!/bin/bash -e
# Run benchmark 13 (benchmark 6's 1000 known species, all below 1x coverage, as Oxford
# Nanopore R10.4.1 reads simulated with Badread) on the aqua queue via the snakemake
# aqua profile. Each rule is submitted as its own PBS job (snakemake_mqsub).
#
# Note the read simulation rule is the long pole: Badread is single-threaded per genome
# (~0.4 Mbp/s), so 1000 genomes takes hours even at 8-way concurrency.
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking
PYTHONPATH=.. pixi run snakemake \
  --directory 13_thousand_species_nanopore \
  -s 13_thousand_species_nanopore/Snakefile \
  --profile aqua \
  --configfile 13_thousand_species_nanopore/config-8threads.yaml \
  --keep-going "$@"
