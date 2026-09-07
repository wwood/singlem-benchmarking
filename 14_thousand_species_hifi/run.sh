#!/bin/bash -e
# Run benchmark 14 (benchmark 6's 1000 known species, all below 1x coverage, as PacBio
# HiFi reads simulated with Badread) on the aqua queue via the snakemake aqua profile.
# Each rule is submitted as its own PBS job (snakemake_mqsub).
#
# As for benchmark 13, the read simulation rule is the long pole: Badread is
# single-threaded per genome, so 1000 genomes takes hours.
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking
PYTHONPATH=.. pixi run snakemake \
  --directory 14_thousand_species_hifi \
  -s 14_thousand_species_hifi/Snakefile \
  --profile aqua \
  --configfile 14_thousand_species_hifi/config-8threads.yaml \
  --keep-going "$@"
