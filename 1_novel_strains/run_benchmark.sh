#!/bin/bash -e

# Note that this script runs everything using 1 thread, assuming only 1 is
# available (as in e.g. a PBS job that has only been assigned one). Accuracy of
# results can be obtained using higher numbers of threads and higher values for
# snakemake -c, but doing so will not give correct runtime results.

rm -rf benchmarks
rm -rf output*

snakemake -c 1 --use-conda --conda-frontend mamba --configfile config-benchmarking.yaml; done
