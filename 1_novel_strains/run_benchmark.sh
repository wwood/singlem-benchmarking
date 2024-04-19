#!/bin/bash

rm -rf benchmarks

# for i in `echo 32 8 1`; do rm -rf output* && snakemake -c 32 --use-conda --conda-frontend mamba --configfile config-${i}threads.yaml; done
for i in `echo 1`; do rm -rf output* && PYTHONPATH=.. snakemake -c 32 --use-conda --conda-frontend mamba --configfile config-${i}threads.yaml; done
