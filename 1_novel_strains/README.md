To run these benchmarks, after setting up the master conda environment (see ../README.md)

First generate the reads and ground truth community profiles:
```
snakemake --use-conda -c 8 generate_all_communities --configfile config-8threads.yaml
```

Then run the benchmarks:
```
./run_benchmark.sh
```
