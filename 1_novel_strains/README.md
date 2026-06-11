To run these benchmarks, after installing the pixi environments (see ../README.md)

First generate the reads and ground truth community profiles:
```
pixi run snakemake -c 8 generate_all_communities --configfile config-8threads.yaml
```

Then run the benchmarks:
```
./run_benchmark.sh
```
