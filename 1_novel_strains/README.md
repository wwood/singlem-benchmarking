To run these benchmarks, first setup the master environment

```
mamba env create -f env.yml -p ./env
conda activate ./env
```

Then run the benchmarks with different numbers of threads
```
./run_benchmark.sh
```
