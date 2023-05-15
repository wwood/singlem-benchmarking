Somewhat repeatable benchmarking of several metagenome community profilers

To run a benchmark, first create a conda env

```
$ mamba env create -n singlem-benchmarking -f env.yml
```

Then activate it

```
$ conda activate singlem-benchmarking
```

Then run a benchmarking, for instance #1

```
$ cd 1_novel_strains
$ ./run_benchmark.sh
```

Results can be viewed by rerunning the plot.ipynb in each benchmark directory.
