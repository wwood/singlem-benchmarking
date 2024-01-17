# Benchmarks of several metagenome community profilers

This repository contains Snakemake workflows for benchmarking tools which report microbial taxonomies and their relative abundance in microbial communities, using metagenome sequencing as input. It focuses particularly on comparing the [SingleM](https://github.com/wwood/singlem) microbial profiler to others, but can be adapted to new profilers so long as they can output GTDB R207-based taxonomy profiles.

The benchmarks are:

1. `1_novel_strains/` (i.e. 'known species benchmark') - benchmark profilers using communities simulated from genomes which have been assigned taxonomies at the species level in GTDB (the genomes chosen are _not_ representative genomes, however).
2. `2_phylogenetic_novelty/` - benchmark profilers on community profiles made up of a novel lineage and a known species, at equal abundance. This benchmark tests the ability of profilers to detect and classify new lineages.

To run a benchmark, first create a conda env

```bash
mamba env create -n singlem-benchmarking -f env.yml
```

Then activate it

```bash
conda activate singlem-benchmarking
```

Then run a benchmarking, for instance #1

```bash
cd 1_novel_strains
./run_benchmark.sh
```

Results can be viewed by rerunning the `plot.ipynb` in each benchmark directory, and then the `plot_overall.ipynb` notebook in the base directory.

## Download genomes for benchmark #2

Using NCBI datasets CLI (on conda as `ncbi-datasets-cli=14.29.0`)

```bash
cd 2_phylogenetic_novelty
cd genomes
datasets download genome accession --inputfile ../genome_accessions.txt
unzip ncbi_dataset.zip

# Rename files to simple names (e.g. GCA_000508305.1_genomic.fna)
parallel --col-sep "\t" cp {1} {2} :::: ../genome_ncbi_names.tsv

cd ../genome_pairs
datasets download genome accession --inputfile ../genome_pairs_accessions.txt
unzip ncbi_dataset.zip

# Rename files to simple names (e.g. GCA_000508305.1_genomic.fna)
parallel --col-sep "\t" cp {1} {2} :::: ../genome_pairs_ncbi_names.tsv
```
