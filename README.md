# Benchmarks of several metagenome community profilers

This repository contains Snakemake workflows for benchmarking tools which report microbial taxonomies and their relative abundance in microbial communities, using metagenome sequencing as input. It focuses particularly on comparing the [SingleM](https://github.com/wwood/singlem) microbial profiler to others, but can be adapted to new profilers so long as they can output GTDB R207-based taxonomy profiles.

The benchmarks are:

1. `1_novel_strains/` (i.e. 'known species benchmark') - benchmark profilers using communities simulated from genomes which have been assigned taxonomies at the species level in GTDB (the genomes chosen are _not_ representative genomes, however).
2. `2_phylogenetic_novelty/` - benchmark profilers on community profiles made up of a novel lineage and a known species, at equal abundance. This benchmark tests the ability of profilers to detect and classify new lineages.
3. `3_cami2_marine` - benchmark profilers on CAMI2 marine datasets, after converting the taxonomy to GTDB R207-based taxonomy.
4. `4_complex_and_novel` - benchmark profilers on a complex community (defined by the CAMI2 marine coverages), where 0-100% of the community is new in GTDB R214 compared to R207.

To get this repository, git clone with recursive option to get the submodules:

```bash
git clone --recursive https://github.com/wwood/singlem-benchmarking
```

To run a benchmark, first create a conda env

```bash
cd singlem-benchmarking
mamba env create -n singlem-benchmarking -f env.yml
```

Then activate it

```bash
conda activate singlem-benchmarking
```

First, download the reference databases for each tool
```bash
snakemake --snakefile gather_tool_databases.smk --use-conda -c 8
```

Then run the benchmarking, for instance #1

```bash

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
