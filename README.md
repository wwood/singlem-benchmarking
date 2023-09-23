# Somewhat repeatable benchmarking of several metagenome community profilers

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

Results can be viewed by rerunning the plot.ipynb in each benchmark directory.

## Download genomes for #2

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
