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

Software is managed with [pixi](https://pixi.sh). The `pixi.toml` defines a
default environment (used to drive the Snakemake workflows, plotting and
notebooks) plus one isolated environment per benchmarked tool. The Snakemake
rules activate the relevant per-tool environment themselves at runtime via
`pixi shell-hook`, so there is no longer any need for `--use-conda`.

To run a benchmark, first install the environments (this solves and downloads
all tool environments up front):

```bash
cd singlem-benchmarking
pixi install --all
```

You can either prefix commands with `pixi run` (as shown below), or enter the
default environment once with `pixi shell` and drop the prefix.

First, download the reference databases for each tool
```bash
pixi run snakemake --snakefile gather_tool_databases.smk -c 8
```

To generate a GTDB v207 database for sylph you will need a folder with all the GTDB v207 genomes (here ALL_GTDBR207_GENOMES_DIR)
```bash
#ALL_GTDBR207_GENOMES_DIR=/work/microbiome/db/gtdb/gtdb_release207/genomic_files_reps/gtdb_genomes_reps_r207
cd tool_reference_data && {
    find $ALL_GTDBR207_GENOMES_DIR | grep .fna > gtdb_all.txt &&
        pixi run --environment sylph sylph sketch -l gtdb_all.txt -t 50 -o gtdb_database
    cd ..
}
```

For MetaKSSD clone the repository to obtain shuf files:
```bash
git clone --branch v2.24 https://github.com/yhg926/MetaKSSD.git tool_reference_data/MetaKSSD-checkout
```

Generate database
```bash
pixi run --environment metakssd metakssd dist \
    -L tool_reference_data/MetaKSSD-checkout/shuf_files/L3K11.shuf \
    -o tool_reference_data/GTDBr207_genomes_L3K11_sketch \
    -l tool_reference_data/gtdb_all.txt
cd tool_reference_data && {
        pixi run --environment metakssd metakssd set \
        -P GTDBr207_genomes_L3K11_sketch > metakssdr207_genome_name.txt &&
        cat ../*_taxonomy_r207.tsv > metakssdr207_genome2taxonomy.tsv &&
        pixi run --environment metakssd perl MetaKSSD-checkout/scripts/genome_species_labeling.pl \
        metakssdr207_genome_name.txt metakssdr207_genome2taxonomy.tsv > metakssdr207_group_name.txt &&
        pixi run --environment metakssd metakssd set \
        -g metakssdr207_group_name.txt -o GTDBr207_genomes_L3K11_sketch_pan \
        GTDBr207_genomes_L3K11_sketch &&
        pixi run --environment metakssd metakssd set -q \
        -o GTDBr207_genomes_L3K11_sketch_pan_union_sp \
        GTDBr207_genomes_L3K11_sketch_pan &&
        pixi run --environment metakssd metakssd set \
        -i GTDBr207_genomes_L3K11_sketch_pan_union_sp \
        -o GTDBr207_genomes_L3K11_sketch_markerdb \
        GTDBr207_genomes_L3K11_sketch_pan
    cd ..
}
```

Generate id to taxonomy mapping
```bash
pixi run --environment metakssd perl \
    tool_reference_data/MetaKSSD-checkout/scripts/gtdb_psid_species2krona_taxonomy.pl \
        tool_reference_data/metakssdr207_group_name.txt \
        tool_reference_data/metakssdr207_genome2taxonomy.tsv \
        > tool_reference_data/gtdbr207_psid2krona_taxonomy.tsv
```

Then run the benchmarking, for instance #1

```bash
cd 1_novel_strains
./run_benchmark.sh
```

Results can be viewed by rerunning the `plot.ipynb` in each benchmark directory, and then the `plot_overall.ipynb` notebook in the base directory.

To run the test benchmark 5 use

```bash
snakemake --snakefile run_benchmarks.smk -c 8 bench5
```

To run the test benchmark 7 use

```bash
snakemake --snakefile run_benchmarks.smk -c 8 bench7
```

or just download with

```bash
snakemake --snakefile run_benchmarks.smk -c 8 download_bench7
```

## Download genomes for benchmark #2

Using the NCBI datasets CLI (provided by the `ncbi-datasets-cli` package in the
default pixi environment). Either run these inside `pixi shell`, or prefix the
`datasets` calls with `pixi run`.

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
