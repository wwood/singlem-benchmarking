# Benchmark 5 — five known species at 0.2× coverage

A small, controlled synthetic community that probes the **coverage floor** seen
in the CAMI2 marine comparison: 0.2× per-genome coverage is well below singlem's
marker-gene sensitivity (~1.5× on 150 bp / ~5 Gbp samples), so this benchmark
isolates whether the whole-genome methods (sylph, and the `singlem-regime3` joint
method) still recover the community where vanilla singlem should struggle.

## The community

Five **known** GTDB r207 species (present in the tool databases), each simulated
at exactly **0.2× coverage** with `art_illumina` (150 bp paired, HSXt error
model). Equal coverage → ~20% relative abundance each.

| genome | species | genome size |
|---|---|---|
| GCF_001953075.1 | *Halioglobus lutimaris* | 5.0 Mbp |
| GCF_003815695.1 | *Chryseobacterium carnipullorum* | 5.5 Mbp |
| GCF_013321855.1 | *Rhizobium rhizogenes* | 7.3 Mbp |
| GCF_003326775.1 | *Thalassospira lucentensis_A* | 4.8 Mbp |
| GCF_000285335.1 | *Komagataeibacter europaeus* | 4.0 Mbp |

FASTAs come from `../reference_genomes/shadow/`; taxonomy and genome size are
looked up in `../bac120_metadata_r207.tsv`.

## Layout

- `genome_list.tsv` — the 5 genomes and their FASTA paths.
- `coverage_definitions/known5.tsv` — per-genome coverage (all 0.2×). The sample
  `known5` maps to `coverage_definitions/known5.tsv`.
- `bench5_setup.py` — tools, output dirs and local DB paths for this run.
- `Snakefile` — defines the per-benchmark variables and `include`s the shared
  rule library; contains no tool rules of its own.
- `run.sh` — submit the whole benchmark to the aqua queue.

## Shared rules (`../rules/`)

The read generation, biobox/OPAL conversion and per-tool profiling rules are
factored into `../rules/*.smk` and shared by `include:`, rather than copy-pasted
as in the older benchmarks:

- `../rules/data_generation.smk` — `generate_community_and_reads` (ART).
- `../rules/common.smk` — `truth_condensed_to_biobox`, `tool_condensed_to_biobox`, `opal`.
- `../rules/profilers.smk` — singlem, metaphlan, sylph, singlem-regime3.

Each `.smk` documents the variables the including Snakefile must define.

## Running

Generate just the reads + ground truth:

```bash
PYTHONPATH=.. pixi run snakemake --directory 5_known_species_low_coverage \
  -s 5_known_species_low_coverage/Snakefile \
  --configfile 5_known_species_low_coverage/config-8threads.yaml \
  -c 4 truths/known5.condensed.biobox
```

Run the full 4-tool comparison on the queue:

```bash
5_known_species_low_coverage/run.sh
```

Outputs land in `output_<tool>/opal/known5.opal_report` (per-rank OPAL metrics),
with runtime/RAM in `benchmarks/<tool>/known5-8threads.benchmark`.
