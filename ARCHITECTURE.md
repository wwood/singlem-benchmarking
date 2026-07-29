# Architecture of the metagenome profiler benchmarks

This repository benchmarks microbial community profilers — tools that take
metagenome reads and report a taxonomic profile (which taxa, at what relative
abundance) — against a GTDB **R207**-based ground truth. The focus is comparing
[SingleM](https://github.com/wwood/singlem) (and its `singlem-regime3` joint
variant) to whole-genome and marker-gene competitors.

This document describes how the benchmarks are wired together, the shared rule
library under `rules/`, and — importantly — **how to interpret each benchmark**,
because they are not all meant to be read at the same taxonomic rank.

For *running* instructions see the top-level `README.md` and each benchmark's own
`README.md`. This file is about structure and interpretation.


## The common pipeline

Every benchmark, whatever its community design, funnels each (tool, sample) pair
through the same four-stage pipeline and scores it with the same metric:

```
   reads (simulated or staged)
        │   per-tool profiling  (rules/<tool>_run.smk)
        ▼
   output_<tool>/<tool>/<sample>.profile      ← SingleM "condensed" format
        │   bin/condensed_profile_to_biobox.py (rules/common.smk)
        ▼
   output_<tool>/biobox/<sample>.biobox        ← CAMI Bioboxes profiling format
        │   opal.py  vs  the truth biobox      (rules/common.smk)
        ▼
   output_<tool>/opal/<sample>.opal_report     ← per-rank OPAL metrics
```

The unifying interchange format is the SingleM **condensed profile**: a
per-sample, per-lineage relative-abundance table. Every tool has a small adapter
in `bin/<tool>_to_condensed.py` that converts its native output into this format,
so all downstream scoring is tool-agnostic. The ground truth is produced in the
same condensed format, converted to a biobox once, and reused as both the OPAL
gold standard and the template that aligns each tool's biobox to it.

Scoring is done by [OPAL](https://github.com/CAMI-challenge/OPAL), which emits
metrics **per taxonomic rank** (domain/kingdom → phylum → … → species). This is
why the interpretation rank matters: the same `opal_report` contains every rank,
and each benchmark is designed to be *read* at a particular one (see below).

### Ground truth and abundance conventions

- Truth and most tools report **relative abundance** (sequence-abundance /
  coverage weighted). SingleM's condensed profiles are coverage-based; the
  biobox conversion normalises to relative abundance.
- Some tools report **filled** profiles (every ancestor rank present for each
  leaf); `tools_with_filled_output_profiles` (`kraken`, `sourmash` in
  `tool_reference_data.py`) marks them so the biobox conversion does not re-fill.
- Benchmark 2 additionally distinguishes **relative-abundance** tools from
  **reads-wise** (read-count) tools and scores them against different truth files
  (`{method}` = `relabund` vs `reads_wise`).


## The shared rule library (`rules/`)

The per-tool profiling, database staging and scoring rules were originally
copy-pasted into every benchmark Snakefile. They now live once in `rules/` and
are pulled in with `include:`. Each benchmark's Snakefile is reduced to: define a
handful of variables (tools, output dirs, local DB paths, dataset list, read
locations), then `include` the shared rules it needs.

| file | contents |
|---|---|
| `rules/common.smk` | tool-agnostic scoring: `truth_condensed_to_biobox`, `tool_condensed_to_biobox`, `opal`. |
| `rules/data_generation.smk` | generic ART read simulation (`generate_community_and_reads`) from a `<sample>.tsv` coverage file. |
| `rules/staging.smk` | database-copy / index rules for every tool, **opt-in** via a `staged_tools` set. |
| `rules/singlem_run.smk` | SingleM `pipe` → condensed. |
| `rules/metaphlan_run.smk` | MetaPhlAn → SGB report → GTDB profile → condensed. |
| `rules/sylph_run.smk` | `sylph profile` → condensed. |
| `rules/kraken_run.smk` | kraken2 → bracken (per rank) → condensed. |
| `rules/sourmash_run.smk` | sourmash sketch/gather/tax → condensed. |
| `rules/metabuli_run.smk` | metabuli classify → condensed. |
| `rules/coverm_run.smk` | strobealign mapping vs all GTDB genomes → per-genome coverage → condensed. |
| `rules/singlem_regime3_run.smk` | the joint method: weebill (sylph fork) `--two-stage` + SingleM `pipe --no-sylph` + `condense --joint`. |

### How a benchmark opts in

Each Snakefile defines, before the includes:

- `output_dirs_dict` — `{tool: "output_<tool>"}` for every tool it references.
- The local DB path variables for the tools it stages (e.g. `sylph_db_local`).
- `staged_tools` — the set of tools whose DB-copy rule should be emitted.
  `rules/staging.smk` guards each rule with `if "<tool>" in staged_tools`, so a
  benchmark only ever stages (and only needs path variables for) the databases it
  actually uses. This is what lets the two structurally-divergent benchmarks (2
  and 4) share the identical `cp` rules while keeping their own run rules.
- Then `include:` the specific `rules/*_run.smk` files for the tools it runs.

### Two conventions the shared rules rely on

1. **`shell.prefix("cd {workflow.basedir} && ...")`.** Cluster jobs submitted via
   the aqua profile (`snakemake_mqsub` → PBS) start in `$HOME`, not the workflow
   directory. The prefix `cd`s into the benchmark dir so workdir-relative paths
   (`output_<tool>/…`, `../bac120_taxonomy_r207.tsv`, `bin/…`) and pixi manifest
   discovery (`../pixi.toml`) resolve on compute nodes. Benchmarks 1, 3, 5, 6 set
   it. (Benchmarks 2 and 4 do not currently set a prefix; their shared rules are
   only the `cp` staging rules, whose sources are absolute.)
2. **Per-tool pixi environments.** Rules activate the relevant environment at
   runtime with `eval "$(pixi shell-hook -e <env>)"` (or `pixi run --environment
   <env>`), so no `--use-conda` is needed. `singlem-regime3` rules additionally
   `unset PYTHONPATH` to avoid the repo-root `singlem/` submodule dir shadowing
   the branch's editable install.

### Software / environments

Managed by pixi (`pixi.toml`): a default environment drives Snakemake, plotting
and notebooks, plus one isolated environment per benchmarked tool. Reference
databases are gathered by `gather_tool_databases.smk` and located via
`tool_reference_data.py`. Three git submodules provide SingleM (`singlem`), the
joint-method SingleM branch (`singlem_sylph_condense_regime`) and the sylph fork
(`weebill`).


## The benchmarks and how to interpret each

All six share the pipeline above; they differ in the **community design** and
therefore in the **rank and the question** each is meant to answer.

### 1 · `1_novel_strains/` — known species, novel strains (species level)

Communities simulated (ART, 150 bp paired) from **non-representative** GTDB R207
genomes whose *species* is nonetheless in every tool's database. 10 samples
(`marine0`–`marine9`) with coverage profiles from `coverage_definitions/`.

**Interpret at the species level.** The question is how well tools recover known
species from strains that are not the database representative. This is also the
**runtime / RAM benchmark**: `run_benchmark.sh` runs everything single-threaded,
and `benchmarks/<tool>/…benchmark` records wall time and memory, so the per-tool
thread counts here are deliberate and should not be "optimised". It additionally
computes an `opal_low_abund` variant (profiles with sub-0.02% lineages removed)
to separate detection of abundant vs trace taxa.

### 2 · `2_phylogenetic_novelty/` — new lineages (KINGDOM level)

Two-component communities: one **novel** lineage (a genome new to GTDB) plus one
known species, at equal abundance. Read generation is bespoke
(`generate_2component_community.py`), and the truth is written twice — a
relative-abundance truth and a reads-wise truth — so relative-abundance tools and
read-count tools are each scored against the matching gold standard (`{method}`
wildcard).

**Interpret at a high rank — kingdom/domain, not species.** By construction one
member has *no* correct species (or often genus/family) label available, so
species-level accuracy is not the point. What is tested is whether a profiler can
**detect that something is present and place it at the right high rank** rather
than mis-assigning it to a known species or dropping it. mOTUs and kaiju are run
to kingdom level only here, and metaphlan uses the genome-pairs mapping — all
reflecting that this benchmark lives at the top of the taxonomy.

### 3 · `3_cami2_marine/` — CAMI2 marine (species level)

The CAMI2 marine challenge datasets, with the CAMI taxonomy remapped to GTDB
R207. Reads are **not** simulated here — they are staged from the CAMI download
(`stage_reads.sh`) — and the ground truth is built from the CAMI genome-to-id
map, per-genome lengths and read-mapping counts
(`generate_truth_condensed_format`, bench-3-specific). 9 samples
(`marine1`–`marine9`; `marine0` is an upstream 0-byte placeholder).

**Interpret at the species level** (with the rank sweep in `plot.ipynb` showing
the fall-off toward strain). This is the realistic, externally-defined community:
a broad abundance range with a high-abundance head, useful for seeing where each
tool's sensitivity floor lies. `coverm` (strobealign mapping vs all GTDB genomes)
is included as a mapping-based reference point.

### 4 · `4_complex_and_novel/` — novelty gradient (species level, swept over % known)

A complex community (marine coverages) where a controlled **0–100 %** of the
community is *new* in GTDB R214 relative to R207 (`known{known_percent}` ∈
{0,10,…,100}). Every output path carries a `known{pct}/` layer, which is why this
benchmark's run rules cannot share the flat-path shared ones. 5 samples × 11
novelty fractions.

**Interpret at the species level, as a function of the known fraction.** The
axis of interest is the novelty gradient: as more of the community becomes
unknown, how gracefully does each tool degrade — does it keep the known fraction
right and correctly withhold the unknown fraction, or does it hallucinate known
species for novel organisms? Read across `known{pct}` rather than at a single
point.

### 5 · `5_known_species_low_coverage/` — the coverage floor (species level)

Five **known** GTDB R207 species, each simulated at exactly **0.2× coverage**
(≈20 % relative abundance each). 0.2× is well below SingleM's ~1.5× marker-gene
sensitivity, so this isolates whether whole-genome methods (sylph, and the
`singlem-regime3` joint method) recover a community where vanilla SingleM should
struggle. Single sample `known5`.

**Interpret at the species level.** It is a small, controlled probe of the
**coverage floor**: detection/quantification when every member is faint but
findable-in-principle. This benchmark was the first to use the shared `rules/`
library.

### 6 · `6_thousand_species_low_coverage/` — the floor at scale (species level)

1000 **distinct** known GTDB R207 species, every one below **1× coverage**
(median ~0.22×) and below 0.4 % relative abundance — a lognormal skew rescaled so
the most abundant genome sits at exactly 1×. Community definition is seeded and
committed (`generate_definitions.py`, `coverage_definitions/known1000.tsv`).
Single sample `known1000`.

**Interpret at the species level.** It carries benchmark 5's floor probe up to a
complex community with **no abundant taxon to coast on and a thousand chances to
invent a species that is not there**. Unlike benchmark 3 there is no
high-abundance head, so per-taxon abundance error cannot be dominated by a few
well-covered genomes. It separates the two things benchmark 5 conflates: *being
below the marker-gene floor* and *being one of only a handful of taxa*. See its
README for a note on the `singlem-regime3` alpha calibration in this regime.

### Interpretation at a glance

| # | Benchmark | Community | Read primary rank | Core question |
|---|---|---|---|---|
| 1 | novel_strains | known species, non-rep strains | **species** | recovery of known species (+ runtime/RAM) |
| 2 | phylogenetic_novelty | novel lineage + known sp. | **kingdom / domain** | detect & correctly place a new lineage high in the tree |
| 3 | cami2_marine | CAMI2 marine (real) | **species** | realistic community; sensitivity floor |
| 4 | complex_and_novel | 0–100 % novel gradient | **species** (vs % known) | graceful degradation as novelty rises |
| 5 | known_species_low_coverage | 5 species @ 0.2× | **species** | the coverage floor, controlled |
| 6 | thousand_species_low_coverage | 1000 species, all <1× | **species** | the floor at scale, no abundant head |


## Viewing results

Each benchmark has a `plot.ipynb` that reads its `output_<tool>/opal/…` reports
and plots the metrics **at the rank appropriate to that benchmark** (per the
table above). `plot_overall.ipynb` in the repo root aggregates across benchmarks.
The raw OPAL reports contain every rank, so the rank shown in a plot is an
interpretation choice, not a limitation of the data.
