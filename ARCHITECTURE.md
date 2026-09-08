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
   reads (simulated or staged; paired short reads, or single-end long reads in 11-16)
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
| `rules/data_generation.smk` | generic ART read simulation (`generate_community_and_reads`) from a `<sample>.tsv` coverage file. Truth taxonomy comes from the GTDB metadata by inner join, so it can only simulate genomes that *are* in the target release. |
| `rules/data_generation_community.smk` | the same ART simulation driven by a committed `community.tsv` (genome, fasta, coverage, role, truth taxonomy) instead of a coverage file plus a metadata join, via `bin/generate_community_from_tsv.py`. Used by benchmarks 8, 9 and 10, whose communities all contain genomes newer than the target release and so have no metadata row to join against. |
| `rules/longread_data_generation.smk` | the long-read counterpart: Badread simulation from the same coverage-file inputs, with `longread_read_tech` ∈ {`nanopore`, `pacbio-hifi`}. Emits **one** `<sample>.fq.gz`, not a pair. Same metadata-join caveat as the ART rule. |
| `rules/staging.smk` | database-copy / index rules for every tool, **opt-in** via a `staged_tools` set. |
| `rules/singlem_run.smk` | SingleM `pipe` → condensed. |
| `rules/metaphlan_run.smk` | MetaPhlAn → SGB report → GTDB profile → condensed. |
| `rules/sylph_run.smk` | `sylph profile` → condensed. |
| `rules/kraken_run.smk` | kraken2 → bracken (per rank) → condensed. |
| `rules/sourmash_run.smk` | sourmash sketch/gather/tax → condensed. |
| `rules/metabuli_run.smk` | metabuli classify → condensed. |
| `rules/coverm_run.smk` | strobealign mapping vs all GTDB genomes → per-genome coverage → condensed. |
| `rules/singlem_regime3_run.smk` | the joint method: weebill (sylph fork) `--two-stage` + SingleM `pipe --no-sylph` + `condense --joint`. |
| `rules/singlem_longread_run.smk` | SingleM on single-end long reads, in the `singlem-longread` env (0.21.3, not 0.18.0). |
| `rules/sylph_longread_run.smk` | `sylph profile -r` (single-end) → condensed. |
| `rules/singlem_regime3_longread_run.smk` | the joint method with single-end reads (`pipe -1`, `weebill -r`). |

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
and notebooks, plus one isolated environment per benchmarked tool. Read simulation has
its own environments: `art` (art_illumina, benchmarks 1–10) and `badread` (Badread,
benchmarks 11–16). `singlem-longread` pins SingleM 0.21.3 for the long-read benchmarks,
alongside the `singlem` environment's 0.18.0 used everywhere else. Reference
databases are gathered by `gather_tool_databases.smk` and located via
`tool_reference_data.py`. Three git submodules provide SingleM (`singlem`), the
joint-method SingleM branch (`singlem_sylph_condense_regime`) and the sylph fork
(`weebill`).


## The benchmarks and how to interpret each

All of them share the pipeline above; they differ in the **community design** and
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

### 7 · `7_r232_metapackage/` — the same community on GTDB R232 (species level)

Benchmark 6's community with **every reference database swapped from R207 to
R232**: the S6.5.0 R232 SingleM metapackage, the `r232.100.syl2db` weebill db and
the `r232.100.syldb` sylph db. Only the three tools that have an R232 database are
run (`singlem`, `singlem-regime3`, `sylph`); the truth is written in **R232**
taxonomy so legitimate R207→R232 reclassifications are not scored as errors. 119
of the 1000 genomes no longer exist in R232, so the community actually simulated
is the 881 survivors — reads and truth are generated in one pass and stay
consistent.

**Interpret at the species level**, as a database-release comparison against
benchmark 6 rather than as a fresh community design. All the R232 wiring is
injected through variables in `bench7_setup.py` (source db paths, metadata,
`sylph_bac_taxonomy`, `sylph_mem_mb`, `sylph_c`); the run rules are unchanged.

### 8 · `8_known50/` — half novel, half congeneric (species **and genus** level)

12 genomes: **six known in R207, six new in R214** (absent from every R207
database), crossed with a second axis — six of the 12 are from a single genus,
*Streptomyces*. Coverages span 0.1×–200× plus a 1013× genus-level *Streptomyces*
lineage. Single sample `sample5`.

**Interpret at species *and* genus level.** The novel half asks whether a tool
withholds the species call and places the organism at the highest rank it can
support; the congeneric block makes the opposing failure mode easy to commit,
since a novel *Streptomyces* has close known relatives to be mistaken for. Where
benchmark 2 tests novelty against one unrelated known species, this tests it
against near neighbours, and at both ends of a four-order-of-magnitude abundance
range. Species level alone would score "genus right, species withheld" — the
correct answer — the same as a miss.

### 9 · `9_abundant_unknown/` — one abundant unknown (species **and genus** level)

10 genomes: **nine known in R207 plus one new in R214**, and the novel one is
*abundant* (~10× against known members at 1.3×–27.6×), appearing in the truth as a
genus-only lineage (`g__Kineosporia`). Single sample `marine0`; small, fast reads
(~54/64 MB).

**Interpret at species *and* genus level.** The question is what a profiler spends
that abundance on: the highest supportable rank, or the nearest known species — an
error that costs proportionally far more on an abundant member than on a faint
one. Because the other nine are ordinary known species comfortably above the
marker-gene floor, any species-level loss is attributable to the one unknown
rather than to sensitivity. It is the complement of benchmarks 5/6, where
everything is below the floor and nothing is novel.

### 10 · `10_congeneric_novelty/` — novelty within the genus, at two abundances (species **and genus** level)

Four genomes in two unrelated genera. Each genus contributes one species **known in
r207** at **50×** plus one congener whose species is **new in r214** (absent from
every r207 database); the novel member is carried at **100×** in *Streptococcus* and
**10×** in *Collinsella*. Truth (`truths/congeneric4.condensed`) has four lineages,
two of them genus-only. Single sample `congeneric4`.

**Interpret at species *and* genus level.** The question is what
`singlem-regime3` does with novelty *once the sylph species layer is accounted
for*: after `condense --joint` subtracts the sylph species, the leftover has a
known answer — exactly the two novel congeners — and it should surface as
genus-level abundance rather than being heaped onto the known congener (the nearest
database neighbour) or dropped. Holding genus, novelty depth and the known congener
fixed while varying only the unknown's abundance isolates the abundance dependence
of that leftover, which benchmarks 8 (one abundance per member) and 9 (no known
congener) do not separate.

Like benchmarks 8 and 9 it simulates its reads via
`rules/data_generation_community.smk` rather than `rules/data_generation.smk`: the
latter takes taxonomy from the GTDB metadata by inner join, and two of the four
genomes have no r207 metadata row, so it would silently drop exactly the members
under test. Instead `community.tsv` states coverage *and* truth taxonomy per genome,
and the coverages are used as written rather than shuffled positionally — here which
genome gets which coverage *is* the experiment.

With only two genera and two species, presence/absence is nearly saturated (three of
four tools score F1 1.000 at both ranks), so this benchmark is read on **abundance
placement**, not detection. In the first run it also showed how a badly chosen
"known" member can void a whole leg of a design: a genome closer to another species'
r207 representative than to its own made every tool wrong at species level for
reasons that had nothing to do with the tools. Known members are now ANI-screened
against `shadow_vs_gtdb.skani.csv`; see the benchmark README.

### 11–16 · the long-read family

Benchmarks 11–16 re-ask three existing questions with **long reads** instead of 150 bp
paired Illumina, in matched Nanopore/HiFi pairs. Reads are simulated with
[Badread](https://github.com/rrwick/Badread) 0.4.2; `nanopore` uses its default
`nanopore2023` models (which its docs describe as R10.4.1 of mediocre quality,
`--identity 95,99,2.5`), and `pacbio-hifi` uses `pacbio2021` with `--identity 30,3` — two
values, which switches Badread to sampling qscores from a normal distribution, the
recommended form for HiFi-accuracy data. Parameters are passed explicitly by
`bin/generate_longread_community.py` rather than left to Badread's defaults, so a future
Badread release cannot silently change what is simulated.

| # | Illumina original | read tech | question inherited |
|---|---|---|---|
| 11 | *(new: control)* | Nanopore R10.4.1 | does long-read input work at all, and what does it cost |
| 12 | *(new: control)* | PacBio HiFi | the same, error-free |
| 13 | 6 (1000 species <1×) | Nanopore R10.4.1 | the coverage floor at scale |
| 14 | 6 (1000 species <1×) | PacBio HiFi | the same, error-free |
| 15 | 2 (novel lineage + known sp.) | Nanopore R10.4.1 | detect & place a new lineage high in the tree |
| 16 | 2 (novel lineage + known sp.) | PacBio HiFi | the same, error-free |

Three design decisions run through all six:

1. **Each pair varies only accuracy.** Within a pair the genomes, coverages and simulated
   fragment lengths are identical, so a 13-vs-14 or 15-vs-16 difference is attributable to
   the ~95% vs ~99.9% per-base accuracy. Across a pair *and* its Illumina original
   (6/13/14 or 2/15/16), the read-*length* effect — a fixed number of bases delivered in
   ~100× fewer, ~100× longer reads — is shared by both long-read arms and so separable
   from the error effect. This matters most for marker-gene methods: SingleM needs reads
   covering a specific 60 bp window, and long reads give far fewer independent chances to
   hit one.
2. **Reads are single-end.** This is why the long-read benchmarks need their own
   `rules/*_longread_run.smk` files rather than a flag: every rule's input signature
   changes from `<sample>.{1,2}.fq.gz` to one `<sample>.fq.gz`. Everything from the
   condensed profile onward is the unchanged shared pipeline.
3. **The tool set is `singlem`, `sylph`, `singlem-regime3`** — the long-read-capable
   subset. Two consequences worth stating plainly, since both weaken the comparison to the
   Illumina originals:
   - **metaphlan is absent** (it is in benchmarks 2, 5 and 6). The pinned 4.0.6 has no
     long-read mode at all; `--long_reads` arrived in MetaPhlAn 4.2, whose ~20 GB database
     is not currently downloaded. Adding it means the `metaphlan42` env plus that
     download plus a new run rule.
   - **SingleM runs at 0.21.3, not 0.18.0.** Long-read support landed in SingleM 0.20.0.
     0.18.0 does not error on long reads but recovers substantially less (0.79× vs 1.93×
     species coverage on a single-genome 10× Nanopore test), so benchmarking it would
     measure the version rather than the read technology. Benchmarks 11–16 therefore use a
     separate `singlem-longread` pixi environment; benchmarks 1–10 keep 0.18.0 untouched.
     A 6-vs-13 or 2-vs-15 comparison consequently carries a SingleM version change as
     well as a read-technology change. `singlem-regime3` needs no swap — its submodule is
     already 0.21.3.

Benchmarks 13/14 reuse benchmark 6's `genome_list.tsv` and coverage file byte-for-byte
(verified by checksum). Benchmarks 15/16 go further and read benchmark 2's committed
community tables and genome directories *in place* rather than copying them, so they
cannot drift. Like benchmark 10, 15/16 must define their own read-generation rule: the
shared generator takes taxonomy from the GTDB metadata by inner join, and the novel member
of every pair is by construction absent from r207, so the join would drop exactly the
genome under test. They also drop benchmark 2's `{method}` wildcard, because all three
tools report relative abundance — which is what lets them use the shared `common.smk`
scoring rules.

**Cost note.** Badread is single-threaded per genome at roughly 0.4 Mbp/s. Benchmarks
11/12 and 15/16 simulate two genomes per sample and finish in minutes per sample, but
13/14 simulate ~0.94 Gbp over 1000 genomes and take several hours even at 8-way
concurrency (their read-generation rule requests 24 h of queue time).

### Benchmarks defined by a `community.tsv`

Benchmarks 8, 9 and 10 differ from 5/6/7 in where the community definition lives.
Instead of `coverage_definitions/<sample>.tsv` plus `genome_list.tsv` plus a join
against the GTDB r207 metadata for the truth taxonomy, each commits a single
`community.tsv` giving genome, fasta path, coverage, role and truth taxonomy per
member, and includes `rules/data_generation_community.smk`. They must: every one of
these communities contains genomes whose species is new in r214, which have no r207
metadata row, so `rules/data_generation.smk`'s inner join would drop exactly the
members under test.

Their known members come from the Zenodo shadow pool (`reference_genomes/shadow/`).
Benchmarks 8 and 9's novel members are not in that tarball and were never packaged,
so `gather_tool_databases.smk`'s `bench89_novel_genomes_download` fetches them from
NCBI by accession into `novel_r214_genomes/`; benchmark 10's come from benchmark 2's
`genomes/`. Everything from the condensed profile onward is the standard shared
pipeline.

Benchmarks 8 and 9 were originally run from reads handed over as files, with no
recorded community definition. Their `community.tsv` was reconstructed from those
reads: ART names every read `<source contig>-<index>`, so each read can be traced to
the genome it came from — the shadow pool by contig name, the novel members by NCBI
lookup of the WGS accession (RefSeq copies prefix contigs `NZ_`, GenBank ones do not,
which decides which of the two was used). All 34.5 M and 0.71 M read pairs are
accounted for, every stated coverage is within 3% of the coverage measured from the
reads, and each file reproduces its recorded `<sample>.condensed` exactly. Because
ART is not seeded, a re-run draws new reads and shifts the recorded numbers slightly.

(Benchmark 3 also uses non-simulated reads, but stages real CAMI reads with its own
`stage_reads.sh` and builds its truth from the CAMI metadata.)

### Interpretation at a glance

| # | Benchmark | Community | Read primary rank | Core question |
|---|---|---|---|---|
| 1 | novel_strains | known species, non-rep strains | **species** | recovery of known species (+ runtime/RAM) |
| 2 | phylogenetic_novelty | novel lineage + known sp. | **kingdom / domain** | detect & correctly place a new lineage high in the tree |
| 3 | cami2_marine | CAMI2 marine (real) | **species** | realistic community; sensitivity floor |
| 4 | complex_and_novel | 0–100 % novel gradient | **species** (vs % known) | graceful degradation as novelty rises |
| 5 | known_species_low_coverage | 5 species @ 0.2× | **species** | the coverage floor, controlled |
| 6 | thousand_species_low_coverage | 1000 species, all <1× | **species** | the floor at scale, no abundant head |
| 7 | r232_metapackage | benchmark 6's community, R232 dbs | **species** | effect of the database release, community held fixed |
| 8 | known50 | 12 genomes: 6 novel, 6 *Streptomyces* | **species + genus** | novelty judged against close known relatives |
| 9 | abundant_unknown | 9 known + 1 *abundant* novel | **species + genus** | where an abundant unknown's abundance is spent |
| 10 | congeneric_novelty | 2 genera × (known 50× + novel congener at 100× / 10×) | **species + genus** | the post-sylph leftover: is novelty placed at genus, and does that depend on its abundance |
| 11 | longread_2component_nanopore | 2 known species @ 10×, Nanopore R10.4.1 | **species** | does long-read input work at all (control for 13/15) |
| 12 | longread_2component_hifi | the same, PacBio HiFi | **species** | benchmark 11 without the error rate |
| 13 | thousand_species_nanopore | benchmark 6's community, Nanopore R10.4.1 | **species** | the coverage floor when bases arrive as few long reads |
| 14 | thousand_species_hifi | benchmark 6's community, PacBio HiFi | **species** | separates the length effect from the error effect (vs 13) |
| 15 | phylogenetic_novelty_nanopore | benchmark 2's pairs, Nanopore R10.4.1 | **kingdom / domain** | novelty placement when error ≈ divergence |
| 16 | phylogenetic_novelty_hifi | benchmark 2's pairs, PacBio HiFi | **kingdom / domain** | novelty placement without that confound (vs 15) |


## Case studies (`case_studies/`)

A benchmark answers "how well does this tool do?"; it does not hand you the individual
failing reads to work on. `case_studies/` holds small, committed, regenerable extracts
that do — derived *from* a benchmark, but usable in seconds without rerunning one.

| directory | derived from | what it isolates |
|---|---|---|
| `nanopore_marker_windows/` | benchmark 15, sample `GCA_019347805.1_genomic` (family-level novelty) | the 1059 Nanopore reads that **cover** a 60 bp marker window, labelled with the window's error content and with SingleM's own per-read verdict. 51.6% are unused; 90.8% of those carry an indel. |
| `single_indel_reads/` | two reads from the above | the minimal case: one read per community member, each with **one indel and no substitutions** in the window, each window ±200 bp of context. Both lost; repairing the single base recovers both. |

The two are the same problem at two scales, and are meant to be used in that order: fix
`single_indel_reads/` in the inner loop (two reads, sub-second, no OPAL), then confirm the
change holds over `nanopore_marker_windows/`' 1059 reads and 546 failures — which is also
where the `singlem_query_based` rows guard against regressions.

Each has a `build.sh` that regenerates every committed file byte-for-byte, and a `README.md`
stating the headline numbers and how to score a candidate algorithm against them. Large
intermediates (alignments) go to scratch and are not committed.

Three general tools support this, and work on any benchmark's reads:

| script | purpose |
|---|---|
| `case_studies/find_marker_windows.py` | given `singlem pipe -f <genome>` output, locate each ground-truth 60 bp window in the genome as reference coordinates (handles alignment gaps and duplicated markers). |
| `case_studies/extract_marker_window_reads.py` | given those coordinates, an `--eqx` alignment of the reads, and a SingleM `--archive-otu-table`, emit every read whose alignment *spans* a window — including reads an indel interrupts, which an exact match cannot see — with per-read error class, frameshift status and `singlem_query_based` / `diamond` / `not_found` verdict. |
| `case_studies/extract_window_with_flanks.py` | narrow that to hand-picked `read:gene` pairs and excerpt each to the window plus `--flank` bases of context, emitting the read, the perfect reference span, the truth window and a printed alignment. The flanks are the point: an indel inside a window has to be *detected* from the fact that the surrounding sequence still aligns while the frame does not. |

**Controls matter more than the failing case.** A read that fails proves nothing on its own
— it could be too short, or the marker could be absent from the database. Both case studies
therefore ship a perfect-sequence control (same span, no errors) and, in
`single_indel_reads/`, a window-repaired counterfactual that fixes only the window's indel
and leaves every flank error in place. Those two recover where the real read does not, which
is what makes the single base attributable.

The `--archive-otu-table` is the key input: it carries `read_names` and
`taxonomy_assignment_method`, so a read can be labelled as recovered via the exact window,
recovered only via the translated (DIAMOND) search, or missed outright. Those are three
different asks of a candidate algorithm — rescue, upgrade, and don't-regress respectively.


## Viewing results

Each benchmark has a `plot.ipynb` that reads its `output_<tool>/opal/…` reports
and plots the metrics **at the rank appropriate to that benchmark** (per the
table above). `plot_overall.ipynb` in the repo root aggregates across benchmarks.
The raw OPAL reports contain every rank, so the rank shown in a plot is an
interpretation choice, not a limitation of the data.
