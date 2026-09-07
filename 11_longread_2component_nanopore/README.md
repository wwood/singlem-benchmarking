# Benchmark 11 — two known species at 10×, Oxford Nanopore R10.4.1

The first long-read benchmark in this repository, and deliberately the easiest one.
Two known GTDB r207 species at equal, comfortable coverage, sequenced as simulated
Nanopore R10.4.1 reads. Nothing about the community is hard: both members have a
correct species answer in every tool's database, both sit well above every tool's
sensitivity floor, there is no novelty, and neither is a close relative of the other.

That is the point. This benchmark does not ask *how accurate* a profiler is — it asks
whether long-read input works at all, and what it costs relative to Illumina. If a tool
loses a species here, the cause is the read technology, because the community offers no
other explanation. It is the control for benchmarks 13 and 15, which carry the same
reads into a genuinely hard community (1000 sub-1× species) and a novel-lineage
community respectively.

Benchmark 12 is this same community and coverage as PacBio HiFi, so the pair isolates
read accuracy (~95% vs ~99.9%) at matched read length and depth.

## The community

| | |
|---|---|
| members | 2 known GTDB r207 species |
| coverage | 10× each (equal abundance, 50% each) |
| genomes | `GCF_001953075.1` (*Halioglobus lutimaris*), `GCF_003815695.1` (*Chryseobacterium carnipullorum*) |
| genome sizes | 5.00 Mbp, 5.46 Mbp |
| simulated sequence | ~105 Mbp, ~7 000 reads |
| read length | mean ~14.8 kb (Badread `--length 15000,13000`) |
| read identity | ~95% (`--identity 95,99,2.5`, beta distribution) |
| error / qscore model | `nanopore2023` (Badread's R10.4.1 models) |

Both genomes come from `../reference_genomes/shadow/`, the pool of non-representative
GTDB r207 strains the other benchmarks draw on, so each species is in every tool's
database represented by a *different* strain — findable in principle, as in benchmarks 5
and 6.

Both were **ANI-screened** against `../shadow_vs_gtdb.skani.csv` before selection: each
one's top-ANI r207 representative is its own species (99.1% and 98.5%), with the next
species more than 6 ANI points away. Benchmark 10's README records why this matters — a
"known" genome closer to another species' representative than to its own makes every
tool wrong for reasons that have nothing to do with the tools.

## Read simulation

Reads come from [Badread](https://github.com/rrwick/Badread) 0.4.2 via
`../bin/generate_longread_community.py`, whose defaults the Badread documentation
describes as corresponding to "Oxford Nanopore R10.4.1 reads of mediocre quality". The
parameters are passed **explicitly** rather than left to Badread's defaults, so a future
Badread release changing its defaults cannot silently change what this benchmark
simulates.

Long reads are **single-end**, which is why this benchmark uses the
`../rules/*_longread_run.smk` rules rather than the paired-end ones: one
`reads/<sample>.fq.gz` instead of `reads/<sample>.{1,2}.fq.gz`. Everything from the
condensed profile onward is the standard shared pipeline.

## Tools

`singlem`, `sylph` and `singlem-regime3` — the long-read-capable subset of the tools the
other benchmarks run.

**metaphlan is absent.** The pinned 4.0.6 has no long-read mode at all (`--long_reads`
arrived in MetaPhlAn 4.2), so running it here would measure the missing feature rather
than the tool. Adding it would mean the `metaphlan42` environment plus its ~20 GB
database, which is not currently downloaded.

**SingleM runs at 0.21.3, not 0.18.0.** Long-read support was added in SingleM 0.20.0.
The `singlem` pixi environment the short-read benchmarks use is pinned to 0.18.0, which
runs on long reads without erroring but recovers substantially less: on a single-genome
10× Nanopore test it reported 0.79× at species level where 0.21.3 reported 1.93×.
Benchmarks 11–16 therefore use a separate `singlem-longread` environment at 0.21.3;
benchmarks 1–10 are untouched and keep 0.18.0. This does mean the SingleM *version*
differs between the long-read and short-read benchmarks, so a 6-vs-13 or 2-vs-15
comparison carries that caveat — the alternative was to handicap SingleM with a version
that predates the feature under test.

`singlem-regime3` needs no such swap: it is built from the
`singlem_sylph_condense_regime` submodule, already at 0.21.3.

## Interpretation

**Read at the species level.** Both members have a correct species answer available, so
anything short of recovering both at roughly equal abundance is a failure.

With only two members presence/absence saturates easily, so the informative readouts are
(a) abundance accuracy and (b) whether the profile is polluted with spurious extra
species. The latter is the long-read-specific risk for a marker-gene method: a ~5%
per-base error rate gives SingleM many more chances to mis-assign a 60 bp window than
Illumina does.

## Layout

- `genome_list.tsv` — the two genomes and their FASTA paths.
- `coverage_definitions/nanopore2.tsv` — per-genome coverage. Sample `nanopore2`.
- `bench11_setup.py` — tools, output dirs, local DB paths, `longread_read_tech`.
- `Snakefile` — per-benchmark variables; `include`s the shared rule library.
- `run.sh` — submit the whole benchmark to the aqua queue.

Shared rules in `../rules/` (`longread_data_generation.smk`, `staging.smk`,
`common.smk`, and one `<tool>_longread_run.smk` per profiler) provide read generation,
database staging, biobox/OPAL conversion and the per-tool profiling; see
`../ARCHITECTURE.md` for the division of labour.

## Running

```bash
./11_longread_2component_nanopore/run.sh          # aqua queue, one PBS job per rule
```

or locally, which is quite feasible for this benchmark (Badread takes ~1 minute per
genome at 10×):

```bash
PYTHONPATH=.. pixi run snakemake \
  --directory 11_longread_2component_nanopore \
  -s 11_longread_2component_nanopore/Snakefile \
  --configfile 11_longread_2component_nanopore/config-8threads.yaml -c 8
```
