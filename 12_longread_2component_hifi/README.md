# Benchmark 12 — two known species at 10×, PacBio HiFi

Benchmark 11's community, sequenced as PacBio HiFi instead of Nanopore R10.4.1. Same two
genomes, same 10× each, same simulated fragment-length distribution; only the error and
qscore models change. Read accuracy goes from ~95% to ~99.9% with read length and depth
held fixed, which is what makes the 11-vs-12 pair an accuracy experiment rather than a
technology-branding one.

## The community

Identical to benchmark 11 (see its README for the ANI screening of both members):

| | |
|---|---|
| members | 2 known GTDB r207 species |
| coverage | 10× each (equal abundance, 50% each) |
| genomes | `GCF_001953075.1` (*Halioglobus lutimaris*), `GCF_003815695.1` (*Chryseobacterium carnipullorum*) |
| simulated sequence | ~105 Mbp |
| read length | mean ~14 kb (Badread `--length 15000,13000`, as benchmark 11) |
| read identity | ~99.9% (`--identity 30,3` — mean qscore 30, stdev 3) |
| error / qscore model | `pacbio2021` (Badread's Sequel IIe HiFi models) |

The `--identity 30,3` form matters: two values (rather than benchmark 11's three) switch
Badread from sampling identities from a beta distribution to sampling **qscores** from a
normal distribution, which the Badread docs recommend for high-accuracy data such as
HiFi or ONT duplex. It is what makes these reads HiFi-accurate rather than merely long.

## Interpretation

**Read at the species level, and primarily against benchmark 11.** Alone this benchmark
says little: both members are known species well above every tool's floor, so a tool
should simply recover both. Read as a pair with benchmark 11, it says how much of
whatever benchmark 11 loses is attributable to Nanopore's error rate.

That distinction bears on the two method families differently. SingleM must recover a
60 bp marker-gene window accurately enough to place it, so per-base errors hit it
directly; sylph's k-mer containment degrades more gracefully but still needs k-mers to
survive intact. A gap between 11 and 12 for a marker-gene method that closes for a
k-mer method would be the expected signature.

## Tools

`singlem` (at 0.21.3, via the `singlem-longread` environment), `sylph` and
`singlem-regime3`. metaphlan is absent because the pinned 4.0.6 has no long-read mode.
Benchmark 11's README explains both choices in full.

## Layout

- `genome_list.tsv` — the two genomes and their FASTA paths (same as benchmark 11).
- `coverage_definitions/hifi2.tsv` — per-genome coverage. Sample `hifi2`.
- `bench12_setup.py` — tools, output dirs, local DB paths, `longread_read_tech`.
- `Snakefile` — per-benchmark variables; `include`s the shared rule library.
- `run.sh` — submit the whole benchmark to the aqua queue.

Shared rules in `../rules/` (`longread_data_generation.smk`, `staging.smk`,
`common.smk`, and one `<tool>_longread_run.smk` per profiler) provide read generation,
database staging, biobox/OPAL conversion and the per-tool profiling; see
`../ARCHITECTURE.md`.

## Running

```bash
./12_longread_2component_hifi/run.sh              # aqua queue, one PBS job per rule
```

or locally (Badread takes ~1–2 minutes per genome at 10×):

```bash
PYTHONPATH=.. pixi run snakemake \
  --directory 12_longread_2component_hifi \
  -s 12_longread_2component_hifi/Snakefile \
  --configfile 12_longread_2component_hifi/config-8threads.yaml -c 8
```
