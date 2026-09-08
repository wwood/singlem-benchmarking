# Benchmark 9 — one abundant unknown among nine known species

A 10-genome community that is mostly known, with a single novel member that is
also one of the more abundant. Nine genomes are known in GTDB r207 and so are in
every tool's r207 database; the tenth is new in r214 and absent from all of them.

The distinguishing feature is that the unknown is **not** a trace member: it is
carried at ~10× against a community whose known members sit between 1.3× and
27.6×. The question is what a profiler spends that abundance on — does it report
the lineage at the highest rank it can actually support (here genus,
`g__Kineosporia`) and leave the species unassigned, or does it assign it to the
nearest known species? A species-level error on an abundant member costs
proportionally far more in abundance-weighted metrics than the same error on a
faint one, which is what separates this benchmark from the novelty tests where the
novel fraction is either small or evenly spread.

Because the other nine members are ordinary known species at comfortable coverage
— every one of them well above singlem's ~1.5× marker-gene floor — any
species-level loss is attributable to the one unknown rather than to sensitivity.
That is the complement of benchmarks 5/6, where every member is below the floor
and nothing is novel.

**Interpret at both the species and the genus level.** Species level measures the
nine known members plus any species invented for the unknown; genus level shows
whether the unknown was placed correctly at the rank it could be placed at.

## The community

The 10 genomes, their coverages and their truth taxonomy are all in `community.tsv`.
The truth it produces (`truths/marine0.condensed`) has 10 lineages — nine resolved to
species, one stopping at genus (the novel member):

| coverage | lineage |
|---|---|
| 27.6× | `s__Veillonella tobetsuensis` |
| 15.8× | `s__SZUA-592 sp011772965` |
| **10.2×** | **`g__Kineosporia` (genus only — the novel member)** |
| 2.42× | `s__Streptococcus mitis_BM` |
| 2.42× | `s__Sphingomonas koreensis` |
| 2.11× | `s__UBA1174 sp900556855` |
| 1.83× | `s__Listeria_A newyorkensis` |
| 1.34× | `s__Streptococcus lactarius` |
| 1.33× | `s__Flavobacterium hydatis` |
| 1.33× | `s__Beduini massiliensis` |

Reads are simulated into `reads/marine0.{1,2}.fq.gz` (2 × 150 bp, ~54/64 MB gzipped
— a small, fast dataset) and are gitignored.

## Layout

Reads and truth are simulated by `../rules/data_generation_community.smk`, not by
`../rules/data_generation.smk`: the latter looks the truth taxonomy up in the GTDB
r207 metadata, and the novel genome has no r207 metadata row, so its inner join would
drop exactly the member under test. `community.tsv` therefore states coverage *and*
truth taxonomy directly, and there is no `coverage_definitions/` or `genome_list.tsv`.

- `community.tsv` — the whole community: genome, fasta, coverage, role, truth taxonomy.
- `marine0.condensed` — the ground truth recorded from the original run, which
  `community.tsv` reproduces exactly. Kept as a reference; the workflow reads and
  writes `truths/marine0.condensed`.
- `bench9_setup.py` — tools, output dirs and local DB paths.
- `Snakefile` — `include`s the shared rule library.
- `run.sh` — submit the whole benchmark to the aqua queue.

The nine known genomes come from `../reference_genomes/shadow/`; the novel one is
downloaded from NCBI into `../novel_r214_genomes/` by
`../gather_tool_databases.smk`'s `bench89_novel_genomes_download`, so run that first
in a fresh checkout.

Four tools are run — `singlem`, `sylph`, `singlem-regime3`, `metaphlan` — matching
benchmarks 5 and 6. Shared rules in `../rules/` provide database staging,
biobox/OPAL conversion and the per-tool profiling; see `../ARCHITECTURE.md`.

## Running

```bash
9_abundant_unknown/run.sh
```

Or locally, ground truth biobox only:

```bash
PYTHONPATH=.. pixi run snakemake --directory 9_abundant_unknown \
  -s 9_abundant_unknown/Snakefile --configfile 9_abundant_unknown/config-8threads.yaml \
  -c 4 truths/marine0.condensed.biobox
```

Outputs land in `output_<tool>/opal/marine0.opal_report`, with runtime/RAM in
`benchmarks/<tool>/marine0-8threads.benchmark`.

## Results

Sample `marine0`, 8 threads. The truth has 9 species and 10 genera; TP/FP/FN are
OPAL's presence/absence counts at that rank.

The numbers below were produced from the original simulated reads. A fresh run draws
new reads from `community.tsv` (ART is not seeded), so expect small shifts in
Bray–Curtis; the design and the truth are unchanged.

**Species level** — the nine known members, plus any species invented for the unknown:

| tool | Bray–Curtis ↓ | F1 ↑ | Purity | Compl. | TP | FP | FN |
|---|---|---|---|---|---|---|---|
| SingleM regime3 | **0.013** | **1.000** | 1.000 | 1.000 | 9 | 0 | 0 |
| sylph | 0.083 | **1.000** | 1.000 | 1.000 | 9 | 0 | 0 |
| SingleM | 0.206 | 0.857 | 0.750 | 1.000 | 9 | 3 | 0 |
| MetaPhlAn | 0.341 | 0.941 | 1.000 | 0.889 | 8 | 0 | 1 |

**Genus level** — whether the abundant unknown was placed at the rank it can be placed at:

| tool | Bray–Curtis ↓ | F1 ↑ | Purity | Compl. | TP | FP | FN |
|---|---|---|---|---|---|---|---|
| SingleM | **0.044** | **1.000** | 1.000 | 1.000 | 9 | 0 | 0 |
| SingleM regime3 | 0.064 | **1.000** | 1.000 | 1.000 | 9 | 0 | 0 |
| sylph | 0.153 | 0.941 | 1.000 | 0.889 | 8 | 0 | 1 |
| MetaPhlAn | 0.392 | 0.875 | 1.000 | 0.778 | 7 | 0 | 2 |

Every tool recovers the nine known species well — they are all comfortably above
the marker-gene floor — so, as designed, the whole spread comes from the one
abundant unknown, and the three tools split three ways on it:

- **SingleM** detects `g__Kineosporia` (5.49×) but also emits
  `s__Kineosporia sp018499875` (1.82×) — it spends part of the unknown's abundance
  on a species that is not there. That plus two spurious *Veillonella* congeners
  gives 3 species-level false positives and purity 0.750, while its genus level is
  perfect (Bray–Curtis 0.044, F1 1.000). The over-commitment is only partial:
  singlem keeps most of the abundance at genus and leaks about a quarter of it into
  the species call.
- **sylph** commits the opposite error: it never reports *Kineosporia* at all, at
  any rank, so it pays a genus-level false negative (F1 0.941) but stays perfectly
  pure at species level. An abundant organism absent from its database is simply
  invisible rather than misassigned.
- **regime3** is the only tool that gets both ranks right — 1.000 F1 at species
  *and* genus — reporting `g__Kineosporia` at 5.55× with no species call beneath
  it, and it has the lowest species-level abundance error of the four (0.013).
- **MetaPhlAn** misses two members outright, including the abundant novel one and
  `s__SZUA-592 sp011772965` (a known species at 15.8×), and badly distorts the
  remainder: it puts 70% of the community on *V. tobetsuensis*, hence the worst
  Bray–Curtis at both ranks.

Runtime, 8 threads: regime3 ~1 min total, singlem 1 min, sylph 18 s, MetaPhlAn
1.8 min. Peak RSS: MetaPhlAn 17.5 GB, sylph 8.3 GB (whole db in RAM), singlem
1.3 GB. This is a small dataset — the whole benchmark runs end to end in minutes,
which makes it a cheap first check when changing the joint method's rank behaviour.
