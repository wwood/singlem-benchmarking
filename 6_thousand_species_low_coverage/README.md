# Benchmark 6 — 1000 known species, none above 1× coverage

A complex synthetic community in which **every** member is below the coverage
floor. Benchmark 5 showed that at 0.2× a five-species community is invisible to
vanilla singlem and recoverable only through the whole-genome tools; this
benchmark asks the same question of a community with a thousand members, where
there is no abundant taxon for a profiler to coast on and a thousand chances to
invent a species that is not there.

It separates the two things benchmark 5 conflates — *being below the marker-gene
floor* and *being one of only a handful of taxa* — and unlike `3_cami2_marine` it
has no high-abundance head, so per-taxon abundance error cannot be dominated by a
few well-covered genomes.

## The community

1000 **distinct** GTDB r207 species, one genome each:

| | |
|---|---|
| species | 1000 (all distinct; one genome per species) |
| coverage | 0.044× – **1.000×**, median 0.221× |
| total coverage | 253.7× |
| relative abundance | max **0.394%** (cap 10%) |
| median genome size | 3.44 Mbp |
| simulated sequence | ~0.94 Gbp (2 × 150 bp, HSXt error model) |

Both requested caps hold, but only one of them binds: with a thousand members the
1× coverage cap is the active constraint, and the 10% relative-abundance cap is
satisfied with a ~25-fold margin. They cannot both bind at this richness — a
species at 10% relative abundance and ≤1× coverage would imply a total community
coverage of 10×, leaving the other 999 species averaging 0.009×, far below any
tool's detection limit. The distribution is instead lognormal (σ=0.5), rescaled so
the most abundant genome sits exactly at 1×; that gives a ~25-fold dynamic range
with a realistic skew while keeping every member under singlem's ~1.5× marker-gene
sensitivity.

Genomes come from `../reference_genomes/shadow/`, the pool of non-representative
GTDB r207 strains the other benchmarks draw on. Each selected genome's *species*
is in the r207 tool databases (represented there by a different strain), so as in
benchmark 5 every species is findable in principle: what is under test is whether
it is found at these coverages, not whether it is in the database.

## Layout

- `generate_definitions.py` — selects the community. Seeded; run once, outputs committed.
- `genome_list.tsv` — the 1000 genomes and their FASTA paths.
- `coverage_definitions/known1000.tsv` — per-genome coverage. Sample `known1000`.
- `bench6_setup.py` — tools, output dirs and local DB paths.
- `Snakefile` — per-benchmark variables; `include`s the shared rule library.
- `run.sh` — submit the whole benchmark to the aqua queue.
- `plot.ipynb` — Bray–Curtis and F1 comparisons (R kernel), built up two → three →
  four tools; writes `accuracy_results.csv` and `species_level_summary.csv`.
- `make_plot_notebook.py` — regenerates `plot.ipynb` from readable R source.

Shared rules in `../rules/` (`data_generation.smk`, `staging.smk`, `common.smk`,
and one `<tool>_run.smk` per profiler) provide read generation, database staging,
biobox/OPAL conversion and the per-tool profiling; see `../ARCHITECTURE.md` for
the division of labour.

To regenerate the community definition (only needed if the selection or the
distribution should change — the committed files are what the benchmark uses):

```bash
cd 6_thousand_species_low_coverage
pixi run -e art python3 generate_definitions.py
```

Note that `generate_community.py` pairs coverages with genomes positionally
against a shuffled metadata merge, so which genome receives which coverage is
randomised at read-generation time and is not the pairing in
`coverage_definitions/known1000.tsv`. Only the multiset of coverages carries
over. Both caps are properties of that multiset — one species per genome — so
they hold regardless; the ground truth in `truths/known1000.condensed` records
the pairing actually simulated.

## Running

Reads and ground truth only:

```bash
PYTHONPATH=.. pixi run snakemake --directory 6_thousand_species_low_coverage \
  -s 6_thousand_species_low_coverage/Snakefile \
  --configfile 6_thousand_species_low_coverage/config-8threads.yaml \
  -c 16 truths/known1000.condensed
```

The full four-tool comparison on the queue:

```bash
6_thousand_species_low_coverage/run.sh
```

Outputs land in `output_<tool>/opal/known1000.opal_report`, with runtime/RAM in
`benchmarks/<tool>/known1000-8threads.benchmark`.

## Results

`plot.ipynb` plots these; re-execute it after a rerun with

```bash
pixi run jupyter nbconvert --to notebook --execute --inplace \
  6_thousand_species_low_coverage/plot.ipynb
```

The notebook plots Bray–Curtis and **recall** rather than F1: purity is 0.94–0.99
for every tool here, so F1 compresses the one axis on which they differ. At the
species level the truth holds exactly 1000 species, so recall is *species found /
1000*. Purity and F1 are still shown in the final four-metric panel and in the
tables.

Species level, sample `known1000` (1000 species present):

| tool | Bray–Curtis ↓ | Recall ↑ | species found | Purity | F1 | FP |
|---|---|---|---|---|---|---|
| SingleM regime3 | 0.084 | 0.954 | 954 / 1000 | 0.986 | 0.970 | 14 |
| sylph | 0.090 | 0.952 | 952 / 1000 | 0.982 | 0.967 | 17 |
| MetaPhlAn | 0.301 | 0.739 | 739 / 1000 | 0.940 | 0.828 | 47 |
| SingleM | 0.838 | 0.049 | 49 / 1000 | 0.942 | 0.093 | 3 |

Purity is high for every tool, so the whole spread is in recall: below the
marker-gene floor vanilla SingleM sees 49 of 1000 species while what little it
reports is almost all correct. The joint method closes that gap almost entirely
(954) and now edges past sylph (952) on both recall and abundance error, with the
fewest false positives of the four. MetaPhlAn sits between vanilla SingleM and the
whole-genome methods and carries the most false positives — the "thousand chances
to invent a species" this benchmark was designed to expose costs it 47.

These regime3 numbers are from the current
`singlem_sylph_condense_regime` submodule (`c9c5675`, joint pinning + novel
budget) run with weebill `-u` and `condense --alpha 1`. Earlier runs of this
benchmark, before those changes, gave regime3 recall 0.870 and Bray–Curtis 0.229.

## A note on the singlem-regime3 configuration

`../rules/singlem_regime3_run.smk` now runs weebill **with** `-u`
(`--estimate-unknown`), so condense receives a `True_cov` profile already on
singlem's scale and is given `--alpha 1` explicitly, together with
`--joint-pin-sylph-species` and `--joint-novel-budget`.

This resolves the calibration problem earlier versions of this benchmark
documented: without `-u`, condense received `Eff_cov` and had to calibrate alpha
itself, which is biased low at exactly these coverages (on benchmark 5 it landed
near 0.09 where the community's own totals imply ~0.75) — leaving relative
abundances intact but inflating absolute genome-equivalent coverage several-fold.
Since this benchmark lives entirely in the regime where that difference shows up,
it is the one most affected by the fix: regime3's species-level Bray–Curtis
improved from 0.229 to 0.084 and recall from 0.870 to 0.954.
