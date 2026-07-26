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

Shared rules in `../rules/` (`data_generation.smk`, `common.smk`, `profilers.smk`)
provide read generation, biobox/OPAL conversion and the per-tool profiling; see
benchmark 5's README for the division of labour.

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

## A note on the singlem-regime3 configuration

`../rules/profilers.smk` runs weebill without `-u` (`--estimate-unknown`), so
condense receives an `Eff_cov` profile and must calibrate alpha itself. That
calibration is biased low at exactly these coverages — on benchmark 5 it lands
near 0.09 where the community's own totals imply ~0.75, which leaves relative
abundances intact but inflates absolute genome-equivalent coverage several-fold.
Running weebill with `-u` reports `True_cov` instead, already on singlem's scale,
and condense then takes alpha as 1 automatically. Since this benchmark lives
entirely in the regime where the difference shows up, it is worth running both
ways; switching the shared rule affects benchmarks 1, 2, 3 and 5 as well, so it
has been left alone here.
