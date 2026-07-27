# Benchmark 7 — the GTDB r232 metapackage as the backing database

Every other benchmark in this repository is backed by GTDB **r207** databases.
This one asks the same question as benchmark 6 — recovery of a complex community
in which every member is below the coverage floor — but with each reference
database swapped to GTDB **r232**. It isolates the effect of the database release
(and, for singlem, the newer S6.5.0 metapackage) from the community design, which
is held fixed.

Only the three tools that have an r232 database are run:

| tool | r232 database |
|---|---|
| `singlem` | `S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb` (metapackage format v6) |
| `singlem-regime3` | the same r232 metapackage + weebill `r232.100.syl2db` (two-stage sylph fork db) |
| `sylph` | `r232.100.syldb` (standalone sylph db) |

`metaphlan` and `coverm` are omitted: they have no r232 database here.

## The community

Reused **verbatim** from benchmark 6: `genome_list.tsv` and
`coverage_definitions/known1000.tsv` are copies of benchmark 6's. See
`../6_thousand_species_low_coverage/README.md` for the community design (1000
distinct GTDB r207 species, lognormal coverages 0.04×–1.0×, median ~0.22×, none
above singlem's ~1.5× marker-gene floor).

### r207 → r232 attrition

The community was selected against GTDB r207. Between r207 and r232, **119 of the
1000 genomes were removed from GTDB** (true removals or superseded accessions) and
so are absent from the r232 metadata. `generate_community.py` builds the community
by an inner join against the metadata it is given, so with r232 metadata the
community actually simulated is the **881 surviving genomes**. Because the reads
*and* the ground truth are produced by that same script in one pass, they are
generated from the same 881-genome multiset and remain consistent — the dropped
genomes are simply never simulated. This is deliberate: a genome that no longer
exists in r232 has no r232 lineage to be scored against, so keeping it would make
the truth incoherent rather than more complete.

The coverage multiset is therefore the 881-genome subset of benchmark 6's, not
identical to it; the two benchmarks are close but not a controlled genome-for-genome
pair. What is controlled is the community *design* (same selection procedure, same
coverage distribution) and the genomes themselves where they survive.

## Ground truth taxonomy

The truth is written in **r232** taxonomy (from `../bac120_metadata_r232.tsv` /
`../ar53_metadata_r232.tsv`), and sylph's genome hits are mapped to lineages via
the r232 taxonomy (`../{bac120,ar53}_taxonomy_r232.tsv`). Tool output and truth
therefore share a taxonomy, so a lineage that was renamed or moved between r207
and r232 (e.g. `p__Proteobacteria` → `p__Pseudomonadota`) is not counted as an
error. **Interpret at the species level**, as in benchmark 6.

## Layout

- `bench7_setup.py` — tools, output dirs, and the r232 overrides of the source db
  paths, metadata and sylph taxonomy.
- `Snakefile` — per-benchmark variables; `include`s the shared rule library.
- `genome_list.tsv`, `coverage_definitions/known1000.tsv` — community (from bench 6).
- `run.sh` — submit the whole benchmark to the aqua queue.

The shared rules in `../rules/` provide read generation, biobox/OPAL conversion
and the per-tool profiling. This benchmark uses the same `singlem_run`,
`sylph_run` and `singlem_regime3_run` rules as the others; the r232 databases and
taxonomy are injected purely through the variables set in `bench7_setup.py`.

## Reference data

The r232 databases are staged outside git:

- singlem metapackage: copied into `../tool_reference_data/` (staged from
  `/work/microbiome/db/singlem/`, which compute nodes cannot see).
- sylph / weebill dbs: already on weka at
  `/scratch/microbiome/woodcrob/non_sensitive/weebill_dbs/r232.100.{syldb,syl2db}`.
- r232 taxonomy / metadata: staged to the repo root next to the r207 copies
  (`../{bac120,ar53}_{taxonomy,metadata}_r232.tsv`).

## Running

Reads and ground truth only:

```bash
PYTHONPATH=.. pixi run snakemake --directory 7_r232_metapackage \
  -s 7_r232_metapackage/Snakefile \
  --configfile 7_r232_metapackage/config-8threads.yaml \
  -c 16 truths/known1000.condensed
```

The full three-tool comparison on the queue:

```bash
7_r232_metapackage/run.sh
```

Outputs land in `output_<tool>/opal/known1000.opal_report`, with runtime/RAM in
`benchmarks/<tool>/known1000-8threads.benchmark`.
