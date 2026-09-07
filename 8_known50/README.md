# Benchmark 8 — half novel, half Streptomyces

A 12-genome community with two crossed axes:

- **six known / six novel.** Six genomes are known in GTDB r207 and so are in
  every tool's r207 database; the other six are new in r214 and therefore absent
  from all of them.
- **half congeneric.** Independently of the above, six of the 12 come from a
  single genus, *Streptomyces*.

The crossing is the point. The novel half asks whether a profiler withholds a
species call and places the organism at the highest rank it can support, rather
than snapping it onto the nearest known species; the *Streptomyces* block makes
that failure mode easy to commit, because a novel *Streptomyces* has several
closely related known congeners available to be mistaken for. Where benchmark 2
tests novelty against a backdrop of one unrelated known species, this one tests it
against near neighbours.

Coverages span four orders of magnitude (0.1× – 200×, plus a 1013× genus-level
*Streptomyces* lineage in the truth), so the question is asked at both ends of the
abundance range rather than at a single operating point.

**Interpret at both the species and the genus level.** Species-level accuracy
measures the known half and penalises species calls invented for the novel half;
genus level shows whether the novel members were nonetheless placed correctly, in
the manner of benchmark 2's high-rank reading. Reading species level alone would
score a correct "genus known, species withheld" answer the same as a miss.

## The community

The ground truth (`sample5.condensed`) has 10 lineages: six resolved to species,
four stopping at genus (the novel members, whose species has no r207/r214-truth
label to be scored against).

| coverage | lineage |
|---|---|
| 1013× | `g__Streptomyces` (genus only) |
| 200× | `s__Streptomyces rimosus` |
| 100× | `s__Dickeya poaceiphila` |
| 30× | `s__Streptomyces lunaelactis` |
| 20× | `s__VadinCA11 sp002498365` (Archaea) |
| 2× | `g__SIG402` (genus only) |
| 1× | `g__Flavobacterium` (genus only) |
| 0.6× | `s__Streptococcus equinus` |
| 0.5× | `g__Paenibacillus` (genus only) |
| 0.1× | `s__Streptomyces sp900091845` |

Reads: `sample5.{1,2}.fq.gz` (2 × 150 bp, ~2.6/3.1 GB gzipped).

## Layout

Unlike benchmarks 5/6/7 the reads and truth are **not** simulated here — they ship
with the dataset. `stage_provided_reads` / `stage_provided_truth` in the Snakefile
symlink them into the `reads/` and `truths/` layout the shared rules address, so
`../rules/data_generation.smk` is not included and there is no
`coverage_definitions/` or `genome_list.tsv`.

- `sample5.{1,2}.fq.gz`, `sample5.condensed` — the provided reads and ground truth.
- `bench8_setup.py` — tools, output dirs and local DB paths.
- `Snakefile` — staging rules for the provided data; `include`s the shared rule library.
- `run.sh` — submit the whole benchmark to the aqua queue.

Four tools are run — `singlem`, `sylph`, `singlem-regime3`, `metaphlan` — matching
benchmarks 5 and 6. Shared rules in `../rules/` provide database staging,
biobox/OPAL conversion and the per-tool profiling; see `../ARCHITECTURE.md`.

## Running

```bash
8_known50/run.sh
```

Or locally, ground truth biobox only:

```bash
PYTHONPATH=.. pixi run snakemake --directory 8_known50 \
  -s 8_known50/Snakefile --configfile 8_known50/config-8threads.yaml \
  -c 4 truths/sample5.condensed.biobox
```

Outputs land in `output_<tool>/opal/sample5.opal_report`, with runtime/RAM in
`benchmarks/<tool>/sample5-8threads.benchmark`.

## Results

Sample `sample5`, 8 threads. The truth has 6 species and 7 genera; TP/FP/FN are
OPAL's presence/absence counts at that rank.

**Species level** — the known half, plus any species invented for the novel half:

| tool | Bray–Curtis ↓ | F1 ↑ | Purity | Compl. | TP | FP | FN |
|---|---|---|---|---|---|---|---|
| SingleM regime3 | **0.095** | **0.923** | 0.857 | 1.000 | 6 | 1 | 0 |
| SingleM | 0.430 | 0.727 | 0.800 | 0.667 | 4 | 1 | 2 |
| sylph | 0.592 | 0.923 | 0.857 | 1.000 | 6 | 1 | 0 |
| MetaPhlAn | 0.593 | 0.571 | 0.500 | 0.667 | 4 | 4 | 2 |

**Genus level** — whether the novel members were nonetheless placed correctly:

| tool | Bray–Curtis ↓ | F1 ↑ | Purity | Compl. | TP | FP | FN |
|---|---|---|---|---|---|---|---|
| SingleM | **0.004** | **1.000** | 1.000 | 1.000 | 7 | 0 | 0 |
| SingleM regime3 | 0.020 | **1.000** | 1.000 | 1.000 | 7 | 0 | 0 |
| sylph | 0.258 | 0.833 | 1.000 | 0.714 | 5 | 0 | 2 |
| MetaPhlAn | 0.297 | 0.833 | 1.000 | 0.714 | 5 | 0 | 2 |

The two ranks tell opposite stories for vanilla SingleM, which is why both are
needed. At genus level it is essentially perfect (Bray–Curtis 0.004, all 7 genera,
no false positives) — the novel members are detected and placed correctly. But at
species level it commits precisely the error the *Streptomyces* block was built to
provoke: the 1013× novel *Streptomyces* is assigned to a known congener,
`s__Streptomyces sp002920635`, at 510× coverage. Reading species level alone would
misdescribe a tool that got the lineages right and only over-committed on rank;
reading genus level alone would hide the over-commitment entirely.

sylph and regime3 both leave the novel *Streptomyces* at genus and recover all 6
known species, so they tie on species-level F1 — but their abundance error differs
sevenfold (0.095 vs 0.592), because sylph systematically under-reports coverage for
the abundant members (142× vs a true 200× for *S. rimosus*, 72× vs 100× for
*Dickeya*). The joint method inherits sylph's correct rank behaviour while keeping
singlem's calibration, which is what the Bray–Curtis gap measures. sylph's weaker
genus-level score is the flip side: it withholds two genus-level lineages that
SingleM's marker genes do detect.

MetaPhlAn is worst on both axes here, and its 4 species-level false positives at
purity 0.500 mean half its species calls are for organisms that are not present.

Runtime, 8 threads: sylph 2.5 min, singlem 8 min, regime3 20 min total (16.5 min
of it `singlem pipe`), MetaPhlAn 46 min. MetaPhlAn's cost is bowtie2 over the
5.8 GB concatenated read file; peak RSS 17.5 GB against sylph's 8.6 GB (whole db in
RAM) and singlem's 3.4 GB.
