# Benchmark 10 — novelty within the genus, at two abundances

Four genomes in two unrelated genera. Each genus contributes **one species known
in GTDB r207** — so present in every tool's database — at **50×**, plus **one
congener whose species is new in r214** and therefore absent from all of them. The
novel member is carried at **100×** in one genus and **10×** in the other, so the
same question is asked once with the unknown dominating its known congener 2:1 and
once with it at a fifth of it.

The question is what `singlem-regime3` does with novelty **once the sylph species
layer has been accounted for**. Regime3 gets a whole-genome species profile from
weebill and marker-gene evidence from `singlem pipe`, and `condense --joint` has to
decide what is left over after the sylph species are subtracted. Here the leftover
has a known answer: it is exactly the two novel congeners, and it should surface as
genus-level abundance under `g__Streptococcus` and `g__Collinsella` — not as extra
abundance heaped onto the known congener (the 50× member is the nearest database
neighbour and the easiest thing to mistake it for), and not as nothing at all.

Holding genus, novelty depth and the known congener fixed while varying **only the
unknown's abundance** is what isolates the abundance dependence of that leftover.
Benchmark 8 crosses novelty with congenericity but reads out at one abundance per
member; benchmark 9 has one abundant unknown with no known congener at all.

**Interpret at both the species and the genus level.** Species level is the two
known members plus anything invented for the two unknowns; genus level is whether
the unknowns were placed at the rank they can be placed at, and with how much
abundance. Because both known members are at 50× — far above singlem's ~1.5×
marker-gene floor — any species-level loss is attributable to the novel members
rather than to sensitivity.

## The community

`community.tsv`, four genomes, 210× total:

| coverage | role | genome | truth lineage |
|---|---|---|---|
| 50× | known | `GCF_000964315.1` | `s__Streptococcus equinus` |
| **100×** | **novel** | `GCA_019794795.1` | **`g__Streptococcus` (genus only)** |
| 50× | known | `GCF_902375985.1` | `s__Collinsella ihuae` |
| **10×** | **novel** | `GCA_019795105.1` | **`g__Collinsella` (genus only)** |

The two novel genomes come from `../2_phylogenetic_novelty/genomes/`, the pool of
post-r207 genomes benchmark 2 draws on; both are confirmed absent from the r207
metadata and classify (by gtdb-tk against r207,
`../2_phylogenetic_novelty/stratified_gtdbtk_r207.tsv`) to a known genus with an
unassigned species — exactly the genus-known / species-novel case this benchmark
needs. The two known genomes come from `../reference_genomes/shadow/`, the
non-representative r207 strain pool the other benchmarks use, so each known species
is in the databases represented by a *different* strain.

**The known members must be unambiguously assignable**, or the benchmark measures
GTDB's species boundaries rather than the tools. Both were screened against
`../shadow_vs_gtdb.skani.csv` for ANI to their own species' r207 representative
minus ANI to the nearest *other* species' representative — `s__Streptococcus
equinus` **+7.88** (97.10 vs 89.22), `s__Collinsella ihuae` **+15.89** (100.00 vs
84.11). The first draft used `GCA_900762955.1` (`s__Collinsella sp900758375`), which
is **closer to another species' representative than to its own** (94.98 vs 95.19,
margin −0.21): all four tools duly called it a different *Collinsella*, so that leg
of the design measured nothing. If the community is ever changed, keep the margin
comfortably positive.

The resulting truth (`truths/congeneric4.condensed`) has **four lineages, two of
them genus-only**. At genus level that is 71.4% *Streptococcus* / 28.6%
*Collinsella*; at species level only the two known members appear, holding 23.8%
each — the other 52.4% of the community legitimately has no species-level answer.

**One novel species per genus, deliberately.** Two novel congeners in the same
genus would share the identical genus-only truth lineage and be summed into one
row, so the 10× and 100× levels would stop being separable in the truth. One per
genus keeps four genomes ↔ four scorable lineages.

## Layout

Reads are simulated here, but **not** by `../rules/data_generation.smk`. That rule
takes coverage from `coverage_definitions/<sample>.tsv` and taxonomy from the GTDB
metadata; two of these four genomes have no r207 metadata row, so its inner join
against the metadata would silently drop exactly the members under test —
simulating a two-genome community and writing a truth that omits the novel half.
This benchmark therefore has:

- `community.tsv` — the whole community: genome, fasta, coverage, role and **truth
  taxonomy**. Coverages are used as written, not shuffled against a randomised
  metadata merge as the shared generator does, because here which genome gets which
  coverage *is* the experiment.
- `generate_community.py` — same ART invocation as
  `../1_novel_strains/generate_community.py` (HSXt, 150 bp paired, `-m 400 -s 10`),
  reading coverage and taxonomy from `community.tsv`.
- `bench10_setup.py` — tools, output dirs and local DB paths.
- `Snakefile` — the local `generate_community_and_reads` rule; `include`s the shared
  rule library for everything downstream.
- `run.sh` — submit the whole benchmark to the aqua queue.

Four tools are run — `singlem`, `sylph`, `singlem-regime3`, `metaphlan` — matching
benchmarks 5/6/8/9. Everything from the condensed profile onward (DB staging,
biobox conversion, OPAL) is the shared pipeline in `../rules/`; see
`../ARCHITECTURE.md`.

## Running

```bash
10_congeneric_novelty/run.sh
```

Reads and ground truth only (fits a 4-CPU allocation; ~1 min):

```bash
PYTHONPATH=.. pixi run snakemake --directory 10_congeneric_novelty \
  -s 10_congeneric_novelty/Snakefile \
  --configfile 10_congeneric_novelty/config-8threads.yaml \
  -c 4 truths/congeneric4.condensed
```

Reads are ~126/150 MB gzipped and are gitignored (regenerable, and the generator is
seed-free only in ART's error sampling — the community itself is fully specified by
`community.tsv`).

Outputs land in `output_<tool>/opal/congeneric4.opal_report`, with runtime/RAM in
`benchmarks/<tool>/congeneric4-8threads.benchmark`.

## Results

Sample `congeneric4`, 8 threads. The truth has 2 species and 2 genera; TP/FP/FN are
OPAL's presence/absence counts at that rank.

| tool | rank | Bray–Curtis ↓ | F1 ↑ | Purity | Compl. | TP | FP | FN |
|---|---|---|---|---|---|---|---|---|
| SingleM | species | **0.031** | 1.000 | 1.000 | 1.000 | 2 | 0 | 0 |
| SingleM regime3 | species | 0.033 | 1.000 | 1.000 | 1.000 | 2 | 0 | 0 |
| sylph | species | 0.355 | 1.000 | 1.000 | 1.000 | 2 | 0 | 0 |
| MetaPhlAn | species | 0.677 | 0.500 | 0.500 | 0.500 | 1 | 1 | 1 |
| SingleM | genus | **0.009** | 1.000 | 1.000 | 1.000 | 2 | 0 | 0 |
| SingleM regime3 | genus | 0.012 | 1.000 | 1.000 | 1.000 | 2 | 0 | 0 |
| sylph | genus | 0.214 | 1.000 | 1.000 | 1.000 | 2 | 0 | 0 |
| MetaPhlAn | genus | 0.317 | 1.000 | 1.000 | 1.000 | 2 | 0 | 0 |

Presence/absence is nearly uninformative here — with two genera and two species,
three of the four tools score F1 1.000 at both ranks. The design's discrimination is
almost entirely in **where the abundance went**, so the genus split is the table to
read:

| tool | `g__Streptococcus` (truth 71.4%) | `g__Collinsella` (truth 28.6%) |
|---|---|---|
| SingleM | 71.44% | 26.86% |
| SingleM regime3 | 69.50% | 28.98% |
| **sylph** | **50.00%** | **50.00%** |
| MetaPhlAn | 39.72% | 60.28% |

And the species level, where 52.4% of the community correctly has *no* answer:

| tool | `s__S. equinus` (truth 23.8%) | `s__C. ihuae` (truth 23.8%) | other |
|---|---|---|---|
| SingleM | 24.96% | 25.66% | — |
| SingleM regime3 | 25.44% | 25.44% | — |
| sylph | 50.00% | 50.00% | — |
| MetaPhlAn | — | 60.28% | 39.72% `s__S. lutetiensis` (FP) |

**Reading the three questions the benchmark was built to ask:**

1. **Is the novel abundance placed at genus?** For singlem and regime3, yes, and
   accurately: regime3 reports `g__Streptococcus` 69.50% / `g__Collinsella` 28.98%
   against a truth of 71.4/28.6, and puts only ~25% on each known species — i.e. it
   holds back the ~52% that has no species answer instead of spending it. Both keep
   perfect species purity: **no tool invented a species for either novel congener**,
   which is the nearest-neighbour error the congeneric design makes easy to commit.

2. **The 100× vs 10× contrast — the actual question about the post-sylph leftover.**
   regime3's genus figures land near the truth in *both* genera, so the leftover is
   **not strongly abundance-dependent** across a 10× span: it is neither too tight at
   100× (which would have shown as `g__Streptococcus` falling short) nor too loose at
   10× (which would have inflated `g__Collinsella`). The residual is a mild ~2-point
   shift of *Streptococcus* → *Collinsella*, i.e. very slightly *over*-crediting the
   faint novel member relative to the abundant one.

3. **What the whole-genome-only baseline does.** sylph's failure is the clean
   illustration of why the joint method exists: it reports exactly the two known
   species at 50/50 and **nothing else at any rank**, so the 52% of the community
   that is novel is simply invisible to it. That single omission is the entire
   difference between its Bray–Curtis (0.355 species / 0.214 genus) and regime3's
   (0.033 / 0.012) — and note that the genus-level error is *not* a detection failure
   (F1 1.000, both genera found) but purely a proportion error, since sylph
   renormalises the missing 52% across the two genomes it does see.

MetaPhlAn is the one tool that commits the nearest-neighbour error: it calls the
known *S. equinus* `s__S. lutetiensis` (a species-level FP and FN), and inverts the
genus proportions (39.7/60.3 against 71.4/28.6).

Notably **vanilla SingleM is at its best in exactly this regime** — marginally the
best Bray–Curtis at both ranks — which is the mirror image of benchmarks 5/6: every
member here is at 10–100×, far above its ~1.5× marker-gene floor, and marker genes
place novel organisms at the highest supportable rank by construction. This
benchmark is therefore *not* where the joint method wins on sensitivity; it is where
it has to avoid *losing* the correctness singlem already has while adding sylph's
species resolution, and it does (0.033 vs 0.031 species, 0.012 vs 0.009 genus).

Runtime, 8 threads: regime3 ~1.4 min total (weebill 3.9 s + pipe 74 s + condense
8 s), singlem 89 s, sylph 19 s, MetaPhlAn 2.9 min. Peak RSS: MetaPhlAn 17.5 GB,
sylph 8.2 GB, singlem 2.5 GB, regime3 1.2 GB (weebill, the largest of its three
steps). A small, fast benchmark — useful as a quick check when changing the joint
method's rank behaviour.
