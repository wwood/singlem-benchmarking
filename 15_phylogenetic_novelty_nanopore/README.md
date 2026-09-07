# Benchmark 15 — novel lineage + known species, Oxford Nanopore R10.4.1

Benchmark 2's phylogenetic-novelty community — one **novel** lineage plus one known
species at equal abundance — sequenced as simulated Nanopore R10.4.1 reads instead of
150 bp paired Illumina.

The community definition is benchmark 2's, read **directly from its committed tables**
(`stratified_genome_metadata.tsv`, `stratified_gtdbtk_r207.tsv`, `genome_pairs.tsv`) and
its genome directories rather than copied into this benchmark, so the two cannot drift
apart.

## The community

| | |
|---|---|
| samples | 120, one per novel genome |
| members per sample | 1 novel genome + 1 known species partner |
| coverage | 10× each (equal abundance, 50% each) |
| novelty depths | species / genus / family / order / class / phylum — **20 samples each** |
| read length | mean ~14 kb (Badread `--length 15000,13000`) |
| read identity | ~95% (`--identity 95,99,2.5`) |
| error / qscore model | `nanopore2023` (Badread's R10.4.1 models) |

Each sample's novel genome is new to GTDB at the stated rank; its ground-truth taxonomy
is its **r207** classification — the lineage truncated at whatever rank still exists in
r207 — because that is the most specific answer any r207-backed tool could legitimately
give. (The r214 taxonomy in benchmark 2's metadata table is the full lineage the genome
eventually received, which no tool here could report.)

The equal-abundance design is what makes this a test of *placement* rather than
quantification: whatever a tool does with the novel member, it cannot be excused by that
member being rare.

## Interpretation

**Read at a high rank — kingdom/domain, not species**, exactly as for benchmark 2. By
construction one member has no correct species label available (often no genus or family
either), so species-level accuracy is not the question. What is tested is whether a
profiler **detects that something is present and places it at the highest rank it can
support**, instead of mis-assigning it to a known species or dropping it.

### The long-read angle

Novelty detection is where read length and read accuracy pull hardest against each other.
Placing a novel lineage means resolving how far a sequence diverges from anything in the
database — and Nanopore's ~5% per-base error rate is the *same order* as the divergence
being measured at the shallow end of the gradient (a novel species is roughly 5%
divergent in ANI terms). Errors can therefore masquerade as novelty and novelty as
errors.

Benchmark 16 repeats this at HiFi accuracy, where that confound largely disappears. The
expected signature is a 15-vs-16 gap that is **largest at species/genus novelty and
narrows toward phylum**, where the true divergence dwarfs any error rate. Reading the two
benchmarks across the six novelty depths is the intended analysis.

## Structural differences from benchmark 2

Both follow from the tool set, and both simplify this benchmark relative to benchmark 2:

1. **No `{method}` wildcard.** Benchmark 2 runs eleven tools split into
   relative-abundance and reads-wise (read-count) families and scores them against two
   different truth files. All three tools here report relative abundance, so there is
   exactly one truth per sample — which is what lets this benchmark use the shared
   `../rules/common.smk` scoring rules instead of benchmark 2's bespoke `{method}` copies.
   If a read-count tool (kraken / kaiju / metabuli) is ever added, the reads-wise truth
   half must be restored; see `../bin/generate_longread_2component_community.py`.
2. **Reads are single-end**, so the `../rules/*_longread_run.smk` rules apply.

Read generation is **local to this benchmark**, not `../rules/longread_data_generation.smk`
— that shared rule takes taxonomy from the GTDB metadata by inner join, and the novel
member of every pair is by construction absent from r207, so the join would silently drop
exactly the genome under test. Coverage and truth taxonomy are therefore stated per
genome, as in benchmark 2. (This is the same reason benchmark 10 has its own generator.)

## Tools

`singlem` (at 0.21.3, via the `singlem-longread` environment), `sylph` and
`singlem-regime3` — the long-read-capable subset of benchmark 2's eleven tools.
metaphlan, motus, kraken, kaiju, sourmash, map2b, metabuli and metaphlan42 are all
absent. Benchmark 11's README explains the metaphlan and SingleM-version choices in full;
note that the 2-vs-15 comparison therefore also carries a SingleM version change
(0.18.0 → 0.21.3), because 0.18.0 predates long-read support entirely.

## Cost

120 samples × (read simulation + 3 tools) = 1804 jobs. Read simulation is ~2 genomes at
10× per sample, so a few minutes each — much cheaper per sample than benchmarks 13/14,
but there are 120 of them. `run.sh` passes `--keep-going` so one failed sample does not
stop the rest.

## Layout

- `bench15_setup.py` — tools, output dirs, local DB paths, `longread_read_tech`,
  `bench2_dir`.
- `Snakefile` — reads benchmark 2's community tables, defines the local read-generation
  rule, `include`s the shared rule library.
- `run.sh` — submit the whole benchmark to the aqua queue.

The community tables and genome FASTAs live in `../2_phylogenetic_novelty/`; if that
benchmark's genomes have not been downloaded yet, see the root `README.md`.

Shared rules in `../rules/` (`staging.smk`, `common.smk`, and one
`<tool>_longread_run.smk` per profiler) provide database staging, biobox/OPAL conversion
and the per-tool profiling; see `../ARCHITECTURE.md`.

## Running

```bash
./15_phylogenetic_novelty_nanopore/run.sh         # aqua queue, one PBS job per rule
```

A single sample can be run locally to smoke-test the pipeline:

```bash
PYTHONPATH=.. pixi run snakemake \
  --directory 15_phylogenetic_novelty_nanopore \
  -s 15_phylogenetic_novelty_nanopore/Snakefile \
  --configfile 15_phylogenetic_novelty_nanopore/config-8threads.yaml -c 8 \
  output_singlem/opal/GCA_013154095.1_genomic.opal_report
```
