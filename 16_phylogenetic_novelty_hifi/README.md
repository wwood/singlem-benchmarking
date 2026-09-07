# Benchmark 16 — novel lineage + known species, PacBio HiFi

Benchmark 2's phylogenetic-novelty community — one **novel** lineage plus one known
species at equal abundance — sequenced as simulated PacBio HiFi reads. This is
benchmark 15 with the error rate removed and everything else held fixed: the same 120
samples read from benchmark 2's own committed tables, the same six novelty depths, the
same coverages and fragment lengths; only Badread's error and qscore models change.

## The community

Identical to benchmark 15 (see its README for how the ground truth is defined):

| | |
|---|---|
| samples | 120, one per novel genome |
| members per sample | 1 novel genome + 1 known species partner |
| coverage | 10× each (equal abundance, 50% each) |
| novelty depths | species / genus / family / order / class / phylum — **20 samples each** |
| read length | mean ~14 kb (Badread `--length 15000,13000`, as benchmark 15) |
| read identity | ~99.9% (`--identity 30,3` — mean qscore 30, stdev 3) |
| error / qscore model | `pacbio2021` (Badread's Sequel IIe HiFi models) |

## Why this benchmark exists: the 15-vs-16 gradient

Placing a novel lineage means measuring how far a sequence diverges from its nearest
database relative. At the shallow end of the novelty gradient that divergence — roughly
5% in ANI terms for a novel species — is the **same order as Nanopore's per-base error
rate**, so error and novelty are confounded in benchmark 15 in a way they are not at HiFi
accuracy.

The expected signature is therefore a 15-vs-16 gap that is:

- **largest at species/genus novelty**, where true divergence and Nanopore error are
  comparable and a tool cannot easily tell them apart;
- **narrowest at class/phylum novelty**, where the true divergence dwarfs any error rate
  and both technologies should reach the same conclusion.

Reading the two benchmarks across all six novelty depths is the intended analysis; a
single depth in isolation says much less.

## Interpretation

**Read at a high rank — kingdom/domain, not species**, as for benchmarks 2 and 15. By
construction one member has no correct species label available, so the question is
whether the profiler places it at the highest rank it can support rather than inventing a
known species for it or dropping it entirely.

Because the reads here are close to error-free, this benchmark is also the fairest
available test of whether a *long* read helps novelty placement in principle: a ~14 kb
read spans whole marker genes and much of their genomic context, which is strictly more
evidence about divergence than a 150 bp fragment carries.

## Structural differences from benchmark 2

The same two as benchmark 15, both following from the tool set: there is no `{method}`
wildcard (all three tools report relative abundance, so one truth per sample, and the
shared `../rules/common.smk` scoring rules apply), and reads are single-end. Read
generation is local to this benchmark rather than the shared long-read rule, because the
novel member of every pair is absent from r207 and a metadata inner join would drop
exactly the genome under test. Benchmark 15's README explains both in full.

## Tools

`singlem` (at 0.21.3, via the `singlem-longread` environment), `sylph` and
`singlem-regime3` — the long-read-capable subset of benchmark 2's eleven tools. See
benchmark 11's README for the metaphlan exclusion and the SingleM version choice; the
2-vs-16 comparison carries the same SingleM version caveat as 2-vs-15.

## Cost

120 samples × (read simulation + 3 tools) = 1804 jobs. Read simulation is 2 genomes at
10× per sample, so a few minutes each. `run.sh` passes `--keep-going` so one failed
sample does not stop the rest.

## Layout

- `bench16_setup.py` — tools, output dirs, local DB paths, `longread_read_tech`,
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
./16_phylogenetic_novelty_hifi/run.sh             # aqua queue, one PBS job per rule
```

A single sample can be run locally to smoke-test the pipeline:

```bash
PYTHONPATH=.. pixi run snakemake \
  --directory 16_phylogenetic_novelty_hifi \
  -s 16_phylogenetic_novelty_hifi/Snakefile \
  --configfile 16_phylogenetic_novelty_hifi/config-8threads.yaml -c 8 \
  output_singlem/opal/GCA_013154095.1_genomic.opal_report
```
