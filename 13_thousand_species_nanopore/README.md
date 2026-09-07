# Benchmark 13 — 1000 species below 1× coverage, Oxford Nanopore R10.4.1

Benchmark 6's community — 1000 distinct known GTDB r207 species, every one below 1×
coverage — sequenced as simulated Nanopore R10.4.1 reads instead of 150 bp paired
Illumina.

`genome_list.tsv` and `coverage_definitions/nanopore1000.tsv` are byte-for-byte copies of
benchmark 6's (verified by checksum), so the community, the per-genome coverages and the
ground truth are held fixed and the read technology is the single variable. Read this
benchmark against benchmark 6 the way benchmark 7 is read against benchmark 6: same
community, one axis changed — there the database release, here the sequencing chemistry.

## The community

Identical to benchmark 6; see its README for how the distribution was chosen.

| | |
|---|---|
| species | 1000 (all distinct; one genome per species) |
| coverage | 0.044× – 1.000×, median 0.221× |
| total coverage | 253.7× |
| relative abundance | max 0.394% |
| median genome size | 3.44 Mbp |
| simulated sequence | ~0.94 Gbp |
| read length | mean ~14 kb (Badread `--length 15000,13000`) |
| read identity | ~95% (`--identity 95,99,2.5`) |
| error / qscore model | `nanopore2023` (Badread's R10.4.1 models) |

### 996 genomes, not 1000

Four of the 1000 shadow-genome FASTAs are **0 bytes** (failed downloads):
`GCA_012275145.1`, `GCF_002846295.1`, `GCF_004337435.2`, `GCF_004402975.1`. Badread
cannot simulate from an empty reference — it divides its progress counter by the target
base count and dies — so `../bin/generate_longread_community.py` drops them, with a
warning in the read-generation log, **from both the reads and the truth**.

Dropping them from the truth as well is the important half. Benchmark 6 (ART) does *not*
fail on an empty reference, so those four species stay in its truth while contributing no
reads — meaning every tool is scored as having missed four species that were never
sequenced. Benchmark 6's truth therefore has 1000 lineages where this benchmark's has
996. That is a small but real difference in the 6-vs-13 comparison, in benchmark 13's
favour; the four affected species are all below 0.4% relative abundance.

## Why the coverage floor is the interesting place to change read length

The same number of bases arrives in ~100× fewer, ~100× longer reads. For the two method
families that cuts in opposite directions:

- **Marker-gene (SingleM).** A 0.22× genome yields very few reads to begin with; making
  each one 100× longer means ~100× fewer independent chances to land on any *particular*
  60 bp marker-gene window. Each read that does land spans whole genes, but coverage of a
  specific short window is what SingleM needs, and that is a Poisson process with far
  fewer draws.
- **K-mer containment (sylph).** Cares about total bases rather than how they are
  packaged, so the repackaging should be close to neutral — but its k-mers must survive a
  ~5% per-base error rate, which at k=31 removes most k-mers from any given read.

Benchmark 14 repeats this at HiFi accuracy, which separates the length effect (shared by
13 and 14) from the error effect (present only in 13).

## Interpretation

**Read at the species level**, as a paired comparison with benchmark 6. As there, no
member reaches SingleM's ~1.5× marker-gene sensitivity and there is no abundant taxon to
coast on, so per-taxon abundance error cannot be dominated by a few well-covered genomes
— and there are a thousand chances to invent a species that is not present.

## Tools

`singlem` (at 0.21.3, via the `singlem-longread` environment), `sylph` and
`singlem-regime3` — the long-read-capable subset of benchmark 6's four tools, and
coincidentally exactly benchmark 7's trio.

**metaphlan is absent**, unlike in benchmark 6: the pinned 4.0.6 has no long-read mode
(`--long_reads` arrived in MetaPhlAn 4.2). Benchmark 11's README explains this and the
SingleM version choice in full. Note that the 6-vs-13 comparison therefore also carries a
SingleM version change (0.18.0 → 0.21.3), because 0.18.0 predates long-read support
entirely.

## Cost

**Badread is the long pole.** It is single-threaded per genome and simulates roughly
0.4 Mbp/s, so ~0.94 Gbp over 1000 genomes takes several hours even at 8-way concurrency
(the `generate_longread_community_and_reads` rule requests 24 h of queue time so a slow
node cannot lose the whole simulation). The resulting reads are also ~10× larger on disk
than benchmark 6's. Genomes are simulated concurrently via `extern.run_many`, and each
gets a distinct derived `--seed`, so a rerun reproduces the same reads.

## Layout

- `genome_list.tsv` — the 1000 genomes and their FASTA paths (copy of benchmark 6's).
- `coverage_definitions/nanopore1000.tsv` — per-genome coverage. Sample `nanopore1000`.
- `bench13_setup.py` — tools, output dirs, local DB paths, `longread_read_tech`.
- `Snakefile` — per-benchmark variables; `include`s the shared rule library.
- `run.sh` — submit the whole benchmark to the aqua queue.

Shared rules in `../rules/` (`longread_data_generation.smk`, `staging.smk`,
`common.smk`, and one `<tool>_longread_run.smk` per profiler) provide read generation,
database staging, biobox/OPAL conversion and the per-tool profiling; see
`../ARCHITECTURE.md`.

## Running

```bash
./13_thousand_species_nanopore/run.sh             # aqua queue, one PBS job per rule
```

Running this one locally is not recommended — the read simulation alone is a multi-hour,
single-node job.
