# Benchmark 14 — 1000 species below 1× coverage, PacBio HiFi

Benchmark 6's community — 1000 distinct known GTDB r207 species, every one below 1×
coverage — sequenced as simulated PacBio HiFi reads. This is benchmark 13 with the error
rate removed and everything else held fixed: the same 1000 genomes, the same per-genome
coverages, the same simulated fragment lengths, only Badread's error and qscore models
change.

## The community

Identical to benchmarks 6 and 13 (`genome_list.tsv` and the coverage file are
byte-for-byte copies of benchmark 6's, verified by checksum).

| | |
|---|---|
| species | 1000 (all distinct; one genome per species) |
| coverage | 0.044× – 1.000×, median 0.221× |
| total coverage | 253.7× |
| simulated sequence | ~0.94 Gbp |
| read length | mean ~14 kb (Badread `--length 15000,13000`, as benchmark 13) |
| read identity | ~99.9% (`--identity 30,3` — mean qscore 30, stdev 3) |
| error / qscore model | `pacbio2021` (Badread's Sequel IIe HiFi models) |

### 996 genomes, not 1000

Four of the 1000 shadow-genome FASTAs are 0 bytes (failed downloads) and are dropped from
both the reads and the truth, exactly as in benchmark 13 — see its README for why this
makes benchmark 6's truth (1000 lineages) differ from these two (996).

## The three-way comparison is the point

|  | read length | accuracy |
|---|---|---|
| benchmark 6 | 150 bp paired | ~99.9% (Illumina) |
| benchmark 13 | ~14 kb | ~95% (Nanopore R10.4.1) |
| **benchmark 14** | ~14 kb | ~99.9% (PacBio HiFi) |

Long reads change two things at once, and 6/13/14 separates them. The *packaging* of a
fixed number of bases into ~100× fewer, ~100× longer reads is shared by 13 and 14, so it
survives the 13-vs-14 comparison; the per-base error rate is present only in 13.

- If a marker-gene method loses ground in 13 and **recovers** it in 14, the cause was
  error.
- If it loses the same ground in **both**, the cause was that a sub-1× genome simply
  produces too few long reads to hit any given 60 bp window — a length effect that no
  improvement in accuracy can fix.

## Interpretation

**Read at the species level**, as the accuracy-controlled arm of the 6/13/14 set. Every
member is below SingleM's ~1.5× marker-gene floor and there is no abundant head, so this
is the coverage floor probed with the most favourable long reads currently available —
close to the best case for long-read profiling at low coverage. A tool that still cannot
recover this community with HiFi reads is limited by read length, not read quality.

## Tools

`singlem` (at 0.21.3, via the `singlem-longread` environment), `sylph` and
`singlem-regime3`. metaphlan is absent because the pinned 4.0.6 has no long-read mode;
see benchmark 11's README for that and for the SingleM version choice.

## Cost

As for benchmark 13: Badread is single-threaded per genome (~0.4 Mbp/s), so ~0.94 Gbp
over 1000 genomes takes several hours even at 8-way concurrency, and the reads are ~10×
larger on disk than benchmark 6's. The read-generation rule requests 24 h of queue time.

## Layout

- `genome_list.tsv` — the 1000 genomes and their FASTA paths (copy of benchmark 6's).
- `coverage_definitions/hifi1000.tsv` — per-genome coverage. Sample `hifi1000`.
- `bench14_setup.py` — tools, output dirs, local DB paths, `longread_read_tech`.
- `Snakefile` — per-benchmark variables; `include`s the shared rule library.
- `run.sh` — submit the whole benchmark to the aqua queue.

Shared rules in `../rules/` (`longread_data_generation.smk`, `staging.smk`,
`common.smk`, and one `<tool>_longread_run.smk` per profiler) provide read generation,
database staging, biobox/OPAL conversion and the per-tool profiling; see
`../ARCHITECTURE.md`.

## Running

```bash
./14_thousand_species_hifi/run.sh                 # aqua queue, one PBS job per rule
```

Running this one locally is not recommended — the read simulation alone is a multi-hour,
single-node job.
