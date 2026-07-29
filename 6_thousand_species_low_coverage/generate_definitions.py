#!/usr/bin/env python3
"""Pick benchmark 6's community: 1000 distinct GTDB r207 species, every one below
1x coverage and below 10% relative abundance.

Writes genome_list.tsv and coverage_definitions/known1000.tsv. Run once from this
directory; both outputs are committed, so the benchmark does not depend on
re-running it. Seeded, so re-running reproduces the same community.

    cd 6_thousand_species_low_coverage && pixi run -e art python3 generate_definitions.py

Genomes come from ../reference_genomes/shadow, the pool of non-representative
GTDB r207 strains the other benchmarks draw on. Each selected genome's *species*
is in the r207 tool databases (the species' representative is a different strain),
so this is a known-species benchmark in the same sense as benchmark 5: every
species is findable, and what is under test is whether it is found at coverages
below the marker-gene floor.
"""

import argparse
import os

import numpy as np
import polars as pl

NUM_SPECIES = 1000
MAX_COVERAGE = 1.0        # hard cap: no genome may exceed 1x
MAX_RELATIVE_ABUNDANCE = 0.10
MIN_COVERAGE = 0.02       # floor, so no genome is simulated at an absurd depth
SIGMA = 0.5               # lognormal spread; see the note on the range below
SEED = 1000

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--genome-directory', default='../reference_genomes/shadow')
parser.add_argument('--bac-metadata', default='../bac120_metadata_r207.tsv')
parser.add_argument('--ar-metadata', default='../ar53_metadata_r207.tsv')
parser.add_argument('--output-genome-list', default='genome_list.tsv')
parser.add_argument('--output-coverage', default='coverage_definitions/known1000.tsv')
args = parser.parse_args()

accessions = sorted(f[:-len('.fna')] for f in os.listdir(args.genome_directory)
                    if f.endswith('.fna'))
print("{} genomes in {}".format(len(accessions), args.genome_directory))

metadata = pl.concat([
    pl.read_csv(f, separator='\t', infer_schema_length=100000, ignore_errors=True)
      .select('accession', 'genome_size', 'gtdb_taxonomy')
    for f in (args.bac_metadata, args.ar_metadata)
]).with_columns(
    # GTDB metadata accessions carry a GB_/RS_ prefix the FASTA names do not.
    pl.col('accession').str.slice(3).alias('genome'),
    pl.col('gtdb_taxonomy').str.split(';').list.last().alias('species'),
).filter(pl.col('genome').is_in(accessions))
print("{} matched GTDB r207 metadata, {} distinct species".format(
    metadata.height, metadata.select('species').n_unique()))

# One genome per species: the taxonomic profile is per-species, so two strains of
# one species would have their coverages summed into a single truth row and the
# per-species caps would no longer be the per-genome caps. Sorted then seeded-
# sampled rather than taking the first of each group, so the choice within a
# species is not an artefact of accession ordering.
rng = np.random.default_rng(SEED)
metadata = metadata.sort('genome')
chosen = (metadata
          .with_columns(pl.Series('r', rng.random(metadata.height)))
          .sort('r')
          .unique(subset=['species'], keep='first')
          .sort('genome'))
if chosen.height < NUM_SPECIES:
    raise SystemExit("only {} distinct species available, need {}".format(
        chosen.height, NUM_SPECIES))
chosen = chosen.sort(pl.Series('r2', rng.random(chosen.height))).head(NUM_SPECIES).sort('genome')

# Coverages: lognormal, rescaled so the most abundant genome sits exactly at the
# 1x cap. sigma=0.5 over 1000 draws spans roughly 25-fold, so the community runs
# from ~0.04x to 1x with a median near 0.2x -- skewed like a real community, but
# with every member below SingleM's ~1.5x marker-gene sensitivity, which is the
# regime the benchmark exists to probe.
coverages = rng.lognormal(mean=0.0, sigma=SIGMA, size=NUM_SPECIES)
coverages = coverages / coverages.max() * MAX_COVERAGE
coverages = np.maximum(coverages, MIN_COVERAGE)
coverages = np.round(coverages, 4)

relative = coverages / coverages.sum()
assert coverages.max() <= MAX_COVERAGE + 1e-9, coverages.max()
assert relative.max() <= MAX_RELATIVE_ABUNDANCE, relative.max()

os.makedirs(os.path.dirname(args.output_coverage) or '.', exist_ok=True)
with open(args.output_genome_list, 'w') as genome_list, open(args.output_coverage, 'w') as coverage_file:
    for genome, coverage in zip(chosen['genome'], coverages):
        genome_list.write("{}\t{}/{}.fna\n".format(genome, args.genome_directory, genome))
        coverage_file.write("{}\t{}\n".format(genome, coverage))

sizes = chosen['genome_size'].to_numpy()
print("wrote {} genomes to {} and {}".format(NUM_SPECIES, args.output_genome_list, args.output_coverage))
print("coverage:  min {:.4f}  median {:.4f}  max {:.4f}  total {:.1f}".format(
    coverages.min(), np.median(coverages), coverages.max(), coverages.sum()))
print("relative abundance: max {:.4%}  (cap {:.0%})".format(relative.max(), MAX_RELATIVE_ABUNDANCE))
print("genome size: median {:.2f} Mbp".format(np.median(sizes) / 1e6))
# generate_community.py pairs coverages with genomes positionally against a
# shuffled metadata merge, so which genome receives which coverage is not the
# pairing written here -- only the multiset of coverages carries over. Both caps
# are properties of that multiset (one species per genome), so they hold either
# way, but the sequence total is an expectation rather than an exact figure.
print("simulated sequence: {:.2f} Gbp (expected; genome/coverage pairing is randomised)".format(
    float(coverages.mean() * sizes.sum()) / 1e9))
