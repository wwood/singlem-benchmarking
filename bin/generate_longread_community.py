#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2026 Ben Woodcroft
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

"""Simulate a synthetic community as **long** reads with Badread.

This is the long-read counterpart of ../1_novel_strains/generate_community.py: same
inputs (a per-genome coverage file plus a genome->fasta list, with genome sizes and
truth taxonomy taken from the GTDB metadata by inner join), same outputs (a single
reads file, a condensed ground-truth profile and a per-genome coverage record), but
reads come from Badread rather than art_illumina.

Two differences from the ART generator are structural rather than incidental:

1. **Reads are single-end.** Long reads are not paired, so this writes ONE fastq
   (`-r/--reads`) where the ART generator writes `-1`/`-2`. Every downstream rule for
   these benchmarks therefore takes a single read file.

2. **Coverage is requested per genome via Badread's `--quantity <N>x`.** Badread
   computes the target base count from the reference it is given, so the coverage
   column is honoured per genome without needing the genome length up front. (The
   GTDB `genome_size` is still read, because the truth taxonomy comes from the same
   metadata join and keeping the join identical to the ART generator's is what makes
   the two benchmark families comparable.)

Read technology is selected with `--read-tech`:

* `nanopore` - Badread's defaults, which the Badread docs state "correspond to Oxford
  Nanopore R10.4.1 reads of mediocre quality": `--error_model nanopore2023
  --qscore_model nanopore2023 --identity 95,99,2.5 --length 15000,13000`. Passed
  explicitly rather than relied upon so a future Badread release changing its
  defaults cannot silently change what this benchmark simulates.
* `pacbio-hifi` - `--error_model pacbio2021 --qscore_model pacbio2021 --identity
  30,3`, the HiFi recipe from the Badread docs. The two-value `--identity` switches
  Badread from its beta identity distribution to sampling qscores from a normal
  distribution (mean 30, stdev 3), which is what makes the reads HiFi-accurate
  (~99.9%) rather than merely long.

Badread is single-threaded, so genomes are simulated concurrently with
`extern.run_many` (--threads controls how many at once), as the ART generators do.
Each genome gets a distinct `--seed` derived from --seed so a rerun reproduces the
same reads.

    pixi run -e badread bin/generate_longread_community.py \
        --read-tech nanopore --coverage-file coverage_definitions/nanopore2.tsv \
        --genome-list genome_list.tsv \
        --gtdb-bac-metadata ../bac120_metadata_r207.tsv \
        --gtdb-ar-metadata ../ar53_metadata_r207.tsv \
        --output-condensed truths/nanopore2.condensed \
        --output-genomewise-coverage truths/nanopore2.genomewise.csv \
        -r reads/nanopore2.fq.gz
"""

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2026"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
import logging
import os
import sys
import tempfile

import extern
import polars as pl

# Badread parameters per read technology. Both entries are the Badread
# documentation's own recommended settings for that chemistry; see the module
# docstring. Stated explicitly (not left to Badread's defaults) so the simulated
# data is pinned to these values regardless of the Badread version installed.
READ_TECH_PARAMS = {
    'nanopore': {
        'error_model': 'nanopore2023',
        'qscore_model': 'nanopore2023',
        'identity': '95,99,2.5',
        'length': '15000,13000',
    },
    'pacbio-hifi': {
        'error_model': 'pacbio2021',
        'qscore_model': 'pacbio2021',
        # Two values => normal qscore distribution (mean,stdev), Badread's
        # recommended form for high-accuracy reads such as HiFi.
        'identity': '30,3',
        'length': '15000,13000',
    },
}


def badread_command(badread, fasta, coverage, tech_params, seed, out_fastq):
    """A single Badread invocation, writing uncompressed fastq to out_fastq."""
    return (
        "{badread} simulate --reference {fasta} --quantity {coverage}x "
        "--error_model {error_model} --qscore_model {qscore_model} "
        "--identity {identity} --length {length} --seed {seed} "
        "> {out} 2>/dev/null".format(
            badread=badread,
            fasta=fasta,
            coverage=coverage,
            seed=seed,
            out=out_fastq,
            **tech_params))


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--read-tech', required=True, choices=sorted(READ_TECH_PARAMS.keys()),
                               help='read technology to simulate')
    parent_parser.add_argument('--coverage-file', required=True, help='Path to coverage file')
    parent_parser.add_argument('--genome-list', required=True, help='Path to genome list')
    parent_parser.add_argument('--gtdb-bac-metadata', required=True,
                               help='Path to GTDB bacterial metadata file, for taxonomy')
    parent_parser.add_argument('--gtdb-ar-metadata', required=True,
                               help='Path to GTDB archaeal metadata file, for taxonomy')
    parent_parser.add_argument('--output-condensed', required=True,
                               help='Path to output file in singlem condensed format')
    parent_parser.add_argument('--output-genomewise-coverage', required=True,
                               help='Path to output file specifying coverage for each strain')
    parent_parser.add_argument('-r', '--reads', required=True,
                               help='Path to output fq.gz file (long reads are single-end)')
    parent_parser.add_argument('--threads', type=int, default=1,
                               help='Number of Badread processes to run concurrently')
    parent_parser.add_argument('--badread', default='badread', help='Path to the badread executable')
    parent_parser.add_argument('--seed', type=int, default=1,
                               help='Base RNG seed; genome i is simulated with --seed (seed + i)')

    args = parent_parser.parse_args()

    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    tech_params = READ_TECH_PARAMS[args.read_tech]
    logging.info("Simulating %s reads with Badread parameters: %s", args.read_tech, tech_params)

    # Badread and the concatenation step run inside a tmpdir below, so every path
    # leaving this process must be absolute first.
    output_reads = os.path.abspath(args.reads)
    output_condensed = os.path.abspath(args.output_condensed)
    output_genomewise = os.path.abspath(args.output_genomewise_coverage)

    # Coverages: same handling as the ART generator, including dropping the RNODE
    # (plasmid etc.) rows that the CAMI-derived coverage definitions can contain.
    coverages = pl.read_csv(args.coverage_file, separator='\t', has_header=False,
                            new_columns=['otu', 'coverage'])
    coverages = coverages.filter(pl.col('coverage') > 0)
    logging.info("Read %d coverages > 0.", coverages.height)
    coverages = coverages.filter(~pl.col('otu').str.contains('RNODE'))
    logging.info("After removing plasmids etc, %d coverages > 0 remain.", coverages.height)

    genomes = pl.read_csv(args.genome_list, separator='\t', has_header=False,
                          new_columns=['genome', 'fasta'])
    logging.info("Read %d genome fasta paths.", genomes.height)

    # Genome sizes and taxonomy come from the GTDB metadata, so the truth is written
    # in whatever GTDB release the caller targets. Identical join to the ART
    # generator's, which is what keeps the long-read benchmarks comparable to their
    # short-read counterparts (and, as there, means a genome absent from the target
    # release is dropped from both the simulation and the truth).
    bac = pl.read_csv(args.gtdb_bac_metadata, separator='\t', infer_schema_length=100000,
                      ignore_errors=True)
    ar = pl.read_csv(args.gtdb_ar_metadata, separator='\t', infer_schema_length=100000,
                     ignore_errors=True)
    metadata = pl.concat([
        bac.select('accession', 'genome_size', 'gtdb_taxonomy'),
        ar.select('accession', 'genome_size', 'gtdb_taxonomy'),
    ])
    logging.info("Read %d GTDB metadata entries.", metadata.height)

    # Strip the GB_/RS_ prefix to match the genome list's accessions.
    metadata = metadata.with_columns(pl.col('accession').str.slice(3).alias('genome'))

    # Shuffle so the coverage <-> genome pairing is randomised, as the ART generator
    # does: the coverage column is positional, so which genome draws which coverage
    # is arbitrary by design in these benchmarks. (A benchmark where the pairing is
    # the experiment must state it per genome instead -- see 10_congeneric_novelty.)
    metadata = metadata.sample(fraction=1.0, shuffle=True, seed=args.seed)

    g2 = genomes.join(metadata, on='genome', how='inner')
    if g2.height < genomes.height:
        logging.warning("%d of %d genomes had no GTDB metadata row and were dropped.",
                        genomes.height - g2.height, genomes.height)

    g2 = g2.with_columns(
        pl.col('fasta').map_elements(os.path.abspath, return_dtype=pl.String))

    # The coverage column is consumed positionally against the (shuffled) genome
    # table, so only as many genomes as there are coverages are simulated.
    n = min(g2.height, coverages.height)
    if n < coverages.height:
        logging.warning("Only %d genomes available for %d coverages; simulating %d.",
                        g2.height, coverages.height, n)

    for d in (output_reads, output_condensed, output_genomewise):
        parent = os.path.dirname(d)
        if parent:
            os.makedirs(parent, exist_ok=True)

    accessions = g2['accession'].to_list()[:n]
    fastas = g2['fasta'].to_list()[:n]
    taxonomies = g2['gtdb_taxonomy'].to_list()[:n]
    cov_values = coverages['coverage'].to_list()[:n]

    for fasta in fastas:
        if not os.path.exists(fasta):
            raise SystemExit("genome fasta does not exist: {}".format(fasta))

    # Drop genomes whose FASTA is empty. A handful of the shadow-genome FASTAs are
    # 0 bytes (failed downloads), and Badread cannot simulate from an empty reference:
    # it divides the progress counter by the target base count and dies with
    # ZeroDivisionError. They must also be dropped from the *truth*, not just from the
    # simulation -- otherwise the truth asserts species that produced no reads, which
    # every tool would then be scored as having missed. (ART does not fail on an empty
    # reference, so benchmark 6 keeps such genomes in its truth; that is why its truth
    # has 1000 lineages where this generator writes fewer.)
    keep = [i for i, fasta in enumerate(fastas) if os.path.getsize(fasta) > 0]
    if len(keep) < len(fastas):
        dropped = [(accessions[i], fastas[i]) for i in range(len(fastas)) if i not in set(keep)]
        logging.warning(
            "Dropping %d genome(s) with a 0-byte FASTA (excluded from both the reads and "
            "the truth): %s", len(dropped), ", ".join(a for a, _ in dropped))
        accessions = [accessions[i] for i in keep]
        fastas = [fastas[i] for i in keep]
        taxonomies = [taxonomies[i] for i in keep]
        cov_values = [cov_values[i] for i in keep]

    # Truth is per lineage, summing genomes that share one (as the ART generator does).
    tax_to_coverage = {}
    with open(output_genomewise, 'w') as genome_wise_f:
        genome_wise_f.write("accession\tcoverage\tfasta\ttaxonomy\n")
        for accession, fasta, taxonomy, coverage in zip(accessions, fastas, taxonomies, cov_values):
            tax_to_coverage[taxonomy] = tax_to_coverage.get(taxonomy, 0) + coverage
            genome_wise_f.write("{}\t{}\t{}\t{}\n".format(accession, coverage, fasta, taxonomy))

    with open(output_condensed, 'w') as f:
        f.write("sample\tcoverage\ttaxonomy\n")
        for taxonomy, coverage in tax_to_coverage.items():
            f.write("{}\t{}\t{}\n".format(os.path.basename(args.coverage_file), coverage, taxonomy))
    logging.info("Wrote truth (%d lineages) to %s", len(tax_to_coverage), output_condensed)

    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        os.makedirs('simulated_reads')
        sim_commands = [
            badread_command(args.badread, fasta, coverage, tech_params, args.seed + i,
                            'simulated_reads/{}.fq'.format(i))
            for i, (fasta, coverage) in enumerate(zip(fastas, cov_values))
        ]
        logging.info("Simulating %d genomes with Badread (%d at a time) ..",
                     len(sim_commands), args.threads)
        extern.run_many(sim_commands, num_threads=args.threads, progress_stream=sys.stderr)

        logging.info("Concatenating simulated reads and compressing ..")
        # Badread read names are already unique per read (a UUID) and carry no mate
        # suffix, so unlike the ART generators there is no /1,/2 to strip. Reads from
        # different genomes can share a name only if two Badread runs drew the same
        # UUID, which does not happen in practice.
        extern.run("cat simulated_reads/*.fq | pigz -p {} >{}".format(args.threads, output_reads))
    logging.info("Done.")
