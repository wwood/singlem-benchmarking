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

"""Simulate benchmark 2's two-component community as **long** reads with Badread.

The long-read counterpart of ../2_phylogenetic_novelty/generate_2component_community.py:
one novel genome plus one known ("example") genome at equal coverage, with the truth
taxonomy stated per genome on the command line rather than looked up in the GTDB
metadata -- which is the whole point of that benchmark, since the novel member is by
construction absent from the target release and so has no metadata row to join against.

Differences from the ART version, both structural:

1. **Reads are single-end** (`-r/--reads`, one fastq) rather than paired `-1`/`-2`.
2. **Only the relative-abundance truth is written.** The ART version also writes a
   reads-wise (genome-length-weighted) truth for the read-count tools kraken / kaiju /
   metabuli. The long-read benchmarks (15/16) run only singlem, sylph and
   singlem-regime3, all of which report relative abundance, so the reads-wise truth
   would be an unused file. If a read-count tool is ever added to a long-read
   benchmark, port the `--output-reads-wise-condensed` half back from the ART script.

Read technology (`--read-tech`) and its Badread parameters are shared with
generate_longread_community.py -- imported from it rather than duplicated, so the two
generators cannot drift apart on what "nanopore" or "pacbio-hifi" means.

    pixi run -e badread bin/generate_longread_2component_community.py \
        --read-tech nanopore \
        --genome-fasta 2_phylogenetic_novelty/genomes/GCA_013154095.1_genomic.fna \
        --taxonomy 'd__Archaea;p__Huberarchaeota;...' \
        --example-fasta 2_phylogenetic_novelty/genome_pairs/GCF_000364725.1_genomic.fna \
        --example-taxonomy 'd__Bacteria;p__Proteobacteria;...' \
        --output-condensed truths/GCA_013154095.1_genomic.condensed \
        -r reads/GCA_013154095.1_genomic.fq.gz
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

sys.path = [os.path.dirname(os.path.realpath(__file__))] + sys.path
from generate_longread_community import READ_TECH_PARAMS, badread_command

# Equal coverage for both members, matching the ART version's `coverage = 10`. The
# equal-abundance design is what makes the benchmark a test of *placement* rather than
# of quantification: whatever a tool does with the novel member, it cannot be excused
# by that member being rare.
COVERAGE = 10

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--read-tech', required=True, choices=sorted(READ_TECH_PARAMS.keys()),
                               help='read technology to simulate')
    parent_parser.add_argument('--output-condensed', required=True,
                               help='Path to relative-abundance truth, singlem condensed format')
    parent_parser.add_argument('-r', '--reads', required=True,
                               help='Path to output fq.gz file (long reads are single-end)')
    parent_parser.add_argument('--example-fasta', required=True, help='Path to known genome fasta')
    parent_parser.add_argument('--example-taxonomy', required=True, help='Taxonomy of known genome')
    parent_parser.add_argument('--genome-fasta', required=True, help='Path to novel genome fasta')
    parent_parser.add_argument('--taxonomy', required=True, help='Taxonomy of novel genome')
    parent_parser.add_argument('--sample', default=None,
                               help='sample name in the condensed profile '
                                    '(default: novel genome basename minus .fna, as the ART version)')
    parent_parser.add_argument('--threads', type=int, default=1,
                               help='Number of Badread processes to run concurrently')
    parent_parser.add_argument('--badread', default='badread', help='Path to the badread executable')
    parent_parser.add_argument('--seed', type=int, default=1,
                               help='Base RNG seed; the two genomes use --seed and --seed+1')

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

    # Badread runs inside a tmpdir below, so absolutise everything leaving this process.
    output_reads = os.path.abspath(args.reads)
    output_condensed = os.path.abspath(args.output_condensed)
    novel_genome_fasta = os.path.abspath(args.genome_fasta)
    example_genome_fasta = os.path.abspath(args.example_fasta)

    # Both members are named explicitly here (unlike the coverage-file generator, which
    # can simply drop a bad genome), so an empty FASTA is a hard error rather than
    # something to skip: silently simulating a one-member community would make the truth
    # wrong. Badread itself dies with ZeroDivisionError on a 0-byte reference, so check
    # for it up front and say why.
    for fasta in (novel_genome_fasta, example_genome_fasta):
        if not os.path.exists(fasta):
            raise SystemExit("genome fasta does not exist: {}".format(fasta))
        if os.path.getsize(fasta) == 0:
            raise SystemExit("genome fasta is empty (0 bytes), cannot simulate: {}".format(fasta))

    sample_name = args.sample
    if sample_name is None:
        sample_name = os.path.basename(args.genome_fasta).replace('.fna', '')

    for d in (output_reads, output_condensed):
        parent = os.path.dirname(d)
        if parent:
            os.makedirs(parent, exist_ok=True)

    # Truth: both members at equal coverage. Written in the order the ART version uses
    # (known first, then novel) so the two benchmarks' truth files are directly
    # diffable. Note the novel member's taxonomy is typically truncated (no species,
    # often no genus), which is exactly what the benchmark asks tools to reproduce.
    logging.info("Writing 2-component truth for sample %s ..", sample_name)
    with open(output_condensed, 'w') as f:
        f.write("sample\tcoverage\ttaxonomy\n")
        f.write("{}\t{}\t{}\n".format(sample_name, COVERAGE, args.example_taxonomy))
        f.write("{}\t{}\t{}\n".format(sample_name, COVERAGE, args.taxonomy))

    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        os.makedirs('simulated_reads')
        sim_commands = [
            badread_command(args.badread, fasta, COVERAGE, tech_params, args.seed + i,
                            'simulated_reads/{}.fq'.format(i))
            for i, fasta in enumerate([example_genome_fasta, novel_genome_fasta])
        ]
        logging.info("Simulating 2 genomes with Badread ..")
        extern.run_many(sim_commands, num_threads=args.threads, progress_stream=sys.stderr)

        logging.info("Concatenating simulated reads and compressing ..")
        # Badread read names are unique UUIDs with no mate suffix, so unlike the ART
        # generators there is no /1,/2 to strip.
        extern.run("cat simulated_reads/*.fq | pigz -p {} >{}".format(args.threads, output_reads))
    logging.info("Done.")
