#!/usr/bin/env python3.11

###############################################################################
#
#    Copyright (C) 2021 Ben Woodcroft
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

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2022"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
import logging
import sys
import os, re

import polars as pl
import tempfile
import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--output-condensed', required=True, help='Path to output file in singlem condensed format')
    parent_parser.add_argument('-1', '--read1', required=True, help='Path to output fq.gz file')
    parent_parser.add_argument('-2', '--read2', required=True, help='Path to output fq.gz file')
    parent_parser.add_argument('--art', required=True, help='Path to ART binary (art_illumina)')

    parent_parser.add_argument('--example-bacteria-fasta', required=True, help='Path to example bacteria genome fasta')
    parent_parser.add_argument('--example-archaea-fasta', required=True, help='Path to example archaea genome fasta')
    parent_parser.add_argument('--example-bacteria-taxonomy', required=True, help='Taxonomy of example bacteria genome')
    parent_parser.add_argument('--example-archaea-taxonomy', required=True, help='Taxonomy of example archaea genome')
    parent_parser.add_argument('--genome-fasta', required=True, help='Path to novel genome fasta file')
    parent_parser.add_argument('--taxonomy', required=True, help='Taxonomy of novel genome')
    
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    output1 = os.path.abspath(args.read1)
    output2 = os.path.abspath(args.read2)
    output_condensed = os.path.abspath(args.output_condensed)

    novel_genome_fasta = os.path.abspath(args.genome_fasta)

    read_length = 150
    coverage = 10
    sample_name = os.path.basename(args.genome_fasta).replace('.fna','')

    logging.info("Simulating 2-component community condensed file ..")
    sim_commands = []
    with open(output_condensed, 'w') as f:
        f.write("sample\tcoverage\ttaxonomy\n")

        # Write example genome
        if 'Bacteria' in args.taxonomy:
            example_genome = os.path.abspath(args.example_archaea_fasta)
            f.write(f"{sample_name}\t{coverage}\t{args.example_archaea_taxonomy}\n")
        elif 'Archaea' in args.taxonomy:
            example_genome = os.path.abspath(args.example_bacteria_fasta)
            f.write(f"{sample_name}\t{coverage}\t{args.example_bacteria_taxonomy}\n")
        else:
            raise Exception(f"Unknown taxonomy {args.taxonomy}")
        
        # Write novel genome
        f.write(f"{sample_name}\t{coverage}\t{args.taxonomy}\n")

    logging.info("Simulating reads ..")
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        os.makedirs('simulated_reads')
        # Simulate each fasta
        for i, fa in enumerate([example_genome, novel_genome_fasta]):
            logging.info("Simulating reads from %s ..", fa)
            extern.run(
                f"{args.art} -ss HSXt -i {fa} -p -l {read_length} -f {coverage} -m 400 -s 10 -o simulated_reads/{i}. &>/dev/null")

        logging.info("Concatenating simulated reads and compressing ..")
        extern.run("cat simulated_reads/*1.fq |sed 's=/= =' |pigz -p 1 >{}".format(output1))
        extern.run("cat simulated_reads/*2.fq |sed 's=/= =' |pigz -p 1 >{}".format(output2))

    logging.info("Done.")
    