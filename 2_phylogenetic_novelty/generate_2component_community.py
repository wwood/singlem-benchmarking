#!/usr/bin/env python3

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

def calculate_genome_length(fasta_path):
    """Calculate the length of a genome fasta file"""
    glen = 0
    with open(fasta_path) as f:
        for name, seq, qual in SeqReader().readfq(f):
            glen += len(seq)
    return glen

class SeqReader:
    # Stolen from https://github.com/lh3/readfq/blob/master/readfq.py
    def readfq(self, fp): # this is a generator function
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in fp: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:].partition(" ")[0], [], None
            for l in fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs); # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None # yield a fasta record instead
                    break

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--output-condensed', required=True, help='Path to relabund-wise output file in singlem condensed format')
    # --output-reads-wise-condensed
    parent_parser.add_argument('--output-reads-wise-condensed', required=True, help='Path to reads-wise output file in singlem condensed format')
    parent_parser.add_argument('-1', '--read1', required=True, help='Path to output fq.gz file')
    parent_parser.add_argument('-2', '--read2', required=True, help='Path to output fq.gz file')
    parent_parser.add_argument('--art', required=True, help='Path to ART binary (art_illumina)')

    parent_parser.add_argument('--example-fasta', required=True, help='Path to example genome fasta')
    parent_parser.add_argument('--example-taxonomy', required=True, help='Taxonomy of example genome')
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
    output_reads_wise_condensed = os.path.abspath(args.output_reads_wise_condensed)

    novel_genome_fasta = os.path.abspath(args.genome_fasta)

    read_length = 150
    coverage = 10
    sample_name = os.path.basename(args.genome_fasta).replace('.fna','')

    logging.info("Simulating 2-component community condensed file ..")
    sim_commands = []
    with open(output_condensed, 'w') as f:
        with open(args.output_reads_wise_condensed, 'w') as r:
            f.write("sample\tcoverage\ttaxonomy\n")
            r.write("sample\tcoverage\ttaxonomy\n")

            # Write example genome
            example_genome = os.path.abspath(args.example_fasta)
            f.write(f"{sample_name}\t{coverage}\t{args.example_taxonomy}\n")

            example_genome_length = calculate_genome_length(example_genome)
            r.write(f"{sample_name}\t{example_genome_length}\t{args.example_taxonomy}\n")
            
            # Write novel genome
            f.write(f"{sample_name}\t{coverage}\t{args.taxonomy}\n")
            novel_genome_length = calculate_genome_length(novel_genome_fasta)
            r.write(f"{sample_name}\t{novel_genome_length}\t{args.taxonomy}\n")

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
    

