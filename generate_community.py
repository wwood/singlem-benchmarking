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

import pandas as pd
import tempfile
import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--coverage-file', required=True, help='Path to coverage file')
    parent_parser.add_argument('--genome-list', required=True, help='Path to genome list')
    # (env)cl5n007:20221031:~/m/msingle/mess/115_camisim_ish_benchmarking$ \ls -f ~/m/msingle/sam/1_gtdb_r207_smpkg/20220513/shadow_GTDB/genomes |grep GC |sed 's=\(.*\).fna=\1\t/work/microbiome/msingle/sam/1_gtdb_r207_smpkg/20220513/shadow_GTDB/genomes\1.fna=' >shadow_genome_paths.csv
    parent_parser.add_argument('--gtdb-metadata', required=True, help='Path to GTDB metadata file, for genome length')
    parent_parser.add_argument('--output-condensed', required=True, help='Path to output file in singlem condensed format')
    parent_parser.add_argument('-1', '--read1', required=True, help='Path to output fq.gz file')
    parent_parser.add_argument('-2', '--read2', required=True, help='Path to output fq.gz file')
    parent_parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')

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

    # Read coverages
    coverages = pd.read_csv(args.coverage_file, sep='\t', header=None, names=['otu', 'coverage'])
    coverages = coverages[coverages['coverage'] > 0]
    logging.info(f"Read {len(coverages)} coverages > 0.")

    genomes = pd.read_csv(args.genome_list, sep='\t', header=None, names=['genome','fasta'])
    logging.info(f"Read {len(genomes)} genome fasta paths.")

    metadata = pd.read_csv(args.gtdb_metadata, sep='\t', header=0)
    logging.info(f"Read {len(metadata)} GTDB metadata entries.")

    metadata['genome'] = [g[3:] for g in metadata['accession']] # get rid of GB_, RS_

    # Shuffle genomes order so we get randomness
    metadata = metadata.sample(frac=1)

    g2 = pd.merge(genomes, metadata, on='genome', how='inner')

    # wgsim each genome and the correct coverage

    read_length = 150

    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        os.makedirs('simulated_reads')
        sim_commands = []
        with open(output_condensed, 'w') as f:
            f.write("sample\tcoverage\ttaxonomy\n")
            tax_to_coverage = {}
            for i, (fasta, genome_size, gtdb_taxonomy, coverage) in enumerate(zip(g2['fasta'], g2['genome_size'],g2['gtdb_taxonomy'], coverages['coverage'])):
                if gtdb_taxonomy not in tax_to_coverage:
                    tax_to_coverage[gtdb_taxonomy] = 0
                tax_to_coverage[gtdb_taxonomy] += coverage

                # num_reads = int(coverage*float(genome_size)/(float(read_length)*2.0))
                sim_commands.append(
                    f"~/bioinfo/art_bin_MountRainier/art_illumina -ss HSXt -i {fasta} -p -l {read_length} -f {coverage} -m 400 -s 10 -o simulated_reads/{i}. &>/dev/null"
                )

            for tax, cov in tax_to_coverage.items():
                f.write(f"{os.path.basename(args.coverage_file)}\t{cov}\t{tax}\n")
                
        logging.info(f"Simulating {len(sim_commands)} genomes ..")
        extern.run_many(sim_commands, num_threads=args.threads, progress_stream=sys.stderr)

        logging.info("Concatenating simulated reads and compressing ..")
        extern.run("cat simulated_reads/*1.fq |sed 's=/= =' |pigz -p {} >{}".format(args.threads, output1))
        extern.run("cat simulated_reads/*2.fq |sed 's=/= =' |pigz -p {} >{}".format(args.threads, output2))
    