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
__copyright__ = "Copyright 2021"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
import logging
import sys
import os
import pandas as pd
import polars as pl
import re

import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input', required=True)
    parent_parser.add_argument('--bacterial-taxonomy', required=True)
    parent_parser.add_argument('--archaeal-taxonomy', required=True)
    parent_parser.add_argument('--max-levels', type=int, default=7)

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Read tax, creating hash of name to fully qualified name
    name_to_taxonomy = []
    for i in range(args.max_levels):
        name_to_taxonomy.append({}) # 7 ranks

    tax_files = [args.bacterial_taxonomy, args.archaeal_taxonomy]
    for tax_file in tax_files:
        logging.info("Reading taxonomy from %s" % tax_file)
        t = pd.read_csv(tax_file, sep='\t', header=None, names=['name', 'taxonomy'])

        for taxon in t['taxonomy']:
            splits = taxon.split(';')
            for i, name in enumerate(splits):
                name_to_taxonomy[i][name] = ';'.join(splits[0:(i+1)])
    
    logging.info("Read %d taxonomic levels" % len(name_to_taxonomy))

    # 42.1946 15131   15131   no rank 0       unclassified
    # 57.8054 20729   185     no rank 1       root
    # 48.5806 17421   322     superkingdom    11484     d__Bacteria
    # 44.5176 15964   73      phylum  11485       p__Proteobacteria
    # 43.6196 15642   706     class   11486         c__Gammaproteobacteria
    # 38.7925 13911   509     order   11487           o__Enterobacterales
    # 35.2705 12648   3235    family  11488             f__Enterobacteriaceae
    # 10.0530 3605    1899    genus   347767              g__Sodalis_A
    # 2.6827  962     0       species 360862                s__Sodalis_A praecaptivus
    # 2.6827  962     962     subspecies      360863                  RS_GCF_000517425.1

    # Read braken report and convert to singlem condense format
    logging.info("Reading and printing metabuli report from %s" % args.input)
    print("\t".join(['sample','coverage','taxonomy']))
    sample_name = os.path.basename(args.input.replace('_profile.tsv',''))

    df = pl.read_csv(args.input, separator='\t')
    df.columns = ['relabund', 'id1', 'id2', 'rank', 'id3', 'taxonomy']
    
    # These profiles are 'unfilled', so need to add the subspecies to the species' totals
    last_species_taxonomy = None
    last_species_relabund = 0

    def print_last_species():
        global last_species_taxonomy, last_species_relabund
        if last_species_taxonomy:
            rank_id = ('superkingdom','phylum','class','order','family','genus','species').index('species')
            tax = name_to_taxonomy[rank_id][last_species_taxonomy.strip()]
            print("\t".join([
                sample_name,
                str(last_species_relabund),
                'Root;'+tax]))

    for row in df.rows(named=True):
        # Ignore taxons that aren't in GTDB
        if row['rank'] not in ('no rank') and row['taxonomy'].strip() not in ('d__Eukaryota','p__Chordata','c__Mammalia','o__Primates','f__Hominidae','g__Homo','s__Homo sapiens'):
            rank = row['rank']
            if rank == 'species':
                print_last_species()
                last_species_taxonomy = row['taxonomy']
                last_species_relabund = row['relabund']
            elif rank == 'subspecies':
                last_species_relabund += row['relabund']
            else:
                last_species_taxonomy = None
                rank_id = ('superkingdom','phylum','class','order','family','genus','species').index(rank)
                
                tax = name_to_taxonomy[rank_id][row['taxonomy'].strip()]
                print("\t".join([
                    sample_name,
                    str(row['relabund']),
                    'Root;'+tax]))
    print_last_species()

    logging.info("Done")

