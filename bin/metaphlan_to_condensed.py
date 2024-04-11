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
import re

import extern
import pandas as pd

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--metaphlan', required=True)
    parent_parser.add_argument('--sample', required=True)
    parent_parser.add_argument('--genome-pairs')

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    print('sample\tcoverage\ttaxonomy')
    d = pd.read_csv(args.metaphlan, sep='\t', names=['clade_name', 'relative_abundance'])
    total_coverage = 0
    taxon_to_coverage = {}
    for i, row in d.iterrows():
        if 's__' in row['clade_name']:
            splits = [s.strip() for s in row['clade_name'].split(';')]
            if len(splits) != 7:
                raise Exception(f'Unexpected number of splits: {len(splits)}')
            # Remove 3 char taxons from the end of the splits
            rev = reversed(splits)
            for s in rev:
                if len(s) == 3:
                    splits.pop()
                else:
                    break
            taxons = '; '.join(splits)
            if taxons in taxon_to_coverage:
                raise Exception(f'Duplicate taxon: {taxons}')
            else:
                taxon_to_coverage[taxons] = row['relative_abundance']
            total_coverage += float(row['relative_abundance'])
    if (total_coverage < 95) & (total_coverage > 0):
        raise Exception(f'Total coverage is less than 95%: {total_coverage}')
    elif total_coverage > 101:
        raise Exception(f'Total coverage is greater than 101%: {total_coverage}')

    if total_coverage == 0:
        if args.genome_pairs:
            logging.warning(f'No coverage for {args.sample}, setting known taxonomy to 100%')
            genome_pairs = pd.read_csv(args.genome_pairs, sep='\t')
            paired_taxonomy = genome_pairs[genome_pairs['genome_ID'] == args.sample]['paired_taxonomy'].values[0]

            print(f'{args.sample}\t100.0\t{paired_taxonomy}')
        else:
            raise Exception(f'No coverage for {args.sample}')

    # Print lines in alphabetical order so less specific names go before more specific
    taxons_sorted = sorted(taxon_to_coverage.keys())
    for taxons in taxons_sorted:
        print(f'{args.sample}\t{taxon_to_coverage[taxons]}\t{taxons}')

    logging.info("Done")

