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

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--report-prefix', required=True)
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
            splits2 = [re.sub(r'^.__', '', s) for s in splits]
            for i, name in enumerate(splits2):
                name_to_taxonomy[i][name] = ';'.join(splits[0:(i+1)])
    
    # logging.info("Read %d taxonomic names" % len(name_to_taxonomy))

    # Read braken report and convert to singlem condense format
    logging.info("Reading and printing kraken report from %s" % args.report_prefix)
    print("\t".join(['sample','coverage','taxonomy']))
    sample_name = os.path.basename(args.report_prefix)
    # print("\t".join([sample_name,'1.0','Root']))
    for tax_level, lvl in enumerate(('D','P','C','O','F','G','S')):
    # for lvl in ('P','C','O','F','G','S'):
        # name    taxonomy_id     taxonomy_lvl    kraken_assigned_reads   added_reads     new_est_reads   fraction_total_reads
        # Methanobacterium        850123699       G       237135  351     237486  0.97176
        df = pd.read_csv(args.report_prefix + '.' + lvl, sep='\t')
        for i, row in df.iterrows():
            tax = name_to_taxonomy[tax_level][row['name']]
            if row['fraction_total_reads'] > 0:
                print("\t".join([
                    sample_name,
                    str(row['fraction_total_reads']),
                    tax]))

    logging.info("Done")

