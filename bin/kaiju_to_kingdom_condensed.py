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
import polars as pl
import re

import pytaxonkit

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    # shell:
    # "../bin/kaiju_to_condensed.py {input} " \
    # "--sample {wildcards.sample} " \
    # " > {output} "
    parent_parser.add_argument('--input', required=True)
    parent_parser.add_argument('--sample', required=True)
    # data-dir
    parent_parser.add_argument('--data-dir', required=True)

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Read input
    logging.info("Reading input from %s" % args.input)
    df = pl.read_csv(args.input, separator='\t', has_header=True, null_values='NA')

    # Get taxid and fraction of reads as dict
    taxid_to_fraction_assigned = {}
    for taxid, fraction in zip(df['taxon_id'].to_list(), df['percent'].to_list()):
        taxid_to_fraction_assigned[taxid] = fraction

    # Get kingdom for each taxid with pytaxonkit
    logging.info("Getting kingdom for each taxid")
    taxid_to_kingdom = {}
    pdf = pytaxonkit.lineage([taxid for taxid in taxid_to_fraction_assigned], data_dir=args.data_dir, formatstr="{k}")
    # iterate pandas by row
    for row in pdf.itertuples():
        taxid_to_kingdom[row.TaxID] = row.Lineage

    # Write out condensed file
    print("\t".join(['sample','coverage','taxonomy']))

    sample_name = args.sample

    totals = {
        'Bacteria': 0,
        'Archaea': 0,
    }
    for taxid, fraction in taxid_to_fraction_assigned.items():
        if taxid in taxid_to_kingdom:
            kingdom = taxid_to_kingdom[taxid]
            if kingdom in totals: # e.g. remove viruses
                totals[kingdom] += fraction
    total_coverage = sum(totals.values())
    for kingdom, fraction in totals.items():
        # Use follu defined taxonomies, even though one will be wrong after the kingdom level, because otherwise OPAL croaks
        if kingdom == 'Bacteria':
            taxonomy = 'd__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus aureus'
        elif kingdom == 'Archaea':
            taxonomy = 'd__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter ruminantium'
            
        print("\t".join([sample_name, str(fraction/total_coverage), taxonomy]))
    
    logging.info("Done")

