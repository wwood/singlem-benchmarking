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

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--summary-csv', required=True)
    parent_parser.add_argument('--sample', required=True)

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
    df = pl.read_csv(args.summary_csv, has_header=True)
    # query_name,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank,query_ani_at_rank,total_weighted_hashes                                                                            G
    # GCA_020052375.1_genomic,superkingdom,0.281,d__Archaea,3189a479,../../local_reads/GCA_020052375.1_genomic.2.fq.gz,0.456,1619000,0.960,25843                                                                       
    # GCA_020052375.1_genomic,superkingdom,0.719,unclassified,3189a479,../../local_reads/GCA_020052375.1_genomic.2.fq.gz,0.544,4148000,,25843                                                                          G
    # GCA_020052375.1_genomic,phylum,0.281,d__Archaea;p__Thermoproteota,3189a479,../../local_reads/GCA_020052375.1_genomic.2.fq.gz,0.456,1619000,0.960,25843                                                           G
    # ..
    for row in df.rows(named=True):
        if row['lineage'] != 'unclassified' and row['rank'] in ['superkingdom','phylum', 'class', 'order', 'family', 'genus', 'species']:
            print('{sample}\t{coverage}\t{taxon}'.format(
                sample=args.sample, coverage=row['f_weighted_at_rank'], taxon=row['lineage']))
            
    logging.info("Done")

