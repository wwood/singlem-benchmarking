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

    parent_parser.add_argument('--motus', required=True)
    parent_parser.add_argument('--gtdb-map', required=True)
    parent_parser.add_argument('--sample')
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

    ## git tag version 3.0.3 |  motus version 3.0.3 | map_tax 3.0.3 | gene database: nr3.0.3 | calc_mgc 3.0.3 -y insert.scaled_counts -l 75 | calc_motu 3.0.3 -k mOTU -C no_CAMI -g 3 | taxonomy: ref_mOTU_3.0.3 meta_mOTU_3.0.3
    ## call: python /home/woodcrob/e/motus-v3.0.3/bin/motus profile -f 1.fq -r 2.fq -o test.motus
    ## consensus_taxonomy     unnamed sample
    # Leptospira alexanderi [ref_mOTU_v3_00001]       0.0000000000
    # Leptospira weilii [ref_mOTU_v3_00002]   0.0000000000
    # Chryseobacterium sp. [ref_mOTU_v3_00004]        0.0000000000
    d = pd.read_csv(args.motus, sep='\t', names=['motu_full','coverage'], header=None, comment='#')
    d['motu'] = [n.split(' ')[-1].replace('[','').replace(']','') for n in d['motu_full']]

    gtdb = pd.read_csv(args.gtdb_map, sep='\t', header=None)

    gtdb2 = pd.DataFrame({
        'motu': gtdb.loc[:,0],
        'taxonomy': ['; '.join(r[1:]) for _, r in gtdb.iterrows()]
    })

    d2 = pd.merge(d[d['coverage']>0], gtdb2, on='motu')

    taxon_to_cov = {}
    for i, row in d2.iterrows():
        tax = row['taxonomy']
        if tax not in taxon_to_cov:
            taxon_to_cov[tax] = 0.
        taxon_to_cov[tax] += row['coverage']

    if len(taxon_to_cov) > 0:
        for taxon, cov in taxon_to_cov.items():
            print('sample_name\t{coverage}\t{taxon}'.format(
                coverage=cov, taxon=taxon))
    else:
        if bool(args.genome_pairs) & bool(args.sample):
            logging.warning(f'No taxa identified in {args.sample}, setting known taxonomy to 100%')
            genome_pairs = pd.read_csv(args.genome_pairs, sep='\t')
            paired_taxonomy = genome_pairs[genome_pairs['genome_ID'] == args.sample]['paired_taxonomy'].values[0]

            print(f'{args.sample}\t100.0\t{paired_taxonomy}')


    logging.info("Done")

