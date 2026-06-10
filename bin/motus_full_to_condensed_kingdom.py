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
__copyright__ = "Copyright 2024"
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
import pandas as pd

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--motus', required=True)
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

    # # git tag version 3.1.0 |  motus version 3.1.0 | map_tax 3.1.0 | gene database: nr3.1.0 | calc_mgc 3.1.0 -y insert.scaled_counts -l 75 | calc_motu 3.1.0 -k kingdom -C no_CAMI -g 3 | taxonomy: ref_mOTU_3.1.0 meta_mOTU_3.1.0
    # # call: python /home/woodcrob/e/motus-v3.1.0/bin/motus profile -s local_reads/GCA_019720975.1_genomic.1.fq.gz -k kingdom -q -t 32 -db ../tool_reference_data/motus/db_mOTU/
    # #consensus_taxonomy     unnamed sample
    # k__Bacteria     1.0000000000
    # k__Archaea      0.0000000000
    # k__Eukaryota    0.0000000000
    # k__Unknown cellular organism    0.0000000000
    # unassigned      0.0000000000
    d = pl.read_csv(args.motus, separator='\t', skip_rows=2)
    d.columns = ['kingdom', 'coverage']
    d = d.filter(pl.col('coverage')>0)
    d = d.select([
        pl.col('kingdom').str.replace('k__', 'd__'),
        'coverage',
    ])

    taxon_to_cov = {}
    for row in d.rows(named=True):
        if row['kingdom'] in ('d__Eukaryota', 'd__Unknown cellular organism', 'unassigned'):
            continue  # Skip these
        if row['kingdom'] in taxon_to_cov:
            logging.warning(f"Duplicate kingdom {row['kingdom']} in {args.motus}")
        taxon_to_cov[row['kingdom']] = row['coverage']

    if len(taxon_to_cov) > 0:
        for taxon, cov in taxon_to_cov.items():
            print('{sample}\t{coverage}\t{taxon}'.format(
                sample=args.sample, coverage=cov, taxon=taxon))
    else:
        if bool(args.genome_pairs) & bool(args.sample):
            logging.warning(f'No taxa identified in {args.sample}, setting known taxonomy to 100%')
            genome_pairs = pd.read_csv(args.genome_pairs, sep='\t')
            paired_taxonomy = genome_pairs[genome_pairs['genome_ID'] == args.sample]['paired_taxonomy'].values[0]

            print(f'{args.sample}\t1.0\t{paired_taxonomy}')

    # Add in a lower level level classification, because otherwise opal croaks. Must be non-zero, but make it tiny since all we care about is Bray-Curtis distance.
    print(f'{args.sample}\t0.000000001\td__Archaea; p__Archaeoglobi; c__Archaeoglobi; o__Archaeoglobales; f__Archaeoglobaceae; g__Archaeoglobus; s__Archaeoglobus fulgidus')

    logging.info("Done")

