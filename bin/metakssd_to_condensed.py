#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2021 Ben Woodcroft, 2025 Tim Lamberton
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

__author__ = "Tim Lamberton"
__copyright__ = "Copyright 2025"
__credits__ = ["Ben Woodcroft", "Tim Lamberton"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
import logging
import sys
import os
import polars as pl

import pdb

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--metakssd-genome', required=True)
    parent_parser.add_argument('--sample', required=True)
    parent_parser.add_argument('--bac-tax', required=True)
    parent_parser.add_argument('--arc-tax', required=True)

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

    SPECIES_REGEX = r';s__([^;]+)'

    taxonomy = (
        pl.concat([
            pl.read_csv(args.bac_tax, separator='\t', has_header=False, new_columns=['genome', 'taxonomy']),
            pl.read_csv(args.arc_tax, separator='\t', has_header=False, new_columns=['genome', 'taxonomy'])
            ])
        .select(pl.col('taxonomy'))
        .unique()
        .with_columns(
            pl.col('taxonomy').str.extract(SPECIES_REGEX).alias('species')
            )
    )

    PSID_REGEX = r'^(\d+)_'

    metakssd_genome = (
        pl.read_csv(args.metakssd_genome, separator="\t", has_header=False, new_columns=['sample', 'psid_species', 'coverage'])
        .with_columns(
            pl.col("psid_species").str.extract(PSID_REGEX).alias("psid"),
            pl.col("psid_species").str.replace(PSID_REGEX, '').alias("species"),
        )
        .filter(pl.col("species").is_not_null())
        .join(
            taxonomy.with_columns(
                pl.col("species").str.replace(' +', '_')
            ),
            on="species",
            how="left"
        )
        .select([
            pl.lit(args.sample).alias("sample"),
            pl.col("coverage"),
            pl.col("taxonomy"),
            ])
    )

    output = metakssd_genome.write_csv(include_header=False, separator='\t')
    print(output, end="")

    logging.info("Done")

