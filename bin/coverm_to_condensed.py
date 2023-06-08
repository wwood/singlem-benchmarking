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

    parent_parser.add_argument('--coverm-genome', required=True)
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

    NCBI_GENOME_ID_REGEX = r'GC[AF]_([0-9]{9})'

    taxonomy = (
        pl.concat([
            pl.read_csv(args.bac_tax, separator='\t', has_header=False, new_columns=['genome', 'taxonomy']),
            pl.read_csv(args.arc_tax, separator='\t', has_header=False, new_columns=['genome', 'taxonomy'])
            ])
        .with_columns(
            pl.col('genome').str.extract(NCBI_GENOME_ID_REGEX).alias('genome')
            )
    )

    coverm_genome = (
        pl.read_csv(args.coverm_genome, separator='\t', new_columns=['genome', 'trimmed_mean'], dtypes=[pl.Utf8, pl.Float64])
        .filter(pl.col('trimmed_mean') != 0)
        .with_columns(pl.col('genome').str.extract(NCBI_GENOME_ID_REGEX))
        .filter(pl.col("genome").is_not_null())
        .join(taxonomy, on="genome", how="left")
        .groupby("taxonomy")
        .agg(pl.sum("trimmed_mean"))
        .select([
            pl.lit(args.sample).alias("sample"),
            pl.col("trimmed_mean").alias("coverage"),
            pl.col("taxonomy"),
            ])
    )

    output = coverm_genome.write_csv(has_header=False, separator='\t')
    print(output, end="")

    logging.info("Done")

