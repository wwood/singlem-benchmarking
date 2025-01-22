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

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--gtdbtk-output-directory', help='GTDB-Tk output directory', required=True)
    parent_parser.add_argument('--genome-to-id', help='Mapping from genome to ID', required=True)
    parent_parser.add_argument('--genome-lengths-directory', help='Directory of genome lengths', required=True)
    parent_parser.add_argument('--read-mapping-counts', help='Read mapping counts', required=True)

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Get the genome lengths of each CAMI genome
    genome_to_length = {}
    for tsv in os.listdir(args.genome_lengths_directory):
        if 'tsv' not in tsv:
            continue
        gdf = pl.read_csv(os.path.join(args.genome_lengths_directory, tsv), separator='\t', has_header=False)
        cami_genome_name = os.path.basename(tsv).replace('.fasta.tsv','')
        genome_length = gdf['column_2'].sum()
        genome_to_length[cami_genome_name] = genome_length
    logging.info("Read in %d genome lengths" % len(genome_to_length))

    # Read genome to Id mapping
    # Otu520.0        /net/sgi/cami/data/CAMI2/nwillassen_data/marine/simulation_short_read/genomes/Filomicrobium_sp_W1_genomic.fasta
    # Otu1947 /net/sgi/cami/data/CAMI2/nwillassen_data/marine/simulation_short_read/genomes/Shewanella_sp_MR-4_genomic.fasta
    # Otu806.0        /net/sgi/cami/data/CAMI2/nwillassen_data/marine/simulation_short_read/genomes/Shewanella_violacea_DSS12_genomic.fasta
    genome_to_id_df = pl.read_csv(args.genome_to_id, separator='\t', has_header=False)
    # Want OTU to basename minus .fasta
    genome_to_id = {}
    for row in genome_to_id_df.rows():
        genome_to_id[row[0]] = os.path.basename(row[1]).replace('.fasta','')

    # (coverm-v0.7.0)cl5n006:20240501:~/m/msingle/mess/124_singlem-benchmarking/3_cami2_marine$ seq 0 9 |parallel pigz -cd ../../114_cami2_benchmarking/data/marine/short_read/simulation_short_read/'*'_sample_{}/reads/reads_mapping.tsv.gz '|'cut -f2 '|'sort '|'uniq -c '|'sed '"s/^  *//g"' '>' read_mapping_counts/marine{}.tsv
    mapping_counts = pl.read_csv(args.read_mapping_counts, separator=' ', has_header=False).filter(pl.col('column_2') != 'genome_id')
    mapping_counts.columns = ['read_count','otu']
    mapping_counts = mapping_counts.with_columns(pl.col('otu').replace(genome_to_id).alias('cami_name'))

    # Read in taxons from GTDB-Tk
    taxonomies = {}
    bac_taxonomy_file = os.path.join(args.gtdbtk_output_directory, 'gtdbtk.bac120.summary.tsv')
    if os.path.exists(bac_taxonomy_file):
        logging.info('Reading taxonomy from %s' % bac_taxonomy_file)
        d = pd.read_csv(bac_taxonomy_file, sep='\t')
        for i, row in d.iterrows():
            taxonomies[row['user_genome']] = row['classification']
    logging.info("Read %d taxonomies after Bacteria" % len(taxonomies))
    # Archaea
    arc_taxonomy_file = os.path.join(args.gtdbtk_output_directory, 'gtdbtk.ar53.summary.tsv')
    if os.path.exists(arc_taxonomy_file):
        logging.info('Reading taxonomy from %s' % arc_taxonomy_file)
        d = pd.read_csv(arc_taxonomy_file, sep='\t')
        for i, row in d.iterrows():
            taxonomies[row['user_genome']] = row['classification']
    logging.info("Read %d taxonomies after Archaea" % len(taxonomies))

    print('sample\tcoverage\ttaxonomy')
    taxon_abundances = {}
    for row in mapping_counts.rows(named=True):
        otu = row['otu']
        abund = row['read_count'] * 150 / genome_to_length[row['cami_name']]
        cami_name = row['cami_name']
        if 'RNODE' in otu:
            continue  # Ignore plasmids etc

        if cami_name in taxonomies:
            tax = taxonomies[cami_name]

            # Remove empty taxon strings because they mess up biobox creation
            tax_splits = [s.strip() for s in tax.split(';')]
            if tax_splits[-1] == 's__':
                tax_splits = tax_splits[0:-1]
            if tax_splits[-1] == 'g__':
                tax_splits = tax_splits[0:-1]
            if tax_splits[-1] == 'f__':
                tax_splits = tax_splits[0:-1]
            if tax_splits[-1] == 'o__':
                tax_splits = tax_splits[0:-1]
            if tax_splits[-1] == 'c__':
                tax_splits = tax_splits[0:-1]
            tax = ';'.join(tax_splits)
            if tax == 'Unclassified':
                raise Exception(row)
            if tax not in taxon_abundances:
                taxon_abundances[tax] = 0
            taxon_abundances[tax] += abund
        else:
            raise Exception('No taxonomy for %s' % cami_name)

    # Sort so biobox creation sees less defined taxons before more defined ones and doesn't croak.
    for tax in sorted(taxon_abundances.keys()):
        abund = taxon_abundances[tax]
        if abund > 0:
            print('%s\t%f\t%s' % (os.path.basename(args.read_mapping_counts).replace('.tsv',''), abund, 'Root;'+tax))


