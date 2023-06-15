#!/usr/bin/env python3

"""
Author: Samuel Aroney
Find pair-genomes (from old GTDB) from opposite domain as novel genomes (from new GTDB) with similar genome sizes

e.g.
python find_appropriate_pair.py \
    --old-gtdb-bac /work/microbiome/db/gtdb/gtdb_release207/bac120_metadata_r207.tsv \
    --old-gtdb-arc /work/microbiome/db/gtdb/gtdb_release207/ar53_metadata_r207.tsv \
    --new-gtdb-bac /work/microbiome/db/gtdb/gtdb_release214/bac120_metadata_r214.tsv \
    --new-gtdb-arc /work/microbiome/db/gtdb/gtdb_release214/ar53_metadata_r214.tsv \
    --genome-metadata stratified_genome_metadata.tsv \
    --output genome_pairs.tsv
"""

import sys
import argparse
import logging
import polars as pl

ACCESSION_REGEX = r"(\d{9})"

def pipeline(genome_metadata, old_gtdb_bac, old_gtdb_arc, new_gtdb_bac, new_gtdb_arc):
    new_genome_size = (
        pl.concat([new_gtdb_bac, new_gtdb_arc])
        .select(
            pl.col("accession").str.extract(ACCESSION_REGEX),
            "genome_size",
            )
    )

    old_genomes_bac = (
        old_gtdb_bac
        .filter(pl.col("gtdb_representative") == "t")
        .select(
            "genome_size",
            bac_accession = pl.col("accession"),
            bac_taxonomy = pl.col("gtdb_taxonomy"),
            )
        .sort("genome_size")
    )

    old_genomes_arc = (
        old_gtdb_arc
        .filter(pl.col("gtdb_representative") == "t")
        .select(
            "genome_size",
            arc_accession = pl.col("accession"),
            arc_taxonomy = pl.col("gtdb_taxonomy"),
            )
        .sort("genome_size")
    )

    paired_genomes = (
        genome_metadata
        .select(
            "genome_ID",
            accession = pl.col("genome_ID").str.extract(ACCESSION_REGEX),
            domain = pl.col("taxonomy").str.replace(";p__.*", ""),
            )
        .join(new_genome_size, on = "accession")
        .sort("genome_size")
        .join_asof(old_genomes_bac, on = "genome_size", strategy = "nearest")
        .join_asof(old_genomes_arc, on = "genome_size", strategy = "nearest")
        .select(
            "genome_ID",
            paired_genome = pl.when(pl.col("domain") == "d__Bacteria")
                .then(pl.col("arc_accession"))
                .otherwise(pl.col("bac_accession")),
            paired_taxonomy = pl.when(pl.col("domain") == "d__Bacteria")
                .then(pl.col("arc_taxonomy"))
                .otherwise(pl.col("bac_taxonomy")),
            )
        .with_columns(
            pl.col("paired_genome").str.replace("^.._", "").str.replace("$", "_genomic")
            )
    )

    return paired_genomes

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--debug", help="output debug information", action="store_true")
    parser.add_argument("--quiet", help="only output errors", action="store_true")

    parser.add_argument("--old-gtdb-bac", help="Bacterial metadata file from old GTDB version", required=True)
    parser.add_argument("--old-gtdb-arc", help="Archaeal metadata file from old GTDB version", required=True)
    parser.add_argument("--new-gtdb-bac", help="Bacterial metadata file from new GTDB version", required=True)
    parser.add_argument("--new-gtdb-arc", help="Archaeal metadata file from new GTDB version", required=True)
    parser.add_argument("--genome-metadata", help="Novel genome metadata tsv file (genome_ID, novelty_category, taxonomy)", required=True)
    parser.add_argument("--output", help="Output tsv file", required=True)

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format="%(asctime)s %(levelname)s: %(message)s", datefmt="%Y/%m/%d %I:%M:%S %p")

    metadata_cols = ["accession", "genome_size", "gtdb_representative", "gtdb_taxonomy"]
    old_gtdb_bac = pl.read_csv(args.old_gtdb_bac, separator="\t", ignore_errors=True).select(metadata_cols)
    old_gtdb_arc = pl.read_csv(args.old_gtdb_arc, separator="\t", ignore_errors=True).select(metadata_cols)
    new_gtdb_bac = pl.read_csv(args.new_gtdb_bac, separator="\t", ignore_errors=True).select(metadata_cols)
    new_gtdb_arc = pl.read_csv(args.new_gtdb_arc, separator="\t", ignore_errors=True).select(metadata_cols)

    genome_metadata = pl.read_csv(args.genome_metadata, separator="\t", ignore_errors=True)

    paired_genomes = pipeline(genome_metadata, old_gtdb_bac, old_gtdb_arc, new_gtdb_bac, new_gtdb_arc)
    paired_genomes.write_csv(args.output, separator="\t")

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
