#! /usr/bin/env python


import argparse
import logging
import sys
import os, re

import polars as pl
import tempfile
import extern



# %%
def read_gtdbtk(output_directory, taxonomy_only=True, remove_empty_ranks=False):
    if not taxonomy_only:
        raise NotImplementedError("Only taxonomy is supported for now")
    taxonomies = {}
    bac_taxonomy_file = os.path.join(output_directory, 'gtdbtk.bac120.summary.tsv')
    logging.debug('Reading taxonomy from %s' % bac_taxonomy_file)
    d = pl.read_csv(bac_taxonomy_file, separator='\t')
    if remove_empty_ranks:
        empty_ranks = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    for row in d.rows(named=True):
        tax = row['classification']
        if remove_empty_ranks:
            tax = ';'.join([x for x in tax.split(';') if x.strip() not in empty_ranks])
        taxonomies[row['user_genome']] = tax
    logging.debug("Read %d taxonomies from Bacteria" % len(taxonomies))

    # Archaea
    arc_taxonomy_file = os.path.join(output_directory, 'gtdbtk.ar53.summary.tsv')
    logging.debug('Reading taxonomy from %s' % arc_taxonomy_file)
    d = pl.read_csv(arc_taxonomy_file, separator='\t')
    num_archaea = 0
    for row in d.rows(named=True):
        tax = row['classification']
        if remove_empty_ranks:
            tax = ';'.join([x for x in tax.split(';') if x.strip() not in empty_ranks])
        taxonomies[row['user_genome']] = tax
        num_archaea += 1
    logging.debug("Read %d new archaeal taxonomies, so %d total" % (num_archaea, len(taxonomies)))
    return taxonomies


# %%
if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--coverage-file', required=True, help='Path to coverage file')
    parent_parser.add_argument('--known-genome-list', required=True, help='Path to genome list')
    # (env)cl5n007:20221031:~/m/msingle/mess/115_camisim_ish_benchmarking$ \ls -f ~/m/msingle/sam/1_gtdb_r207_smpkg/20220513/shadow_GTDB/genomes |grep GC |sed 's=\(.*\).fna=\1\t/work/microbiome/msingle/sam/1_gtdb_r207_smpkg/20220513/shadow_GTDB/genomes\1.fna=' >shadow_genome_paths.csv
    parent_parser.add_argument('--novel-genome-gtdbtk-output', required=True, help='Path to where new genomes have been run through gtdbtk')
    parent_parser.add_argument('--novel-genome-list', required=True, help='Path to genome list')
    parent_parser.add_argument('--percent-known', required=True, help='Percent of known genomes to use (percent, not fraction)')
    parent_parser.add_argument('--gtdb-bac-metadata', required=True, help='Path to GTDB metadata file, for genome length')
    parent_parser.add_argument('--gtdb-ar-metadata', required=True, help='Path to GTDB metadata file, for genome length')
    parent_parser.add_argument('--output-condensed', required=True, help='Path to output file in singlem condensed format')
    parent_parser.add_argument('-1', '--read1', required=True, help='Path to output fq.gz file')
    parent_parser.add_argument('-2', '--read2', required=True, help='Path to output fq.gz file')
    parent_parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    parent_parser.add_argument('--art', required=True, help='Path to ART binary (art_illumina)')

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


    # %%
    # class Args:
    #     read1 = 'r1.fq.gz'
    #     read2 = 'r2.fq.gz'
    #     coverage_file = 'coverage_definitions/coverage0.tsv'
    #     known_genome_list = 'shadow_genome_paths.csv'
    #     novel_genome_gtdbtk_output = 'gtdbtk_batchfile.random1000.gtdbtk_r207'
    #     novel_genome_list = 'gtdbtk_batchfile.random1000.csv'
    #     gtdb_metadata = '/work/microbiome/db/gtdb/gtdb_release207/bac120_metadata_r207.tsv'
    #     output_condensed = 'output_condensed.tsv'
    #     threads = 1
    #     art = 'art_illumina'
    # args = Args()

    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # %%

    output1 = os.path.abspath(args.read1)
    output2 = os.path.abspath(args.read2)
    output_condensed = os.path.abspath(args.output_condensed)


# %%

    # Read coverages
    coverages = pl.read_csv(args.coverage_file, separator='\t')
    coverages.columns = ['otu', 'coverage']
    coverages = coverages.filter(pl.col('coverage') > 0)
    logging.info(f"Read {len(coverages)} coverages > 0.")


# %%
    # Remove RNODE ones which are plasmids etc.
    coverages = coverages.filter(~pl.col('otu').str.contains('RNODE'))
    logging.info(f"After removing plasmids etc, {len(coverages)} coverages > 0 remain.")


# %%

    genomes = pl.read_csv(args.known_genome_list, separator='\t', has_header=False)
    genomes.columns = ['genome', 'path']
    logging.info(f"Read {len(genomes)} genome fasta paths.")


# %%
    bac = pl.read_csv('../bac120_metadata_r207.tsv', separator='\t', infer_schema_length=100000, ignore_errors=True)
    ar = pl.read_csv('../ar53_metadata_r207.tsv', separator='\t', infer_schema_length=100000, ignore_errors=True)
    metadata = pl.concat([
        bac.select('accession', 'genome_size', 'gtdb_taxonomy'),
        ar.select('accession', 'genome_size', 'gtdb_taxonomy'),
    ])
    logging.info(f"Read {len(metadata)} GTDB metadata entries.")
    # metadata[:3]


# %%

    # get rid of GB_, RS_
    metadata = metadata.with_columns(pl.col('accession').str.slice(3).alias('genome'))
    metadata[:3]

# %%

    # Shuffle genomes order so we get randomness
    metadata = metadata.sample(fraction=1)

    known_info = genomes.join(metadata, on='genome', how='inner').select('path','genome',pl.col('gtdb_taxonomy').alias('taxonomy'))
    known_info[:3]

    r207_taxonomy = read_gtdbtk(args.novel_genome_gtdbtk_output, remove_empty_ranks=True)

# %%

    # Read list of genome paths
    novel_genome_list = pl.read_csv(args.novel_genome_list, separator='\t', has_header=False)
    novel_genome_list.columns = ['fasta', 'genome']


# %%

    # Merge with GTDBTK output
    # novel_genome_list.filter(~(pl.col('genome').is_in(r207_taxonomy.keys())))
    novel_info = novel_genome_list.with_columns(pl.col('genome').map_dict(r207_taxonomy).alias('taxonomy'))
    novel_info = novel_info.rename({'fasta': 'path'})
    novel_info[:3]
    # r207_taxonomy['GB_GCA_023301765.1']


# %%

    logging.info(f"Read {len(novel_info)} novel genomes.")

    # Choose 50% of the coverages to be from the new genomes, 50% from the known genomes
    fraction_new = float(args.percent_known) / 100
    n_new = round(len(coverages) * fraction_new)
    n_known = len(coverages) - n_new

    logging.info(f"Choosing {n_new} novel genomes and {n_known} known genomes.")
    chosen_df = pl.concat(
        [known_info.sample(n_known),
        novel_info.sample(n_new)])
    
    # Add coverage column
    chosen_df = chosen_df.with_columns(pl.lit(coverages.select(pl.col('coverage').shuffle())['coverage']).alias('coverage'))
    
    chosen_df.shape, chosen_df[:3], chosen_df[-3:]

    # Make paths absolute so that they work in a tempdir
    chosen_df = chosen_df.with_columns(pl.col('path').map_elements(lambda x: os.path.abspath(x)))

    read_length = 150

    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        os.makedirs('simulated_reads')
        sim_commands = []
        with open(output_condensed, 'w') as f:
            f.write("sample\tcoverage\ttaxonomy\n")
            tax_to_coverage = {}

            for i, (fasta, genome_id, taxonomy, coverage) in enumerate(chosen_df.rows()):
                if taxonomy not in tax_to_coverage:
                    tax_to_coverage[taxonomy] = 0
                tax_to_coverage[taxonomy] += coverage

                sim_commands.append(
                    f"{args.art} -ss HSXt -i {fasta} -p -l {read_length} -f {coverage} -m 400 -s 10 -o simulated_reads/{i}. &>/dev/null"
                )

            for tax, cov in tax_to_coverage.items():
                f.write(f"{os.path.basename(args.coverage_file)}\t{cov}\t{tax}\n")
                
        logging.info(f"Simulating {len(sim_commands)} genomes ..")
        extern.run_many(sim_commands, num_threads=args.threads, progress_stream=sys.stderr)

        logging.info("Concatenating simulated reads and compressing ..")
        extern.run("cat simulated_reads/*1.fq |sed 's=/= =' |pigz -p {} >{}".format(args.threads, output1))
        extern.run("cat simulated_reads/*2.fq |sed 's=/= =' |pigz -p {} >{}".format(args.threads, output2))
    


