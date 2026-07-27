
##### Constants related to where the reference data is for each tool.

from os.path import join, abspath

output_directory = abspath('../tool_reference_data')

singlem_metapackage = join(output_directory, 'S4.1.0.GTDB_r207.metapackage_20240502.smpkg')
singlem_metapackage_tgz = singlem_metapackage + '.zb.tar.gz'

metaphlan_db = join(output_directory, 'metaphlan_bowtiedb')
metaphlan42_db = join(output_directory, 'metaphlan42_bowtiedb')
metaphlan_index = 'mpa_vOct22_CHOCOPhlAnSGB_202212'

motus_db = join(output_directory, 'motus', 'db_mOTU')
motus_gtdb_tsv = join(output_directory, 'motus', 'mOTUs_3.0.0_GTDB_tax.tsv')

kraken_db = join(output_directory, "struo2-kraken-GTDB_release207")

sourmash_db_taxonomy = join(output_directory, 'sourmash', 'gtdb-rs207.taxonomy.sqldb')
sourmash_db_dna = join(output_directory, 'sourmash', 'gtdb-rs207.genomic-reps.dna.k31.zip')

kaiju_db_original_dir = join(output_directory, 'kaiju')
kaiju_db_original_progenomes_fmi = kaiju_db_original_dir + '/kaiju_db_progenomes.fmi'
kaiju_db_original_progenomes_nodes = kaiju_db_original_dir + '/nodes.dmp'
kaiju_db_original_progenomes_names = kaiju_db_original_dir + '/names.dmp'

map2b_checkout_dir = join(output_directory, 'MAP2B-checkout')
map2b_db = join(map2b_checkout_dir, 'database/GTDB')

## metabuli has no current available R207 database, so cannot be used generally for now. See https://github.com/steineggerlab/Metabuli/issues/63
metabuli_db_dir = '/work/microbiome/db/metabuli/gtdb207'
# metabuli_db_dir = join(output_directory, 'metabuli') #'metabuli-gtdb207'

sylph_db = join(output_directory, "gtdb_database.syldb")

# coverm: map reads against all GTDB r207 genomes with strobealign and report
# per-genome coverage. The reference is a single concatenated multi-FASTA whose
# contigs are named `<genome>~<contig>`, so `coverm genome --separator '~'`
# groups contigs back into genomes. (Staged from the historical minimap2 coverm
# db dir; we use the .fna, not the .mmi, because strobealign builds its own
# index.)
coverm_reference_fna = join(output_directory, "gtdb_r207_genomes.fna")

metakssd_checkout_dir = join(output_directory, 'MetaKSSD-checkout')
metakssd_markerdb = join(output_directory, "GTDBr207_genomes_L3K11_sketch_markerdb")

# singlem-regime3: the sylph-condense-regime3 SingleM branch (checked out as the
# singlem_sylph_condense_regime submodule, env `singlem-regime3`). It reuses the
# standard r207 SingleM metapackage, but its --joint condense also needs a sylph
# profile, which is produced by the `weebill` sylph fork (built from the weebill
# submodule) against a GTDB r207 two-stage sylph database built at -c 100.
weebill_binary = abspath(join('..', 'weebill', 'target', 'release', 'weebill'))
weebill_db = '/scratch/microbiome/woodcrob/non_sensitive/weebill_dbs/r207.100.syl2db'

# ---------------------------------------------------------------------------
# GTDB r232-backed databases (benchmark 7). Everything above targets GTDB r207;
# these are the r232 equivalents. The metapackage was fetched fresh with
# `singlem data` (the /work copy was an unreadable HSM-released stub); it extracts
# the S6.5.0 R232 ZenodoBackpack under ../tool_reference_data/singlem_data_r232/.
# The sylph / weebill databases already live on weka scratch.
singlem_metapackage_r232 = join(
    output_directory, 'singlem_data_r232', 'S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb')
# Standalone sylph db, sketched at c=100 (sylph's -c is ignored for pre-sketched
# databases, so it profiles fine against the default c=200 query sketch).
sylph_db_r232 = '/scratch/microbiome/woodcrob/non_sensitive/weebill_dbs/r232.100.syldb'
# weebill (sylph fork) two-stage r232 db, the -c 100 analogue of weebill_db above.
weebill_db_r232 = '/scratch/microbiome/woodcrob/non_sensitive/weebill_dbs/r232.100.syl2db'
# r232 GTDB taxonomy / metadata, staged to the repo root next to the r207 copies.
gtdb_bac_taxonomy_r232 = '../bac120_taxonomy_r232.tsv'
gtdb_ar_taxonomy_r232 = '../ar53_taxonomy_r232.tsv'
gtdb_bac_metadata_r232 = '../bac120_metadata_r232.tsv'
gtdb_ar_metadata_r232 = '../ar53_metadata_r232.tsv'

tools = ['singlem', 'metaphlan', 'motus', 'kraken', 'sourmash', 'kaiju', 'map2b', 'metabuli', 'sylph', 'metaphlan42', 'metakssd', 'singlem-regime3']

tools_with_filled_output_profiles = ('kraken','sourmash')

r207_tools = ['singlem', 'metaphlan', 'kraken', 'sourmash', 'metabuli', 'sylph', 'metaphlan']

gtdb_bac_metadata = '../bac120_metadata_r207.tsv'
gtdb_ar_metadata = '../ar53_metadata_r207.tsv'
