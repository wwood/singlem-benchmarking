
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

metakssd_checkout_dir = join(output_directory, 'MetaKSSD-checkout')
metakssd_markerdb = join(output_directory, "GTDBr207_genomes_L3K11_sketch_markerdb")

tools = ['singlem', 'metaphlan', 'motus', 'kraken', 'sourmash', 'kaiju', 'map2b', 'metabuli', 'sylph', 'metaphlan42', 'metakssd']

tools_with_filled_output_profiles = ('kraken','sourmash')

r207_tools = ['singlem', 'metaphlan', 'kraken', 'sourmash', 'metabuli', 'sylph', 'metaphlan']

gtdb_bac_metadata = '../bac120_metadata_r207.tsv'
gtdb_ar_metadata = '../ar53_metadata_r207.tsv'
