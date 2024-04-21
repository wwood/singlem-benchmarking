import os

from tool_reference_data import *

output_prefix = 'output_'
output_dirs = list([output_prefix+tool for tool in tools])
output_dirs_dict = dict(zip(tools, output_dirs))

benchmark_dir = 'benchmarks'

datasets = list(['marine'+str(i) for i in range(10)])

fastq_dir = 'local_reads'
coverage_definitions_folder = 'coverage_definitions'

# Restrict to r207-consistent tools
tools = r207_tools

# debug
# tools = ['singlem']
# datasets = [datasets[6]] # debug

##################################################################### reference databases

# metapackage = '/work/microbiome/msingle/mess/163_S3.1.1_metapackage/output_dir/metapackage/S3.1.1.GTDB_r207.metapackage_20240302.smpkg'
singlem_metapackage_local = join(output_dirs_dict['singlem'], 'data', os.path.basename(singlem_metapackage))

# metaphlan_db_original1 = '/work/microbiome/msingle/mess/115_camisim_ish_benchmarking/metaphlan_bowtiedb'
metaphlan_db_local1 = output_dirs_dict['metaphlan'] + '/metaphlan/data/metaphlan_bowtiedb'

# motus_db_path_original = '/work/microbiome/msingle/mess/124_singlem-benchmarking/db_mOTU'
# motus_db_path_local = output_dirs_dict['motus'] + '/motus/data/db_mOTU'
# motus_gtdb_tsv = '/work/microbiome/msingle/mess/115_camisim_ish_benchmarking/motus/mOTUs_3.0.0_GTDB_tax.tsv'

# kraken_db = "/work/microbiome/db/struo2/GTDB_release207"
kraken_db_local = output_dirs_dict['kraken'] + "/kraken/data/GTDB_release207"

# sourmash_db1_original = '/work/microbiome/msingle/mess/105_novelty_testing/gtdb-rs207.taxonomy.sqldb'
# sourmash_db2_original = '/work/microbiome/msingle/mess/105_novelty_testing/gtdb-rs207.genomic-reps.dna.k31.zip'
sourmash_db_taxonomy_local = output_dirs_dict['sourmash'] + '/sourmash/data/gtdb-rs207.taxonomy.sqldb'
sourmash_db_dna_local = output_dirs_dict['sourmash'] + '/sourmash/data/gtdb-rs207.genomic-reps.dna.k31.zip'

# kaiju_db_original_dir = '/work/microbiome/msingle/mess/128_kaiju/progenomes'
# kaiju_db_original_progenomes_fmi = kaiju_db_original_dir + '/kaiju_db_progenomes.fmi'
# kaiju_db_original_progenomes_nodes = kaiju_db_original_dir + '/nodes.dmp'
# kaiju_db_original_progenomes_names = kaiju_db_original_dir + '/names.dmp'
kaiju_db = output_dirs_dict['kaiju'] + '/kaiju/data/kaiju_db_progenomes.fmi'
kaiju_db_progenomes_nodes = output_dirs_dict['kaiju'] + '/kaiju/data/nodes.dmp'
kaiju_db_progenomes_names = output_dirs_dict['kaiju'] + '/kaiju/data/names.dmp'

# map2b_checkout_dir = '/home/woodcrob/bioinfo/MAP2B'
# map2b_db = os.path.join(map2b_checkout_dir, 'database.bak/GTDB')
# map2b_db_local = output_dirs_dict['map2b'] + '/map2b/data/GTDB'

# metabuli_db = '/work/microbiome/db/metabuli/gtdb207'
metabuli_db_local = output_dirs_dict['metabuli'] + '/metabuli/data/gtdb207'