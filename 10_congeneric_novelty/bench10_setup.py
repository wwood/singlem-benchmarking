import os
from os.path import join

from tool_reference_data import *

# Four-tool comparison, matching benchmarks 5/6/8/9, over a four-genome community
# holding two genera at three levels of novelty x abundance.
tools = ['singlem', 'sylph', 'singlem-regime3', 'metaphlan']

output_prefix = 'output_'
output_dirs = list([output_prefix + tool for tool in tools])
output_dirs_dict = dict(zip(tools, output_dirs))

benchmark_dir = 'benchmarks'
# No coverage_definitions_folder: the community (coverage *and* truth taxonomy) is
# stated in community.tsv and simulated by ../rules/data_generation_community.smk,
# not by ../rules/data_generation.smk. See the Snakefile.
# fastq_dir and truth_dir are set in the Snakefile; see the note there.

# Single synthetic community: two genera x (one known species at 50x + one
# r214-novel congener, at 10x in one genus and 100x in the other).
datasets = ['congeneric4']

##################################################################### reference databases
# Local (per-run) copies of the tool databases, staged by ../rules/staging.smk.

singlem_metapackage_local = join(output_dirs_dict['singlem'], 'data', os.path.basename(singlem_metapackage))

metaphlan_db_local1 = output_dirs_dict['metaphlan'] + '/metaphlan/data/metaphlan_bowtiedb'

sylph_db_local = output_dirs_dict['sylph'] + '/sylph/data/gtdb_database.syldb'

# singlem-regime3 reuses the standard r207 metapackage plus the weebill sylph db.
singlem_regime3_metapackage_local = join(
    output_dirs_dict['singlem-regime3'], 'data', os.path.basename(singlem_metapackage))
