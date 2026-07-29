import os
from os.path import join

from tool_reference_data import *

# Four-tool comparison at very low (0.2x) coverage.
tools = ['singlem', 'sylph', 'singlem-regime3', 'metaphlan']

output_prefix = 'output_'
output_dirs = list([output_prefix + tool for tool in tools])
output_dirs_dict = dict(zip(tools, output_dirs))

benchmark_dir = 'benchmarks'
coverage_definitions_folder = 'coverage_definitions'
# fastq_dir and truth_dir are set (absolute) in the Snakefile; see the note there.

# Single synthetic community: 5 known GTDB r207 species, each at 0.2x coverage.
datasets = ['known5']

##################################################################### reference databases
# Local (per-run) copies of the tool databases, staged by the copy rules.

singlem_metapackage_local = join(output_dirs_dict['singlem'], 'data', os.path.basename(singlem_metapackage))

metaphlan_db_local1 = output_dirs_dict['metaphlan'] + '/metaphlan/data/metaphlan_bowtiedb'

sylph_db_local = output_dirs_dict['sylph'] + '/sylph/data/gtdb_database.syldb'

# singlem-regime3 reuses the standard r207 metapackage plus the weebill sylph db.
singlem_regime3_metapackage_local = join(
    output_dirs_dict['singlem-regime3'], 'data', os.path.basename(singlem_metapackage))
