import os
from os.path import join

from tool_reference_data import *

# Four-tool comparison, matching benchmarks 5 and 6, over a 12-genome community
# that is half novel (new in GTDB r214, so absent from every r207 database) and
# half Streptomyces.
tools = ['singlem', 'sylph', 'singlem-regime3', 'metaphlan']

output_prefix = 'output_'
output_dirs = list([output_prefix + tool for tool in tools])
output_dirs_dict = dict(zip(tools, output_dirs))

benchmark_dir = 'benchmarks'
# No coverage_definitions_folder: reads and truth are provided with the dataset
# rather than simulated here, so ../rules/data_generation.smk is not included.
# fastq_dir and truth_dir are set in the Snakefile; see the note there.

# Single provided community. The reads (sample5.{1,2}.fq.gz) and ground truth
# (sample5.condensed) ship with the benchmark directory.
datasets = ['sample5']

##################################################################### reference databases
# Local (per-run) copies of the tool databases, staged by ../rules/staging.smk.

singlem_metapackage_local = join(output_dirs_dict['singlem'], 'data', os.path.basename(singlem_metapackage))

metaphlan_db_local1 = output_dirs_dict['metaphlan'] + '/metaphlan/data/metaphlan_bowtiedb'

sylph_db_local = output_dirs_dict['sylph'] + '/sylph/data/gtdb_database.syldb'

# singlem-regime3 reuses the standard r207 metapackage plus the weebill sylph db.
singlem_regime3_metapackage_local = join(
    output_dirs_dict['singlem-regime3'], 'data', os.path.basename(singlem_metapackage))
