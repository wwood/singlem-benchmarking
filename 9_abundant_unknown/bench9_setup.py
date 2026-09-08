import os
from os.path import join

from tool_reference_data import *

# Four-tool comparison, matching benchmarks 5 and 6, over a 10-genome community
# that is nine-tenths known in GTDB r207 plus one abundant genome that is new in
# r214 (and so absent from every r207 database).
tools = ['singlem', 'sylph', 'singlem-regime3', 'metaphlan']

output_prefix = 'output_'
output_dirs = list([output_prefix + tool for tool in tools])
output_dirs_dict = dict(zip(tools, output_dirs))

benchmark_dir = 'benchmarks'
# No coverage_definitions_folder: coverage and truth taxonomy come from
# community.tsv via ../rules/data_generation_community.smk rather than from the
# GTDB metadata. fastq_dir and truth_dir are set in the Snakefile; see the note
# there.

# Single community, defined in full in community.tsv and simulated into
# reads/marine0.{1,2}.fq.gz with truths/marine0.condensed.
datasets = ['marine0']

##################################################################### reference databases
# Local (per-run) copies of the tool databases, staged by ../rules/staging.smk.

singlem_metapackage_local = join(output_dirs_dict['singlem'], 'data', os.path.basename(singlem_metapackage))

metaphlan_db_local1 = output_dirs_dict['metaphlan'] + '/metaphlan/data/metaphlan_bowtiedb'

sylph_db_local = output_dirs_dict['sylph'] + '/sylph/data/gtdb_database.syldb'

# singlem-regime3 reuses the standard r207 metapackage plus the weebill sylph db.
singlem_regime3_metapackage_local = join(
    output_dirs_dict['singlem-regime3'], 'data', os.path.basename(singlem_metapackage))
