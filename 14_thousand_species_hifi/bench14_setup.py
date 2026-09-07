import os
from os.path import join

from tool_reference_data import *

# Three-tool comparison over benchmark 6's 1000 species, none above 1x coverage,
# sequenced as PacBio HiFi. Same tool set and same reasoning as benchmark 13:
# metaphlan is excluded because the pinned 4.0.6 has no long-read mode.
tools = ['singlem', 'sylph', 'singlem-regime3']

output_prefix = 'output_'
output_dirs = list([output_prefix + tool for tool in tools])
output_dirs_dict = dict(zip(tools, output_dirs))

benchmark_dir = 'benchmarks'
coverage_definitions_folder = 'coverage_definitions'
# fastq_dir and truth_dir are set in the Snakefile; see the note there.

# Badread parameters for PacBio HiFi (pacbio2021 models, normal qscore distribution);
# consumed by ../rules/longread_data_generation.smk.
longread_read_tech = 'pacbio-hifi'

# Benchmark 6's community, reused byte-for-byte (genome_list.tsv and the coverage file
# are copies of benchmark 6's, verified identical), and the same coverages as benchmark
# 13 -- so 13-vs-14 isolates read accuracy at matched length and depth.
datasets = ['hifi1000']

##################################################################### reference databases
# Local (per-run) copies of the tool databases, staged by the copy rules.

singlem_metapackage_local = join(output_dirs_dict['singlem'], 'data', os.path.basename(singlem_metapackage))

sylph_db_local = output_dirs_dict['sylph'] + '/sylph/data/gtdb_database.syldb'

# singlem-regime3 reuses the standard r207 metapackage plus the weebill sylph db.
singlem_regime3_metapackage_local = join(
    output_dirs_dict['singlem-regime3'], 'data', os.path.basename(singlem_metapackage))
