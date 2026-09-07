import os
from os.path import join

from tool_reference_data import *

# Three-tool comparison on a deliberately easy two-member community, sequenced as
# Oxford Nanopore R10.4.1 long reads. metaphlan is absent: the pinned 4.0.6 has no
# long-read mode at all (--long_reads arrived in 4.2), so there is no honest way to
# run it here. The three tools listed are all verified to accept single-end long
# reads.
tools = ['singlem', 'sylph', 'singlem-regime3']

output_prefix = 'output_'
output_dirs = list([output_prefix + tool for tool in tools])
output_dirs_dict = dict(zip(tools, output_dirs))

benchmark_dir = 'benchmarks'
coverage_definitions_folder = 'coverage_definitions'
# fastq_dir and truth_dir are set in the Snakefile; see the note there.

# Badread parameters for Oxford Nanopore R10.4.1; consumed by
# ../rules/longread_data_generation.smk.
longread_read_tech = 'nanopore'

# Single synthetic community: two known GTDB r207 species at 10x coverage each.
datasets = ['nanopore2']

##################################################################### reference databases
# Local (per-run) copies of the tool databases, staged by the copy rules. All
# unchanged from the short-read benchmarks: the long-read rules differ only in how
# reads are passed, not in the databases they read.

singlem_metapackage_local = join(output_dirs_dict['singlem'], 'data', os.path.basename(singlem_metapackage))

sylph_db_local = output_dirs_dict['sylph'] + '/sylph/data/gtdb_database.syldb'

# singlem-regime3 reuses the standard r207 metapackage plus the weebill sylph db.
singlem_regime3_metapackage_local = join(
    output_dirs_dict['singlem-regime3'], 'data', os.path.basename(singlem_metapackage))
