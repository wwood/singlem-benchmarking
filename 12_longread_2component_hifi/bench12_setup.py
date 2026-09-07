import os
from os.path import join

from tool_reference_data import *

# Three-tool comparison on the same two-member community as benchmark 11, sequenced as
# PacBio HiFi instead of Nanopore. metaphlan is absent for the same reason as there:
# the pinned 4.0.6 has no long-read mode (--long_reads arrived in 4.2).
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

# Single synthetic community: two known GTDB r207 species at 10x coverage each --
# the same genomes and coverages as benchmark 11, so the only difference between the
# two benchmarks is read accuracy.
datasets = ['hifi2']

##################################################################### reference databases
# Local (per-run) copies of the tool databases, staged by the copy rules.

singlem_metapackage_local = join(output_dirs_dict['singlem'], 'data', os.path.basename(singlem_metapackage))

sylph_db_local = output_dirs_dict['sylph'] + '/sylph/data/gtdb_database.syldb'

# singlem-regime3 reuses the standard r207 metapackage plus the weebill sylph db.
singlem_regime3_metapackage_local = join(
    output_dirs_dict['singlem-regime3'], 'data', os.path.basename(singlem_metapackage))
