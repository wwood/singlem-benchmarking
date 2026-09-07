import os
from os.path import join

from tool_reference_data import *

# Three-tool comparison over benchmark 6's 1000 species, none above 1x coverage,
# sequenced as Oxford Nanopore R10.4.1 long reads. metaphlan (which benchmark 6 runs)
# is absent: the pinned 4.0.6 has no long-read mode (--long_reads arrived in 4.2), so
# comparing it here would measure the missing feature rather than the tool. The
# singlem / sylph / singlem-regime3 trio is the subset of benchmark 6's tools that
# accepts long reads, and is also exactly benchmark 7's trio.
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

# Benchmark 6's community, reused byte-for-byte (genome_list.tsv and the coverage
# file are copies of benchmark 6's, verified identical): 1000 distinct known GTDB
# r207 species, every one below 1x coverage and below 0.4% relative abundance. Only
# the read technology differs, which is what makes 6-vs-13 a read-technology
# comparison rather than a new community design.
datasets = ['nanopore1000']

##################################################################### reference databases
# Local (per-run) copies of the tool databases, staged by the copy rules. Identical
# to benchmark 6's: the long-read rules change only how reads are passed.

singlem_metapackage_local = join(output_dirs_dict['singlem'], 'data', os.path.basename(singlem_metapackage))

sylph_db_local = output_dirs_dict['sylph'] + '/sylph/data/gtdb_database.syldb'

# singlem-regime3 reuses the standard r207 metapackage plus the weebill sylph db.
singlem_regime3_metapackage_local = join(
    output_dirs_dict['singlem-regime3'], 'data', os.path.basename(singlem_metapackage))
