import os
from os.path import join, basename

from tool_reference_data import *

# Three-tool comparison, identical community to benchmark 6 (1000 known GTDB
# species, none above 1x coverage) but every database backed by GTDB **r232**
# instead of r207. Only the tools that have an r232 database are run: singlem,
# the singlem-regime3 joint method, and sylph.
tools = ['singlem', 'sylph', 'singlem-regime3']

output_prefix = 'output_'
output_dirs = list([output_prefix + tool for tool in tools])
output_dirs_dict = dict(zip(tools, output_dirs))

benchmark_dir = 'benchmarks'
coverage_definitions_folder = 'coverage_definitions'
# fastq_dir and truth_dir are set in the Snakefile; see the note there.

# Single synthetic community, reused verbatim from benchmark 6. 119 of its 1000
# genomes were dropped from GTDB between r207 and r232 and so are absent from the
# r232 metadata; generate_community.py inner-joins on the metadata, so the r232
# community is the 881 surviving genomes (reads and truth stay consistent). See
# the README.
datasets = ['known1000']

##################################################################### reference databases
# Point the *source* database paths (imported from tool_reference_data) at the
# r232 builds. The staging rules and *_local paths below then stage/consume these.
singlem_metapackage = singlem_metapackage_r232
sylph_db = sylph_db_r232
weebill_db = weebill_db_r232

# The ground truth and the sylph->condensed mapping must use the r232 taxonomy so
# that tool output (r232 lineages) is scored against an r232 truth. These names are
# read by rules/data_generation.smk (via the Snakefile) and rules/sylph_run.smk.
gtdb_bac_metadata = gtdb_bac_metadata_r232
gtdb_ar_metadata = gtdb_ar_metadata_r232
sylph_bac_taxonomy = gtdb_bac_taxonomy_r232
sylph_ar_taxonomy = gtdb_ar_taxonomy_r232

# The r232 sylph databases are far larger than the r207 ones, and sylph/weebill
# load the whole db into RAM. Raise the per-rule memory requests accordingly
# (read by rules/sylph_run.smk and rules/singlem_regime3_run.smk); each rule is a
# separate aqua/PBS job via the snakemake aqua profile, so these become -m to mqsub.
sylph_mem_mb = 100000            # standalone sylph db is 52GB
regime3_weebill_mem_mb = 80000   # weebill two-stage db is 38GB

# The r232 sylph db is sketched at c=100; `sylph profile` must sketch the reads at
# a -c no larger than that, but its default is 200 (which makes sylph refuse the
# sample). Match the db. (The r207 db was c=200 = default, so no override there.)
sylph_c = 100

# Local (per-run) copies of the tool databases, staged by the copy rules.
singlem_metapackage_local = join(output_dirs_dict['singlem'], 'data', basename(singlem_metapackage))

sylph_db_local = output_dirs_dict['sylph'] + '/sylph/data/' + basename(sylph_db)

# singlem-regime3 reuses the standard r232 metapackage plus the weebill r232 db.
singlem_regime3_metapackage_local = join(
    output_dirs_dict['singlem-regime3'], 'data', basename(singlem_metapackage))
