import os
from os.path import join

from tool_reference_data import *

# Three-tool comparison over benchmark 2's novel-lineage pairs, sequenced as Oxford
# Nanopore R10.4.1 long reads.
#
# Benchmark 2 runs eleven tools split into relative-abundance and reads-wise
# (read-count) families, scored against two different truth files via a {method}
# wildcard. Here only the three long-read-capable tools remain, and all three report
# relative abundance -- so there is exactly one truth per sample and no {method} layer,
# which is what lets this benchmark use the shared ../rules/common.smk scoring rules
# instead of benchmark 2's bespoke {method} copies. (kraken / kaiju / metabuli are the
# reads-wise tools; if one is ever added here, the reads-wise truth half must be
# restored -- see bin/generate_longread_2component_community.py.)
tools = ['singlem', 'sylph', 'singlem-regime3']

output_prefix = 'output_'
output_dirs = list([output_prefix + tool for tool in tools])
output_dirs_dict = dict(zip(tools, output_dirs))

benchmark_dir = 'benchmarks'
# fastq_dir and truth_dir are set in the Snakefile; see the note there.

# Badread parameters for Oxford Nanopore R10.4.1; consumed by the read-generation rule.
longread_read_tech = 'nanopore'

# The community definition is benchmark 2's, read from its committed tables (see the
# Snakefile): 120 samples, each one novel genome plus one known species at 10x each,
# spanning six novelty depths (species/genus/family/order/class/phylum, 20 each).
bench2_dir = '../2_phylogenetic_novelty'

##################################################################### reference databases
# Local (per-run) copies of the tool databases, staged by the copy rules.

singlem_metapackage_local = join(output_dirs_dict['singlem'], 'data', os.path.basename(singlem_metapackage))

sylph_db_local = output_dirs_dict['sylph'] + '/sylph/data/gtdb_database.syldb'

# singlem-regime3 reuses the standard r207 metapackage plus the weebill sylph db.
singlem_regime3_metapackage_local = join(
    output_dirs_dict['singlem-regime3'], 'data', os.path.basename(singlem_metapackage))
