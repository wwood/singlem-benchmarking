#!/usr/bin/env python3
"""Simulate benchmark 10's community and write its ground-truth condensed profile.

Same ART invocation as ../1_novel_strains/generate_community.py (HSXt, 150bp
paired, -m 400 -s 10), but the coverage *and* the truth taxonomy come from the
community file rather than from the GTDB metadata. That is the whole reason this
benchmark cannot use the shared ../rules/data_generation.smk: two of the four
genomes are new in GTDB r214 and have no r207 metadata row, so the shared path's
inner join against the metadata would silently drop them -- simulating a
two-genome community and writing a truth that omits exactly the members under
test.

Coverages are also taken as written, not shuffled positionally against a
randomised metadata merge as the shared generator does. Here which genome gets
which coverage *is* the experiment (novel@10x vs novel@100x), so the pairing in
community.tsv has to be the pairing simulated.

    cd 10_congeneric_novelty && pixi run -e art python3 generate_community.py \
        --community community.tsv --sample congeneric4 \
        --output-condensed truths/congeneric4.condensed \
        --output-genomewise truths/congeneric4.genomewise.csv \
        -1 reads/congeneric4.1.fq.gz -2 reads/congeneric4.2.fq.gz --art art_illumina
"""

import argparse
import logging
import os
import sys
import tempfile
from collections import OrderedDict

import extern
import polars as pl

READ_LENGTH = 150

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--community', required=True,
                    help='TSV of genome/fasta/coverage/role/taxonomy (comment lines start with #)')
parser.add_argument('--sample', required=True, help='sample name written into the condensed profile')
parser.add_argument('--output-condensed', required=True, help='ground truth, singlem condensed format')
parser.add_argument('--output-genomewise', required=True, help='per-genome coverage/role/taxonomy record')
parser.add_argument('-1', '--read1', required=True, help='output fq.gz')
parser.add_argument('-2', '--read2', required=True, help='output fq.gz')
parser.add_argument('--art', required=True, help='path to the ART binary (art_illumina)')
parser.add_argument('--threads', type=int, default=1)
parser.add_argument('--debug', action='store_true')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG if args.debug else logging.INFO,
                    format='%(asctime)s %(levelname)s: %(message)s')

# ART chdirs into a tmpdir below, so every path leaving this process must be
# absolute before then.
output1 = os.path.abspath(args.read1)
output2 = os.path.abspath(args.read2)
output_condensed = os.path.abspath(args.output_condensed)
output_genomewise = os.path.abspath(args.output_genomewise)

community = pl.read_csv(args.community, separator='\t', comment_prefix='#')
community = community.with_columns(
    pl.col('fasta').map_elements(os.path.abspath, return_dtype=pl.String))
logging.info("Read %d community members from %s", community.height, args.community)

for fasta in community['fasta']:
    if not os.path.exists(fasta):
        raise SystemExit("community fasta does not exist: {}".format(fasta))

# Sum coverage per truth lineage. Two genomes sharing a lineage would be summed
# into one row -- which is why the community is one novel species per genus: two
# novel congeners would collapse into a single genus-only row and the 10x and 100x
# levels would stop being separable in the truth.
tax_to_coverage = OrderedDict()
for taxonomy, coverage in zip(community['taxonomy'], community['coverage']):
    tax_to_coverage[taxonomy] = tax_to_coverage.get(taxonomy, 0) + coverage
if len(tax_to_coverage) != community.height:
    logging.warning("%d genomes collapsed into %d truth lineages",
                    community.height, len(tax_to_coverage))

for d in (output_condensed, output_genomewise, output1, output2):
    os.makedirs(os.path.dirname(d), exist_ok=True)

with open(output_genomewise, 'w') as f:
    f.write("genome\tcoverage\trole\tfasta\ttaxonomy\n")
    for row in community.iter_rows(named=True):
        f.write("{genome}\t{coverage}\t{role}\t{fasta}\t{taxonomy}\n".format(**row))

with open(output_condensed, 'w') as f:
    f.write("sample\tcoverage\ttaxonomy\n")
    for taxonomy, coverage in tax_to_coverage.items():
        f.write("{}\t{}\t{}\n".format(args.sample, coverage, taxonomy))
logging.info("Wrote truth (%d lineages) to %s", len(tax_to_coverage), output_condensed)

with tempfile.TemporaryDirectory() as tmpdir:
    os.chdir(tmpdir)
    os.makedirs('simulated_reads')
    sim_commands = [
        "{} -ss HSXt -i {} -p -l {} -f {} -m 400 -s 10 -o simulated_reads/{}. &>/dev/null".format(
            args.art, fasta, READ_LENGTH, coverage, i)
        for i, (fasta, coverage) in enumerate(zip(community['fasta'], community['coverage']))
    ]
    logging.info("Simulating %d genomes ..", len(sim_commands))
    extern.run_many(sim_commands, num_threads=args.threads, progress_stream=sys.stderr)

    # `sed 's=/= ='` strips ART's /1 /2 mate suffix, as the other generators do.
    logging.info("Concatenating and compressing ..")
    extern.run("cat simulated_reads/*1.fq |sed 's=/= =' |pigz -p {} >{}".format(args.threads, output1))
    extern.run("cat simulated_reads/*2.fq |sed 's=/= =' |pigz -p {} >{}".format(args.threads, output2))
logging.info("Done.")
