
import extern
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), '..')))
from tool_reference_data import *

truth_dir = os.path.join(workflow.basedir, 'truths')
genome_fasta_paths = os.path.join(workflow.basedir, 'shadow_genome_paths.csv')

sys.path.append(os.path.abspath(os.path.dirname(os.path.abspath(workflow.snakefile))))
from bench1_setup import *

fastq_dir = os.path.join(workflow.basedir, 'generated_reads')

#####################################################################

rule all:
    input:
        expand(truth_dir + "/{sample}.condensed.biobox", sample=datasets)

rule generate_community_and_reads:
    output:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        condensed = truth_dir + "/{sample}.condensed",
        genomewise = truth_dir + "/{sample}.genomewise.csv",
        done = touch(truth_dir + "/{sample}.finished"),
        done2 = touch(fastq_dir + "/{sample}.finished"),
    params:
        coverage_number = lambda wildcards: wildcards.sample.replace('marine', ''),
    threads: num_threads
    conda:
        'envs/art.yml'
    shell:
        "{workflow.basedir}/generate_community.py --art art_illumina --threads {threads} --coverage-file {workflow.basedir}/coverage_definitions/coverage{params.coverage_number}.tsv --gtdb-bac-metadata {gtdb_bac_metadata} --gtdb-ar-metadata {gtdb_ar_metadata} --genome-list {genome_fasta_paths} --output-condensed {output.condensed} -1 {fastq_dir}/{wildcards.sample}.1.fq.gz -2 {fastq_dir}/{wildcards.sample}.2.fq.gz --output-genomewise-coverage {output.genomewise}"

rule truth_condensed_to_biobox:
    input:
        condensed = truth_dir + "/{sample}.condensed",
    output:
        biobox = truth_dir + "/{sample}.condensed.biobox",
    conda:
        "envs/singlem.yml"
    shell:
        "{workflow.basedir}/../bin/condensed_profile_to_biobox.py --input-condensed-table {input.condensed} " \
        "--output-biobox {output.biobox}"
