# Shared synthetic-community read generation via ART (art_illumina).
#
# Wraps 1_novel_strains/generate_community.py: given a per-sample coverage file
# and a genome list, it simulates 150bp paired reads at the requested per-genome
# coverage and writes the ground-truth condensed profile.
#
# This is the generic form of 1_novel_strains' generate_community_and_reads: the
# coverage file for a sample is coverage_definitions_folder/<sample>.tsv (rather
# than the marine<N> -> coverage<N>.tsv mapping bench1 hardcodes).
#
# Required names (defined by the includer):
#   fastq_dir, truth_dir, num_threads
#   generate_community_script   path to generate_community.py
#   coverage_definitions_folder folder of <sample>.tsv coverage files
#   genome_fasta_paths          TSV of <genome>\t<fasta path>
#   gtdb_bac_metadata, gtdb_ar_metadata  (from tool_reference_data)
#
# generate_community.py reads ../bac120_metadata_r207.tsv relative to CWD, so the
# includer's shell.prefix must cd into the benchmark folder (one level under the
# repo root).


rule generate_community_and_reads:
    output:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        condensed=truth_dir + "/{sample}.condensed",
        genomewise=truth_dir + "/{sample}.genomewise.csv",
        done=touch(truth_dir + "/{sample}.finished"),
        done2=touch(fastq_dir + "/{sample}.finished"),
    input:
        coverage=coverage_definitions_folder + "/{sample}.tsv",
        genome_list=genome_fasta_paths,
    params:
        # generate_community.py chdirs into a tmpdir mid-run and abspaths all of
        # its outputs EXCEPT --output-genomewise-coverage, so pass that one as an
        # absolute path (resolved against snakemake's working dir).
        genomewise_abs=lambda wildcards, output: os.path.abspath(output.genomewise),
    threads: num_threads
    resources:
        runtime=120,  # ART simulation; small communities finish in minutes
    log:
        truth_dir + "/{sample}.readgen.log",
    shell:
        'eval "$(pixi shell-hook -e art)" && '
        "mkdir -p {fastq_dir} {truth_dir} && "
        "{generate_community_script} --art art_illumina --threads {threads} "
        "--coverage-file {input.coverage} "
        "--gtdb-bac-metadata {gtdb_bac_metadata} --gtdb-ar-metadata {gtdb_ar_metadata} "
        "--genome-list {input.genome_list} "
        "--output-condensed {output.condensed} "
        "--output-genomewise-coverage {params.genomewise_abs} "
        "-1 {output.r1} -2 {output.r2} &> {log}"
