# Shared synthetic-community read generation from a fully-specified community file.
#
# The sibling of rules/data_generation.smk. That rule takes coverage from
# coverage_definitions/<sample>.tsv and looks the truth taxonomy up in the GTDB r207
# metadata; this one takes both coverage and truth taxonomy from a committed
# community.tsv. Benchmarks whose community includes genomes that are new in GTDB
# r214 must use this rule: those genomes have no r207 metadata row, so the shared
# rule's inner join against the metadata would silently drop exactly the members
# under test. See bin/generate_community_from_tsv.py.
#
# Required names (defined by the includer):
#   fastq_dir, truth_dir, num_threads
#   community_definition        TSV of genome/fasta/coverage/role/taxonomy
# Optional:
#   community_runtime           queue minutes for the rule (default 120)
#
# Paths inside community.tsv are resolved relative to CWD, so the includer's
# shell.prefix must cd into the benchmark folder (one level under the repo root),
# exactly as for the other generation rules.

community_runtime = globals().get("community_runtime", 120)
generate_community_from_tsv_script = os.path.join(
    workflow.basedir, '../bin/generate_community_from_tsv.py')


rule generate_community_and_reads:
    input:
        community=community_definition,
        script=generate_community_from_tsv_script,
    output:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        condensed=truth_dir + "/{sample}.condensed",
        genomewise=truth_dir + "/{sample}.genomewise.csv",
        done=touch(truth_dir + "/{sample}.finished"),
        done2=touch(fastq_dir + "/{sample}.finished"),
    threads: num_threads
    resources:
        runtime=community_runtime,
    log:
        truth_dir + "/{sample}.readgen.log",
    shell:
        'eval "$(pixi shell-hook -e art)" && '
        "mkdir -p {fastq_dir} {truth_dir} && "
        "{input.script} --art art_illumina --threads {threads} "
        "--community {input.community} --sample {wildcards.sample} "
        "--output-condensed {output.condensed} --output-genomewise {output.genomewise} "
        "-1 {output.r1} -2 {output.r2} &> {log}"
