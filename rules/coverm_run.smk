# coverm profiling: map reads against the concatenated GTDB r207 genomes FASTA
# (contigs named <genome>~<contig>) with strobealign, report per-genome trimmed
# mean coverage, then convert to a condensed profile. The strobealign index is
# built once in rules/staging.smk (rule coverm_strobealign_index).
#
# Required names (includer): num_threads, fastq_dir, benchmark_dir,
#   output_dirs_dict, coverm_reference_fna.
# Reads GTDB r207 taxonomy at ../{bac120,ar53}_taxonomy_r207.tsv.

rule coverm_run:
    input:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        ref=coverm_reference_fna,
        index=coverm_reference_fna + ".sti",
    output:
        tsv=output_dirs_dict['coverm'] + "/coverm/{sample}.tsv",
    threads: num_threads
    resources:
        mem_mb=500000,
        runtime=180,  # strobealign map vs large ref; 3h ceiling
    benchmark:
        benchmark_dir + "/coverm/{sample}-" + str(num_threads) + "threads.benchmark"
    log:
        output_dirs_dict['coverm'] + "/logs/coverm/{sample}.log",
    shell:
        'eval "$(pixi shell-hook -e coverm)" && '
        "coverm genome -1 {input.r1} -2 {input.r2} "
        "--reference {input.ref} --mapper strobealign --strobealign-use-index "
        "--separator '~' -m trimmed_mean -t {threads} > {output.tsv} 2> {log}"

rule coverm_to_condensed:
    input:
        output_dirs_dict['coverm'] + "/coverm/{sample}.tsv",
    output:
        profile=output_dirs_dict['coverm'] + "/coverm/{sample}.profile",
    resources:
        runtime=60,  # trivial conversion
    shell:
        'eval "$(pixi shell-hook -e art)" && '  # for polars
        "{workflow.basedir}/../bin/coverm_to_condensed.py "
        "--coverm-genome {input} --sample {wildcards.sample} "
        "--bac-tax ../bac120_taxonomy_r207.tsv "
        "--arc-tax ../ar53_taxonomy_r207.tsv > {output.profile}"
