# sylph profiling (standard bioconda, single-stage `sylph profile`):
# reads -> sylph genome TSV -> condensed profile.
# Database staging is in rules/staging.smk (rule sylph_copy_db).
#
# Required names (includer): num_threads, fastq_dir, benchmark_dir,
#   output_dirs_dict, sylph_db_local.
# Reads the GTDB r207 taxonomy at ../{bac120,ar53}_taxonomy_r207.tsv (relative to
# the benchmark dir, so the includer's shell.prefix must cd there).

rule sylph_run:
    input:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        db=sylph_db_local,
        done=output_dirs_dict['sylph'] + "/sylph/data/done",
    output:
        report=output_dirs_dict['sylph'] + "/output/{sample}.tsv",
        done=touch(output_dirs_dict['sylph'] + "/output/{sample}.done"),
    threads: num_threads
    resources:
        mem_mb=32000,
        runtime=120,  # sylph is fast; 2h is ample
    benchmark:
        benchmark_dir + "/sylph/{sample}-" + str(num_threads) + "threads.benchmark"
    log:
        output_dirs_dict['sylph'] + "/logs/sylph/{sample}.log",
    shell:
        "pixi run --environment sylph sylph profile {input.db} "
        "-1 {input.r1} -2 {input.r2} -t {threads} > {output.report} 2> {log}"

rule sylph_report_to_condensed:
    input:
        report=output_dirs_dict['sylph'] + "/output/{sample}.tsv",
    output:
        profile=output_dirs_dict['sylph'] + "/sylph/{sample}.profile",
    resources:
        runtime=60,  # trivial conversion
    shell:
        'eval "$(pixi shell-hook -e art)" && '  # for polars
        "python3 {workflow.basedir}/../bin/sylph_to_condensed.py --sylph-genome {input} "
        "--sample {wildcards.sample} "
        "--bac-tax ../bac120_taxonomy_r207.tsv "
        "--arc-tax ../ar53_taxonomy_r207.tsv > {output.profile}"
