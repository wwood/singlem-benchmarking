# sylph profiling of **long, single-end** reads: reads -> sylph genome TSV ->
# condensed profile. The long-read counterpart of rules/sylph_run.smk; the only
# change is `-r` (single-end raw reads) in place of `-1`/`-2`. sylph documents
# support for long reads and was independently benchmarked on Oxford Nanopore data;
# the profiling algorithm and database format are the same, so the r207 db and the
# shared staging rule (rules/staging.smk, sylph_copy_db) are reused unchanged.
#
# Required names (includer): num_threads, fastq_dir, benchmark_dir,
#   output_dirs_dict, sylph_db_local.
# Reads at {fastq_dir}/{sample}.fq.gz.
#
# The taxonomy / memory / -c overrides behave exactly as in the short-read rule; see
# rules/sylph_run.smk for why each exists.
sylph_bac_taxonomy = globals().get("sylph_bac_taxonomy", "../bac120_taxonomy_r207.tsv")
sylph_ar_taxonomy = globals().get("sylph_ar_taxonomy", "../ar53_taxonomy_r207.tsv")
sylph_mem_mb = globals().get("sylph_mem_mb", 32000)
sylph_c = globals().get("sylph_c", None)
sylph_c_arg = "" if sylph_c is None else "-c {} ".format(sylph_c)

rule sylph_longread_run:
    input:
        reads=fastq_dir + "/{sample}.fq.gz",
        db=sylph_db_local,
        done=output_dirs_dict['sylph'] + "/sylph/data/done",
    output:
        report=output_dirs_dict['sylph'] + "/output/{sample}.tsv",
        done=touch(output_dirs_dict['sylph'] + "/output/{sample}.done"),
    threads: num_threads
    resources:
        mem_mb=sylph_mem_mb,
        runtime=120,  # sylph is fast; 2h is ample
    benchmark:
        benchmark_dir + "/sylph/{sample}-" + str(num_threads) + "threads.benchmark"
    log:
        output_dirs_dict['sylph'] + "/logs/sylph/{sample}.log",
    shell:
        "pixi run --environment sylph sylph profile {input.db} "
        "{sylph_c_arg}-r {input.reads} -t {threads} > {output.report} 2> {log}"

rule sylph_longread_report_to_condensed:
    input:
        report=output_dirs_dict['sylph'] + "/output/{sample}.tsv",
    output:
        profile=output_dirs_dict['sylph'] + "/sylph/{sample}.profile",
    resources:
        runtime=60,  # trivial conversion
    shell:
        'eval "$(pixi shell-hook -e badread)" && '  # for polars
        "python3 {workflow.basedir}/../bin/sylph_to_condensed.py --sylph-genome {input} "
        "--sample {wildcards.sample} "
        "--bac-tax {sylph_bac_taxonomy} "
        "--arc-tax {sylph_ar_taxonomy} > {output.profile}"
