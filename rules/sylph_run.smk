# sylph profiling (standard bioconda, single-stage `sylph profile`):
# reads -> sylph genome TSV -> condensed profile.
# Database staging is in rules/staging.smk (rule sylph_copy_db).
#
# Required names (includer): num_threads, fastq_dir, benchmark_dir,
#   output_dirs_dict, sylph_db_local.
# Maps sylph's genome hits to lineages via the GTDB taxonomy. Defaults to the r207
# taxonomy at ../{bac120,ar53}_taxonomy_r207.tsv (relative to the benchmark dir,
# so the includer's shell.prefix must cd there). A benchmark backed by a different
# GTDB release overrides this by defining sylph_bac_taxonomy / sylph_ar_taxonomy
# before the include (e.g. benchmark 7 points them at the r232 taxonomy).
sylph_bac_taxonomy = globals().get("sylph_bac_taxonomy", "../bac120_taxonomy_r207.tsv")
sylph_ar_taxonomy = globals().get("sylph_ar_taxonomy", "../ar53_taxonomy_r207.tsv")

# `sylph profile` loads the whole database into memory, so peak RAM scales with the
# db size. The r207 db (8.5GB) fits comfortably in 32GB; the r232 db is 52GB and
# needs ~100GB. A benchmark backed by a large db overrides sylph_mem_mb before the
# include (benchmark 7 sets 100000).
sylph_mem_mb = globals().get("sylph_mem_mb", 32000)

# sylph requires the reads sketched at a -c no larger than the database's -c.
# `sylph profile` sketches reads at the default -c 200; the r207 db was built at
# c=200 (default, so no -c needed), but the r232 db is c=100 -- passing the default
# 200 makes sylph refuse the sample and report nothing. A benchmark whose db uses a
# smaller -c overrides sylph_c before the include (benchmark 7 sets 100).
sylph_c = globals().get("sylph_c", None)
sylph_c_arg = "" if sylph_c is None else "-c {} ".format(sylph_c)

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
        mem_mb=sylph_mem_mb,
        runtime=120,  # sylph is fast; 2h is ample
    benchmark:
        benchmark_dir + "/sylph/{sample}-" + str(num_threads) + "threads.benchmark"
    log:
        output_dirs_dict['sylph'] + "/logs/sylph/{sample}.log",
    shell:
        "pixi run --environment sylph sylph profile {input.db} "
        "{sylph_c_arg}-1 {input.r1} -2 {input.r2} -t {threads} > {output.report} 2> {log}"

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
        "--bac-tax {sylph_bac_taxonomy} "
        "--arc-tax {sylph_ar_taxonomy} > {output.profile}"
