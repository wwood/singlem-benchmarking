# metabuli profiling: reads -> report -> condensed profile.
# DB staging is in rules/staging.smk (rule metabuli_copy_db).
#
# Required names (includer): num_threads, fastq_dir, benchmark_dir,
#   output_dirs_dict, metabuli_db_local.
# Optional (includer may override before include):
#   metabuli_threads   threads for metabuli classify (default: num_threads).
#   metabuli_mem_mb    RAM request (default 256000; metabuli failed at 128GB).
# Reads GTDB r207 taxonomy at ../{bac120,ar53}_taxonomy_r207.tsv.

metabuli_threads = globals().get("metabuli_threads", num_threads)
metabuli_mem_mb = globals().get("metabuli_mem_mb", 256000)

rule metabuli_run:
    input:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        db=metabuli_db_local,
        done=output_dirs_dict['metabuli'] + "/metabuli/data/done",
    output:
        report=output_dirs_dict['metabuli'] + "/output/{sample}_report.tsv",
        done=touch(output_dirs_dict['metabuli'] + "/output/{sample}.done"),
    threads: metabuli_threads
    resources:
        mem_mb=metabuli_mem_mb,
    benchmark:
        benchmark_dir + "/metabuli/{sample}-" + str(num_threads) + "threads.benchmark"
    params:
        output_dir=output_dirs_dict['metabuli'],
    shell:
        'eval "$(pixi shell-hook -e metabuli)" && '
        "metabuli classify --threads {threads} {input.r1} {input.r2} {input.db} "
        "{params.output_dir}/output {wildcards.sample}"

rule metabuli_report_to_condensed:
    input:
        report=output_dirs_dict['metabuli'] + "/output/{sample}_report.tsv",
    output:
        profile=output_dirs_dict['metabuli'] + "/metabuli/{sample}.profile",
    resources:
        runtime=60,  # trivial conversion
    shell:
        'eval "$(pixi shell-hook -e singlem)" && '
        "{workflow.basedir}/../bin/metabuli_to_condensed.py --input {input} "
        "--bacterial-taxonomy ../bac120_taxonomy_r207.tsv "
        "--archaeal-taxonomy ../ar53_taxonomy_r207.tsv > {output.profile}"
