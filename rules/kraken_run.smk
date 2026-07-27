# kraken2 + bracken profiling: reads -> kraken report -> per-rank bracken
# reports -> condensed profile. DB staging is in rules/staging.smk
# (rule braken_copy_data).
#
# Required names (includer): num_threads, fastq_dir, benchmark_dir, truth_dir,
#   output_dirs_dict, kraken_db_local.
# Optional (includer may override before include):
#   kraken_num_threads   threads for kraken2/kaiju (default: 1; kraken exhausts
#                        RAM if over-parallelised).
#   kraken_mem_mb        RAM request for kraken_run (default 512000).
# Reads GTDB r207 taxonomy at ../{bac120,ar53}_taxonomy_r207.tsv.

kraken_num_threads = globals().get("kraken_num_threads", 1)
kraken_mem_mb = globals().get("kraken_mem_mb", 512000)

rule kraken_run:
    input:
        reads_copied1=fastq_dir + "/{sample}.1.fq.gz",
        reads_copied2=fastq_dir + "/{sample}.2.fq.gz",
        db=kraken_db_local,
        copy_braken_data_done=output_dirs_dict['kraken'] + "/kraken/data/done",
    benchmark:
        benchmark_dir + "/kraken/{sample}-" + str(num_threads) + "threads.benchmark"
    threads: kraken_num_threads
    resources:
        mem_mb=kraken_mem_mb,
    output:
        report=output_dirs_dict['kraken'] + "/kraken/{sample}.kraken",
        done=touch(output_dirs_dict['kraken'] + "/kraken/{sample}.kraken.done"),
    log:
        output_dirs_dict['kraken'] + "/logs/kraken/{sample}.log",
    shell:
        'eval "$(pixi shell-hook -e kraken)" && '
        "kraken2 --db {input.db} --threads {num_threads} --output /dev/null "
        "--report {output.report} --paired {input.reads_copied1} {input.reads_copied2} &> {log}"

rule braken_run:
    input:
        db=kraken_db_local,
        kraken_report=output_dirs_dict['kraken'] + "/kraken/{sample}.kraken",
        kraken_done=output_dirs_dict['kraken'] + "/kraken/{sample}.kraken.done",
    output:
        s=output_dirs_dict['kraken'] + "/braken/{sample}.report.S",
        g=output_dirs_dict['kraken'] + "/braken/{sample}.report.G",
        f=output_dirs_dict['kraken'] + "/braken/{sample}.report.F",
        o=output_dirs_dict['kraken'] + "/braken/{sample}.report.O",
        c=output_dirs_dict['kraken'] + "/braken/{sample}.report.C",
        p=output_dirs_dict['kraken'] + "/braken/{sample}.report.P",
        d=output_dirs_dict['kraken'] + "/braken/{sample}.report.D",
        done=touch(output_dirs_dict['kraken'] + "/braken/{sample}.done"),
    log:
        output_dirs_dict['kraken'] + "/logs/bracken/{sample}.log",
    shell:
        'eval "$(pixi shell-hook -e bracken)" && '
        "bash -c 'bracken -d {input.db} -r 150 -l S -t 10 -o {output.s} -i {input.kraken_report} && "
        "bracken -d {input.db} -r 150 -l G -t 10 -o {output.g} -i {input.kraken_report} && "
        "bracken -d {input.db} -r 150 -l F -t 10 -o {output.f} -i {input.kraken_report} && "
        "bracken -d {input.db} -r 150 -l O -t 10 -o {output.o} -i {input.kraken_report} && "
        "bracken -d {input.db} -r 150 -l C -t 10 -o {output.c} -i {input.kraken_report} && "
        "bracken -d {input.db} -r 150 -l P -t 10 -o {output.p} -i {input.kraken_report} && "
        "bracken -d {input.db} -r 150 -l D -t 10 -o {output.d} -i {input.kraken_report}' &> {log}"

rule bracken_to_profile:
    input:
        output_dirs_dict['kraken'] + "/braken/{sample}.report.S",
        output_dirs_dict['kraken'] + "/braken/{sample}.report.G",
        output_dirs_dict['kraken'] + "/braken/{sample}.report.F",
        output_dirs_dict['kraken'] + "/braken/{sample}.report.O",
        output_dirs_dict['kraken'] + "/braken/{sample}.report.C",
        output_dirs_dict['kraken'] + "/braken/{sample}.report.P",
        output_dirs_dict['kraken'] + "/braken/{sample}.report.D",
        output_dirs_dict['kraken'] + "/braken/{sample}.done",
    params:
        report_prefix=output_dirs_dict['kraken'] + "/braken/{sample}.report",
        biobox_dir=output_dirs_dict['kraken'] + "/biobox",
    output:
        profile=output_dirs_dict['kraken'] + "/kraken/{sample}.profile",
    log:
        output_dirs_dict['kraken'] + "/logs/bracken_to_profile/{sample}.log",
    shell:
        # Convert reports to singlem condense format
        'eval "$(pixi shell-hook -e singlem)" && '
        "mkdir -p {params.biobox_dir} && "
        "{workflow.basedir}/../bin/kraken_to_condensed.py --report-prefix {params.report_prefix} "
        "--bacterial-taxonomy ../bac120_taxonomy_r207.tsv "
        "--archaeal-taxonomy ../ar53_taxonomy_r207.tsv > {output.profile} 2> {log}"
