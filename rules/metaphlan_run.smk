# MetaPhlAn profiling: reads -> SGB report -> GTDB profile -> condensed profile.
# Database staging is in rules/staging.smk (rule metaphlan_copy_db).
#
# Required names (includer): num_threads, fastq_dir, benchmark_dir,
#   output_dirs_dict, metaphlan_db_local1, metaphlan_index.
# Reads at {fastq_dir}/{sample}.{1,2}.fq.gz.

rule metaphlan_profile:
    input:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        done=output_dirs_dict['metaphlan'] + "/metaphlan/data/done",
    benchmark:
        benchmark_dir + "/metaphlan/{sample}-" + str(num_threads) + "threads.benchmark"
    output:
        sgb_report=output_dirs_dict['metaphlan'] + "/metaphlan/{sample}.sgb_report",
        done=touch(output_dirs_dict['metaphlan'] + "/metaphlan/{sample}.profile.done"),
    threads: num_threads
    resources:
        runtime=180,  # bowtie2 vs SGB db; 3h ceiling
    params:
        cat_reads=output_dirs_dict['metaphlan'] + "/metaphlan/{sample}.cat.fq.gz",
    log:
        output_dirs_dict['metaphlan'] + "/logs/metaphlan/{sample}.log",
    shell:
        # metaphlan can't take multiple input files, so concatenate the mates first.
        'eval "$(pixi shell-hook -e metaphlan)" && '
        "rm -f {output.sgb_report} {params.cat_reads}.bowtie2out.txt; "
        "cat {input.r1} {input.r2} > {params.cat_reads} && "
        "metaphlan {params.cat_reads} --index {metaphlan_index} --nproc {threads} "
        "--input_type fastq --bowtie2db {metaphlan_db_local1} -o {output.sgb_report} &> {log}"

rule metaphlan_convert_profile_to_GTDB:
    input:
        report=output_dirs_dict['metaphlan'] + "/metaphlan/{sample}.sgb_report",
    output:
        gtdb_report=output_dirs_dict['metaphlan'] + "/metaphlan/{sample}.gtdb_profile",
        done=touch(output_dirs_dict['metaphlan'] + "/metaphlan/{sample}.gtdb_report.done"),
    resources:
        runtime=60,  # trivial conversion
    log:
        output_dirs_dict['metaphlan'] + "/logs/metaphlan/{sample}-convert.log",
    shell:
        'eval "$(pixi shell-hook -e metaphlan)" && '
        "sgb_to_gtdb_profile.py -i {input.report} -o {output.gtdb_report} "
        "-d {metaphlan_db_local1}/mpa_vOct22_CHOCOPhlAnSGB_202212.pkl &> {log}"

rule metaphlan_profile_to_condensed:
    input:
        report=output_dirs_dict['metaphlan'] + "/metaphlan/{sample}.gtdb_profile",
    output:
        profile=output_dirs_dict['metaphlan'] + "/metaphlan/{sample}.profile",
    resources:
        runtime=60,  # trivial conversion
    shell:
        'eval "$(pixi shell-hook -e singlem)" && '
        "{workflow.basedir}/../bin/metaphlan_to_condensed.py "
        "--metaphlan {input} --sample {wildcards.sample} > {output.profile}"
