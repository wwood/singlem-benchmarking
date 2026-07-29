# SingleM profiling: reads -> condensed taxonomic profile.
# Database staging is in rules/staging.smk (rule singlem_copy_metapackage).
#
# Required names (includer): num_threads, fastq_dir, benchmark_dir,
#   output_dirs_dict, singlem_metapackage_local.
# Reads at {fastq_dir}/{sample}.{1,2}.fq.gz.

rule singlem_run_to_profile:
    input:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        db=singlem_metapackage_local,
        data_done=output_dirs_dict['singlem'] + "/singlem/data/done",
    benchmark:
        benchmark_dir + "/singlem/{sample}-" + str(num_threads) + "threads.benchmark"
    output:
        report=output_dirs_dict['singlem'] + "/singlem/{sample}.profile",
        done=touch(output_dirs_dict['singlem'] + "/singlem/{sample}.profile.done"),
    threads: num_threads
    resources:
        runtime=180,  # 3h ceiling covers larger samples
    log:
        output_dirs_dict['singlem'] + "/logs/singlem/{sample}.log",
    shell:
        'eval "$(pixi shell-hook -e singlem)" && '
        "singlem pipe --threads {threads} -1 {input.r1} -2 {input.r2} "
        "-p {output.report} --metapackage {input.db} &> {log}"
