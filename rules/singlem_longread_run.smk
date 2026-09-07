# SingleM profiling of **long, single-end** reads: reads -> condensed profile.
# The long-read counterpart of rules/singlem_run.smk. Two differences:
#
# 1. Reads are single-end, so only -1 is passed (SingleM's -1/--forward is documented
#    as taking "forward or unpaired" reads, which is the long-read entry point).
# 2. It runs in the `singlem-longread` pixi environment (SingleM 0.21.3), not
#    `singlem` (0.18.0). Long-read support was added in SingleM 0.20.0; 0.18.0 runs
#    on long reads without erroring but recovers substantially less of the community,
#    so using it here would benchmark the version rather than the read technology.
#    Benchmarks 1-10 are unaffected and keep 0.18.0.
#
# Database staging is in rules/staging.smk (rule singlem_copy_metapackage) -- the
# metapackage format is unchanged, so the same staging rule and the same
# `singlem_metapackage_local` path serve both.
#
# Required names (includer): num_threads, fastq_dir, benchmark_dir,
#   output_dirs_dict, singlem_metapackage_local.
# Reads at {fastq_dir}/{sample}.fq.gz.

rule singlem_longread_run_to_profile:
    input:
        reads=fastq_dir + "/{sample}.fq.gz",
        db=singlem_metapackage_local,
        data_done=output_dirs_dict['singlem'] + "/singlem/data/done",
    benchmark:
        benchmark_dir + "/singlem/{sample}-" + str(num_threads) + "threads.benchmark"
    output:
        report=output_dirs_dict['singlem'] + "/singlem/{sample}.profile",
        done=touch(output_dirs_dict['singlem'] + "/singlem/{sample}.profile.done"),
    threads: num_threads
    resources:
        runtime=360,  # long reads mean more bases per sample than the short-read rule
    log:
        output_dirs_dict['singlem'] + "/logs/singlem/{sample}.log",
    shell:
        'eval "$(pixi shell-hook -e singlem-longread)" && '
        "singlem pipe --threads {threads} -1 {input.reads} "
        "-p {output.report} --metapackage {input.db} &> {log}"
