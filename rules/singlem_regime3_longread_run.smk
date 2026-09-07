# singlem-regime3 joint method on **long, single-end** reads: weebill (sylph fork)
# --two-stage + singlem pipe --no-sylph, combined by `singlem condense --joint`.
# The long-read counterpart of rules/singlem_regime3_run.smk. The only changes are
# the read arguments -- `-1 <reads>` for singlem pipe and `-r <reads>` for weebill,
# in place of the paired -1/-2 -- so the two files are otherwise line-for-line the
# same and the same alpha / --joint flags are used.
#
# Unlike `singlem`, no version swap is needed here: the singlem-regime3 environment
# is built from the singlem_sylph_condense_regime submodule, which is already at
# 0.21.3 and so has long-read support.
#
# Metapackage staging is in rules/staging.smk (singlem_regime3_copy_metapackage).
#
# Required names (includer): num_threads, fastq_dir, benchmark_dir,
#   output_dirs_dict, singlem_regime3_metapackage_local, weebill_binary,
#   weebill_db, singlem_regime3_git_base_directory, regime3_weebill_scripts.
# Reads at {fastq_dir}/{sample}.fq.gz.
#
# `unset PYTHONPATH`: the top-level `PYTHONPATH=.. snakemake` puts the repo root
# (which has a `singlem/` submodule dir, no __init__.py) on the path, shadowing this
# branch's editable-installed `singlem`. The regime3 env resolves it itself.
regime3_weebill_mem_mb = globals().get("regime3_weebill_mem_mb", 16000)

rule singlem_regime3_longread_pipe_to_archive:
    input:
        reads=fastq_dir + "/{sample}.fq.gz",
        db=singlem_regime3_metapackage_local,
        data_done=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/data/done",
    benchmark:
        benchmark_dir + "/singlem-regime3-pipe/{sample}-" + str(num_threads) + "threads.benchmark"
    output:
        archive=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/{sample}.archive.json",
        done=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/{sample}.archive.done",
    threads: num_threads
    resources:
        runtime=360,  # long reads mean more bases per sample than the short-read rule
    log:
        output_dirs_dict['singlem-regime3'] + "/logs/singlem-regime3/{sample}.pipe.log",
    shell:
        'unset PYTHONPATH && eval "$(pixi shell-hook -e singlem-regime3)" && '
        "singlem pipe --threads {threads} -1 {input.reads} "
        "--no-sylph --archive-otu-table {output.archive} --metapackage {input.db} &> {log} && "
        "touch {output.done}"

rule singlem_regime3_longread_weebill_profile:
    input:
        reads=fastq_dir + "/{sample}.fq.gz",
        database=weebill_db,
        weebill=weebill_binary,
    benchmark:
        benchmark_dir + "/singlem-regime3-weebill/{sample}-" + str(num_threads) + "threads.benchmark"
    output:
        profile=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/{sample}.weebill.tsv",
        done=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/{sample}.weebill.done",
    threads: num_threads
    resources:
        mem_mb=regime3_weebill_mem_mb,
        runtime=120,  # weebill profile is fast; 2h is ample
    log:
        output_dirs_dict['singlem-regime3'] + "/logs/singlem-regime3/{sample}.weebill.log",
    shell:
        "{input.weebill} profile --two-stage -t {threads} -c 100 "
        "-r {input.reads} -o {output.profile} {input.database} -u &> {log} && "
        # weebill writes only a header when nothing is detected; ensure >=1 data row
        "awk 'END {{ exit NR < 2 }}' {output.profile} && touch {output.done}"

rule singlem_regime3_longread_annotate_weebill:
    input:
        profile=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/{sample}.weebill.tsv",
        db=singlem_regime3_metapackage_local,
    output:
        annotated=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/{sample}.weebill.annotated.tsv",
    resources:
        runtime=60,  # trivial annotation step
    log:
        output_dirs_dict['singlem-regime3'] + "/logs/singlem-regime3/{sample}.annotate.log",
    shell:
        'unset PYTHONPATH && eval "$(pixi shell-hook -e singlem-regime3)" && '
        "PYTHONPATH={singlem_regime3_git_base_directory} python3 {regime3_weebill_scripts}/annotate_weebill.py "
        "--profile {input.profile} --metapackage {input.db} --output {output.annotated} &> {log}"

rule singlem_regime3_longread_condense_joint:
    input:
        archive=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/{sample}.archive.json",
        sylph_profile=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/{sample}.weebill.annotated.tsv",
        db=singlem_regime3_metapackage_local,
    benchmark:
        benchmark_dir + "/singlem-regime3-condense/{sample}-" + str(num_threads) + "threads.benchmark"
    output:
        profile=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/{sample}.profile",
        done=touch(output_dirs_dict['singlem-regime3'] + "/singlem-regime3/{sample}.profile.done"),
    threads: 1
    resources:
        runtime=60,  # condense is quick
    log:
        output_dirs_dict['singlem-regime3'] + "/logs/singlem-regime3/{sample}.condense.log",
    shell:
        'unset PYTHONPATH && eval "$(pixi shell-hook -e singlem-regime3)" && '
        "singlem condense --input-archive-otu-table {input.archive} --metapackage {input.db} "
        "--sylph-profile {input.sylph_profile} --joint --taxonomic-profile {output.profile} "
        "--joint-pin-sylph-species --joint-novel-budget --joint --alpha 1 "
        "&> {log}"
