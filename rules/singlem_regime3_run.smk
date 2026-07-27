# singlem-regime3 joint method: weebill (sylph fork) --two-stage + singlem pipe
# --no-sylph, combined by `singlem condense --joint`. Ends at
# output_singlem-regime3/singlem-regime3/{sample}.profile (condensed format).
# Metapackage staging is in rules/staging.smk (singlem_regime3_copy_metapackage).
#
# Required names (includer): num_threads, fastq_dir, benchmark_dir,
#   output_dirs_dict, singlem_regime3_metapackage_local, weebill_binary,
#   weebill_db, singlem_regime3_git_base_directory, regime3_weebill_scripts.
# Reads at {fastq_dir}/{sample}.{1,2}.fq.gz.
#
# `unset PYTHONPATH`: the top-level `PYTHONPATH=.. snakemake` puts the repo root
# (which has a `singlem/` submodule dir, no __init__.py) on the path, shadowing
# this branch's editable-installed `singlem`. The regime3 env resolves it itself.
#
# weebill's --two-stage profile loads the whole sylph db into memory, so peak RAM
# scales with the db. The r207 two-stage db (12GB) fits in 16GB; the r232 one is
# 38GB and needs more. Overridable via regime3_weebill_mem_mb (benchmark 7 raises it).
regime3_weebill_mem_mb = globals().get("regime3_weebill_mem_mb", 16000)

rule singlem_regime3_pipe_to_archive:
    input:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        db=singlem_regime3_metapackage_local,
        data_done=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/data/done",
    benchmark:
        benchmark_dir + "/singlem-regime3-pipe/{sample}-" + str(num_threads) + "threads.benchmark"
    output:
        archive=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/{sample}.archive.json",
        done=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/{sample}.archive.done",
    threads: num_threads
    resources:
        runtime=180,  # singlem pipe --no-sylph; 3h ceiling
    log:
        output_dirs_dict['singlem-regime3'] + "/logs/singlem-regime3/{sample}.pipe.log",
    shell:
        'unset PYTHONPATH && eval "$(pixi shell-hook -e singlem-regime3)" && '
        "singlem pipe --threads {threads} -1 {input.r1} -2 {input.r2} "
        "--no-sylph --archive-otu-table {output.archive} --metapackage {input.db} &> {log} && "
        "touch {output.done}"

rule singlem_regime3_weebill_profile:
    input:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
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
        "-1 {input.r1} -2 {input.r2} -o {output.profile} {input.database} -u &> {log} && "
        # weebill writes only a header when nothing is detected; ensure >=1 data row
        "awk 'END {{ exit NR < 2 }}' {output.profile} && touch {output.done}"

rule singlem_regime3_annotate_weebill:
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

rule singlem_regime3_condense_joint:
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
