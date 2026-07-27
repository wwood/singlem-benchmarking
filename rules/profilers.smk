# Shared per-tool profiling rules for singlem, metaphlan, sylph and
# singlem-regime3. Each tool ends at output_<tool>/<tool>/<sample>.profile
# (singlem condensed format), which rules/common.smk then turns into a biobox
# and scores with OPAL.
#
# Factored out of the individual benchmark Snakefiles to remove duplication;
# `include: "../rules/profilers.smk"` after defining the names below.
#
# Required names (defined by the includer, e.g. via setup.py + tool_reference_data):
#   num_threads, fastq_dir, benchmark_dir, output_dirs_dict
#   singlem_metapackage, singlem_metapackage_local
#   metaphlan_db, metaphlan_db_local1, metaphlan_index
#   sylph_db, sylph_db_local
#   weebill_binary, weebill_db
#   singlem_regime3_metapackage_local, singlem_regime3_git_base_directory,
#   regime3_weebill_scripts
#
# Reads are expected at {fastq_dir}/{sample}.{1,2}.fq.gz.


###############################################################################
######### singlem
###############################################################################

rule singlem_copy_metapackage:
    input:
        singlem_metapackage,
    output:
        db=directory(singlem_metapackage_local),
        done=output_dirs_dict['singlem'] + "/singlem/data/done",
    resources:
        runtime=60,  # trivial cp; 1h so it isn't stuck at the 48h default
    shell:
        "cp -rL {input} {output.db} && touch {output.done}"

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


###############################################################################
######### metaphlan
###############################################################################

rule metaphlan_copy_db:
    output:
        done=touch(output_dirs_dict['metaphlan'] + "/metaphlan/data/done"),
    resources:
        runtime=60,  # db cp; 1h so it isn't stuck at the 48h default
    shell:
        "cp -rvL {metaphlan_db} {metaphlan_db_local1}"

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


###############################################################################
######### sylph (standard bioconda, single-stage `sylph profile`)
###############################################################################

rule sylph_copy_db:
    input:
        db=sylph_db,
    output:
        db=sylph_db_local,
        done=touch(output_dirs_dict['sylph'] + "/sylph/data/done"),
    resources:
        runtime=60,  # db cp; 1h so it isn't stuck at the 48h default
    shell:
        "mkdir -p $(dirname {output.db}) && cp -r {input.db} {output.db}"

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


###############################################################################
######### singlem-regime3 (weebill --two-stage + singlem pipe + condense --joint)
###############################################################################

rule singlem_regime3_copy_metapackage:
    input:
        singlem_metapackage,
    output:
        db=directory(singlem_regime3_metapackage_local),
        done=output_dirs_dict['singlem-regime3'] + "/singlem-regime3/data/done",
    resources:
        runtime=60,  # trivial cp; 1h so it isn't stuck at the 48h default
    shell:
        "cp -rL {input} {output.db} && touch {output.done}"

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
        mem_mb=16000,
        runtime=120,  # weebill profile is fast; 2h is ample
    log:
        output_dirs_dict['singlem-regime3'] + "/logs/singlem-regime3/{sample}.weebill.log",
    shell:
        "{input.weebill} profile --two-stage -t {threads} -c 100 "
        "-1 {input.r1} -2 {input.r2} -o {output.profile} {input.database} &> {log} && "
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
        "--profile {input.profile} --metapackage {input.db} --output {output.annotated} "
        " &> {log}"

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
