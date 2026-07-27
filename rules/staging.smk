# Shared database-staging rules: copy each tool's reference database from the
# read-only shared location (../tool_reference_data, defined in
# tool_reference_data.py) into a per-run local copy that the profiling rules
# consume. These `cp`/index rules are identical across benchmarks, so they were
# previously copy-pasted into every Snakefile.
#
# Opt-in model: each rule is only defined if its tool is listed in `staged_tools`
# (a set the includer defines *before* the include). A benchmark therefore only
# ever stages the databases it actually uses, and only needs to define the
# *_local / *_db path variables for those tools. This avoids NameErrors at parse
# time for tools a given benchmark does not run.
#
# Required names (defined by the includer, per tool it stages):
#   staged_tools                      set/list of tool names to stage
#   output_dirs_dict                  {tool: output_<tool>} (from *_setup.py)
#   singlem:        singlem_metapackage, singlem_metapackage_local
#   metaphlan:      metaphlan_db, metaphlan_db_local1
#   metaphlan42:    metaphlan42_db, metaphlan42_db_local1
#   sylph:          sylph_db, sylph_db_local
#   kraken:         kraken_db, kraken_db_local
#   sourmash:       sourmash_db_taxonomy(_local), sourmash_db_dna(_local)
#   metabuli:       metabuli_db_dir, metabuli_db_local
#   coverm:         coverm_reference_fna  (builds a strobealign .sti index)
#   singlem-regime3:singlem_metapackage, singlem_regime3_metapackage_local

if "singlem" in staged_tools:

    rule singlem_copy_metapackage:
        input:
            singlem_metapackage,
        output:
            db=directory(singlem_metapackage_local),
            done=output_dirs_dict["singlem"] + "/singlem/data/done",
        resources:
            runtime=60,  # trivial cp; 1h so it isn't stuck at the 48h default
        shell:
            "cp -rL {input} {output.db} && touch {output.done}"


if "metaphlan" in staged_tools:

    rule metaphlan_copy_db:
        output:
            done=touch(output_dirs_dict["metaphlan"] + "/metaphlan/data/done"),
        resources:
            runtime=60,  # db cp; 1h so it isn't stuck at the 48h default
        shell:
            "cp -rvL {metaphlan_db} {metaphlan_db_local1}"


if "metaphlan42" in staged_tools:

    rule metaphlan42_copy_db:
        output:
            done=touch(output_dirs_dict["metaphlan42"] + "/metaphlan42/data/done"),
        resources:
            runtime=60,
        shell:
            "cp -rvL {metaphlan42_db} {metaphlan42_db_local1}"


if "sylph" in staged_tools:

    rule sylph_copy_db:
        input:
            db=sylph_db,
        output:
            db=sylph_db_local,
            done=touch(output_dirs_dict["sylph"] + "/sylph/data/done"),
        resources:
            runtime=60,  # db cp; 1h so it isn't stuck at the 48h default
        shell:
            "mkdir -p $(dirname {output.db}) && cp -r {input.db} {output.db}"


if "kraken" in staged_tools:

    rule braken_copy_data:
        output:
            db=directory(kraken_db_local),
            done=output_dirs_dict["kraken"] + "/kraken/data/done",
        resources:
            runtime=60,
        shell:
            "cp -r {kraken_db} {output.db}/ && touch {output.done}"


if "sourmash" in staged_tools:

    rule sourmash_copy_db:
        input:
            db1=sourmash_db_taxonomy,
            db2=sourmash_db_dna,
        output:
            db1=sourmash_db_taxonomy_local,
            db2=sourmash_db_dna_local,
            done=output_dirs_dict["sourmash"] + "/sourmash/data/done",
        params:
            output_dir=output_dirs_dict["sourmash"],
        resources:
            runtime=60,
        shell:
            "mkdir -p {params.output_dir}/sourmash/data && "
            "cp -rvL {input.db1} {output.db1} && cp -rvL {input.db2} {output.db2} && "
            "touch {output.done}"


if "metabuli" in staged_tools:

    rule metabuli_copy_db:
        input:
            db=metabuli_db_dir,
        output:
            db=directory(metabuli_db_local),
            done=touch(output_dirs_dict["metabuli"] + "/metabuli/data/done"),
        resources:
            runtime=60,
        shell:
            "cp -r {input.db} {output.db}"


if "coverm" in staged_tools:

    # Build the strobealign index of the (large) GTDB genome FASTA once and reuse
    # it across all samples; indexing is the expensive part, so re-indexing per
    # sample would dominate runtime.
    rule coverm_strobealign_index:
        input:
            ref=coverm_reference_fna,
        output:
            index=coverm_reference_fna + ".sti",
        threads: num_threads
        resources:
            mem_mb=500000,
            runtime=180,  # one-time strobealign index of the ~219GB reference
        log:
            output_dirs_dict["coverm"] + "/logs/coverm/strobealign_index.log",
        shell:
            'eval "$(pixi shell-hook -e coverm)" && '
            "strobealign --create-index -t {threads} {input.ref} &> {log}"


if "singlem-regime3" in staged_tools:

    rule singlem_regime3_copy_metapackage:
        input:
            singlem_metapackage,
        output:
            db=directory(singlem_regime3_metapackage_local),
            done=output_dirs_dict["singlem-regime3"] + "/singlem-regime3/data/done",
        resources:
            runtime=60,  # trivial cp; 1h so it isn't stuck at the 48h default
        shell:
            "cp -rL {input} {output.db} && touch {output.done}"
