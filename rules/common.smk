# Shared, tool-agnostic rules: truth/tool -> biobox conversion and OPAL scoring.
#
# These rules were previously copy-pasted into every benchmark Snakefile
# (1_novel_strains, 2_phylogenetic_novelty, 3_cami2_marine, ...). They are
# factored out here so a benchmark can `include: "../rules/common.smk"` instead.
#
# The includer must define these names *before* the include statement:
#   output_prefix                    e.g. "output_"
#   truth_dir                        directory holding <sample>.condensed[.biobox]
#   tools_with_filled_output_profiles  (from tool_reference_data)
# and is expected to have set `shell.prefix` to cd into the workflow dir so the
# ../bin helper scripts and pixi manifest resolve on cluster nodes.


rule truth_condensed_to_biobox:
    input:
        condensed=truth_dir + "/{sample}.condensed",
    output:
        biobox=truth_dir + "/{sample}.condensed.biobox",
    resources:
        runtime=60,  # trivial table conversion
    shell:
        'eval "$(pixi shell-hook -e singlem)" && '
        "{workflow.basedir}/../bin/condensed_profile_to_biobox.py "
        "--input-condensed-table {input.condensed} --output-biobox {output.biobox}"


def get_condensed_to_biobox_extra_args(tool):
    # Tools whose condensed profiles are already filled (every ancestor row
    # present) must not be re-filled when converting to biobox.
    if tool in tools_with_filled_output_profiles:
        return " --no-fill"
    else:
        return ""


rule tool_condensed_to_biobox:
    input:
        profile=output_prefix + "{tool}/{tool}/{sample}.profile",
        truth=truth_dir + "/{sample}.condensed.biobox",
    params:
        extra_args=lambda wildcards: get_condensed_to_biobox_extra_args(wildcards.tool),
    output:
        biobox=output_prefix + "{tool}/biobox/{sample}.biobox",
    resources:
        runtime=60,  # trivial table conversion
    shell:
        'eval "$(pixi shell-hook -e singlem)" && '
        "{workflow.basedir}/../bin/condensed_profile_to_biobox.py {params.extra_args} "
        "--input-condensed-table {input.profile} "
        "--output-biobox {output.biobox} --template-biobox {input.truth}"


rule opal:
    input:
        biobox=output_prefix + "{tool}/biobox/{sample}.biobox",
        truth=truth_dir + "/{sample}.condensed.biobox",
    params:
        output_opal_dir=output_prefix + "{tool}/opal/{sample}.opal_output_directory",
    output:
        report=output_prefix + "{tool}/opal/{sample}.opal_report",
        done=output_prefix + "{tool}/opal/{sample}.opal_report.done",
    resources:
        runtime=60,  # opal is quick; request 1h (not the 48h default) for backfill priority
    shell:
        'eval "$(pixi shell-hook -e opal)" && '
        "opal.py -g {input.truth} -o {params.output_opal_dir} {input.biobox} "
        "|| echo 'expected opal non-zero existatus'; "
        "mv {params.output_opal_dir}/results.tsv {output.report} && "
        "rm -rf {params.output_opal_dir} && touch {output.done}"
