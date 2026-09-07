# Shared synthetic-community **long**-read generation via Badread.
#
# The long-read counterpart of rules/data_generation.smk: same inputs (a per-sample
# coverage file plus a genome->fasta list, taxonomy from the GTDB metadata) and same
# outputs (truth condensed profile + per-genome coverage record), but reads come from
# Badread and are SINGLE-ENDED -- one {sample}.fq.gz, not {sample}.{1,2}.fq.gz. That
# is why this cannot be a flag on the shared ART rule: every downstream rule's input
# signature changes, which is what rules/*_longread_run.smk exist for.
#
# Read technology is chosen by the includer via `longread_read_tech`:
#   'nanopore'    -> Badread's defaults, i.e. Oxford Nanopore R10.4.1
#   'pacbio-hifi' -> pacbio2021 models with a normal qscore distribution
# See bin/generate_longread_community.py for the exact parameters and why they are
# passed explicitly rather than left to Badread's defaults.
#
# Required names (defined by the includer):
#   fastq_dir, truth_dir, num_threads
#   longread_read_tech          'nanopore' or 'pacbio-hifi'
#   coverage_definitions_folder folder of <sample>.tsv coverage files
#   genome_fasta_paths          TSV of <genome>\t<fasta path>
#   gtdb_bac_metadata, gtdb_ar_metadata  (from tool_reference_data)
# Optional:
#   longread_seed               base RNG seed (default 1); genome i uses seed+i, so
#                               a rerun reproduces the same reads
#   longread_runtime            queue minutes for the rule (default 720). Badread is
#                               single-threaded per genome and much slower than ART
#                               (~0.4 Mbp/s), so a 1000-genome community needs hours
#                               rather than the ART rule's 120 minutes.
#
# generate_longread_community.py resolves the GTDB metadata paths relative to CWD, so
# the includer's shell.prefix must cd into the benchmark folder (one level under the
# repo root), exactly as for the ART rule.

longread_seed = globals().get("longread_seed", 1)
longread_runtime = globals().get("longread_runtime", 720)


rule generate_longread_community_and_reads:
    output:
        reads=fastq_dir + "/{sample}.fq.gz",
        condensed=truth_dir + "/{sample}.condensed",
        genomewise=truth_dir + "/{sample}.genomewise.csv",
        done=touch(truth_dir + "/{sample}.finished"),
        done2=touch(fastq_dir + "/{sample}.finished"),
    input:
        coverage=coverage_definitions_folder + "/{sample}.tsv",
        genome_list=genome_fasta_paths,
    params:
        # The script chdirs into a tmpdir before running Badread, so hand it absolute
        # output paths (resolved against snakemake's working dir).
        reads_abs=lambda wildcards, output: os.path.abspath(output.reads),
        condensed_abs=lambda wildcards, output: os.path.abspath(output.condensed),
        genomewise_abs=lambda wildcards, output: os.path.abspath(output.genomewise),
    threads: num_threads
    resources:
        runtime=longread_runtime,
    log:
        truth_dir + "/{sample}.readgen.log",
    shell:
        'eval "$(pixi shell-hook -e badread)" && '
        "mkdir -p {fastq_dir} {truth_dir} && "
        "{workflow.basedir}/../bin/generate_longread_community.py "
        "--read-tech {longread_read_tech} --threads {threads} --seed {longread_seed} "
        "--coverage-file {input.coverage} "
        "--gtdb-bac-metadata {gtdb_bac_metadata} --gtdb-ar-metadata {gtdb_ar_metadata} "
        "--genome-list {input.genome_list} "
        "--output-condensed {params.condensed_abs} "
        "--output-genomewise-coverage {params.genomewise_abs} "
        "-r {params.reads_abs} &> {log}"
