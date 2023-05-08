
import extern
import pandas as pd

benchmark_dir = 'benchmarks_115'

datasets = extern.run('ls ~/m/msingle/mess/115_camisim_ish_benchmarking/simulated_reads/*.1.fq.gz').split()
datasets = [os.path.basename(x).split('.')[0] for x in datasets]

original_reads_dir = '/home/woodcrob/m/msingle/mess/115_camisim_ish_benchmarking/simulated_reads'
fastq_dir = output_dir + '/reads'

# kraken_db = output_dir + "/kraken/data/GTDB_release207"
# metapackage = '/home/woodcrob/git/singlem/db/S3.0.5.metapackage20220806.smpkg'
# singlem_metapackage = output_dir + '/singlem/data/S3.0.5.singlem20220806.smpkg'
metapackage = '/home/woodcrob/git/singlem/db/S3.1.0.metapackage_20221209.smpkg.uploaded'
singlem_metapackage = output_dir + '/singlem/data/S3.1.0.metapackage_20221209.smpkg'

truth_dir = '~/m/msingle/mess/115_camisim_ish_benchmarking/truths'

final_results_directory = "~/m/msingle/mess/115_camisim_ish_benchmarking/snakemake_results/"+output_dir

num_threads = config['num_threads']

tools = ['singlem', 'metaphlan', 'motus', 'bracken']
tools = ['singlem']

output_prefix = 'output_115_'
output_dirs = list([output_prefix+tool for tool in tools])
output_dirs_dict = dict(zip(tools, output_dirs))

rule all:
    input:
        expand("{output_dir}/{sample}.finished_all", sample=datasets, output_dir=output_dirs)

rule copy_reads:
    input:
        in1=original_reads_dir + "/{sample}.1.fq.gz",
        in2=original_reads_dir + "/{sample}.2.fq.gz"
    output:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
    shell:
        "pwd && mkdir -pv {fastq_dir} && cp -vL {input} {fastq_dir}/"


rule copy_metapackage:
    input:
        metapackage
    output:
        db=directory(singlem_metapackage),
        done=output_prefix + "_singlem/singlem/data/done",
        output_dir = output_prefix+'singlem',
    shell:
        "mkdir -p {output_dir}/singlem/data && cp -rL {input} {output.db} && touch {output.done}"

rule singlem_run_to_profile:
    input:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        db=singlem_metapackage,
        copy_braken_data_done=output_dir + "/singlem/data/done"
    benchmark:
        benchmark_dir + "/singlem/{sample}.singlem.benchmark"
    output:
        report=output_dir + "/singlem/{sample}.profile",
        done=output_dir + "/singlem/{sample}.profile.done"
    conda:
        "singlem-dev"
    threads:
        num_threads
    shell:
        "~/git/singlem/bin/singlem pipe -t {threads} -1 {input.r1} -2 {input.r2} -p {output.report} --metapackage {input.db} && touch {output.done}"


rule profile_to_biobox:
    input:
        profile = output_dirs_dict['singlem'] + "/singlem/{sample}.profile",
    params:
        truth = truth_dir + "/{sample}.condensed.biobox",
    output:
        biobox = output_dirs_dict['singlem'] + "/biobox/{sample}.biobox",
    conda:
        "singlem-dev"
    shell:
        # Convert reports to singlem condense format
        # run opam against the truth
        "mkdir -p {output_dir}/biobox && " \
        "~/git/singlem/extras/condensed_profile_to_biobox.py --input-condensed-table {input.profile} " \
        "--output-biobox {output.biobox} --template-biobox {params.truth} " \
        " > {output.biobox}"

rule opal:
    input:
        biobox = output_dir_prefix+"{tool}/biobox/{sample}.biobox",
    params:
        truth = truth_dir + "/{sample}.condensed.biobox",
        output_dir = output_dir_prefix+"{tool}",
    output:
        output_dir_prefix+"{tool}/opal/{sample}.opal_report",
    conda:
        'envs/opal.yml'
    shell:
        "opal.py -g {params.truth} -o {output_dir}/opal/{wildcards.sample}.opal_output_directory {input.biobox} || echo 'expected opal non-zero existatus'; mv {output_dir}/opal/{wildcards.sample}.opal_output_directory/results.tsv {output} && rm -rf {output_dir}/opal/{wildcards.sample}.opal_output_directory"


# rule copy_out_results:
#     group: "sample-wise-group"
#     input:
#         output_dir + "/opal/{sample}.opal_report",
#     output:
#         output_dir + "/{sample}.finished_all",
#     shell:
#         "mkdir -p {final_results_directory}/{wildcards.sample}/ && cp -Rv {input} {benchmark_dir} {final_results_directory}/{wildcards.sample}/ && touch {output}"