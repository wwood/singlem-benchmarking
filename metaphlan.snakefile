
import extern
import pandas as pd

output_dir = 'output_115_metaphlan'
benchmark_dir = 'benchmarks_115_metaphlan'

datasets = extern.run('ls ~/m/msingle/mess/115_camisim_ish_benchmarking/simulated_reads/*.1.fq.gz').split()[:1]
datasets = [os.path.basename(x).split('.')[0] for x in datasets]

original_reads_dir = '/home/woodcrob/m/msingle/mess/115_camisim_ish_benchmarking/simulated_reads'
fastq_dir = output_dir + '/reads'

metaphlan_db_original1 = '~/m/msingle/mess/115_camisim_ish_benchmarking/metaphlan_bowtiedb'
metaphlan_db_local1 = output_dir + '/metaphlan/data/metaphlan_bowtiedb'

truth_dir = '~/m/msingle/mess/115_camisim_ish_benchmarking/truths'

final_results_directory = "~/m/msingle/mess/115_camisim_ish_benchmarking/snakemake_results/"+output_dir

num_threads = 1

rule all:
    input:
        expand(output_dir + "/{sample}.finished_all", sample=datasets) 


rule copy_reads:
    group: "sample-wise-group"
    input:
        in1=original_reads_dir + "/{sample}.1.fq.gz",
        in2=original_reads_dir + "/{sample}.2.fq.gz",
    output:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
    shell:
        "pwd && mkdir -pv {fastq_dir} && cp -vL {input} {fastq_dir}/"

rule copy_db:
    input:
        # Cannot use the directory as input/output because humann complains when
        # there's a snakemake hidden file in the dir

        # db1=directory(metaphlan_db_original1),
        # db2=directory(metaphlan_db_original2),
    output:
        # db1=directory(metaphlan_db_local1),
        # db2=directory(metaphlan_db_local2),
        done=output_dir + "/metaphlan/data/done"
    shell:
        "mkdir -p {metaphlan_db_local1} && rmdir {metaphlan_db_local1} && cp -rvL {metaphlan_db_original1} {metaphlan_db_local1} && touch {output.done}"

rule profile:
    input:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        done=output_dir + "/metaphlan/data/done"
    benchmark:
        benchmark_dir + "/metaphlan/{sample}.benchmark"
    output:
        sgb_report=output_dir + "/metaphlan/{sample}.sgb_report",
        done=output_dir + "/metaphlan/{sample}.profile.done"
    conda:
        "envs/metaphlan.yml"
    threads: num_threads
    shell:
        # Concatenate input files because metaphlan can't handle multiple input files
        "rm -f {output.sgb_report} {wildcards.sample}.cat.fq.gz.bowtie2out.txt; cat {input.r1} {input.r2} > {wildcards.sample}.cat.fq.gz && metaphlan {wildcards.sample}.cat.fq.gz --nproc {threads} --input_type fastq --bowtie2db {metaphlan_db_local1} -o {output.sgb_report} && touch {output.done}"

rule convert_profile_to_GTDB:
    input:
        report=output_dir + "/metaphlan/{sample}.sgb_report"
    output:
        gtdb_report=output_dir + "/metaphlan/{sample}.profile",
        done=output_dir + "/metaphlan/{sample}.gtdb_report.done"
    conda:
        "envs/metaphlan.yml"
    shell:
        "sgb_to_gtdb_profile.py -i {input.report} -o {output.gtdb_report} -d {metaphlan_db_local1}/mpa_vOct22_CHOCOPhlAnSGB_202212.pkl && touch {output.done}"

rule profile_to_biobox:
    input:
        report=output_dir + "/metaphlan/{sample}.profile"
    params:
        truth = truth_dir + "/{sample}.condensed.biobox",
    output:
        biobox = output_dir + "/biobox/{sample}.biobox",
    conda:
        "singlem-dev"
    shell:
        "~/m/msingle/mess/105_novelty_testing/metaphlan_to_condensed.py --metaphlan {input} --sample {wildcards.sample}" \
        " > {output_dir}/biobox/{wildcards.sample}.metaphlan_as_condense.csv " \ 
        "&& ~/git/singlem/extras/condensed_profile_to_biobox.py --input-condensed-table {output_dir}/biobox/{wildcards.sample}.metaphlan_as_condense.csv " \
        "--output-biobox {output.biobox} --template-biobox {params.truth} " \
        " > {output.biobox}"

rule opal:
    input:
        biobox = output_dir + "/biobox/{sample}.biobox",
    params:
        truth = truth_dir + "/{sample}.condensed.biobox",
    output:
        output_dir + "/opal/{sample}.opal_report",
    conda:
        'cami-opal-v1.0.11'
    shell:
        "opal.py -g {params.truth} -o {output_dir}/opal/{wildcards.sample}.opal_output_directory {input.biobox} || echo 'expected opal non-zero existatus'; mv {output_dir}/opal/{wildcards.sample}.opal_output_directory/results.tsv {output} && rm -rf {output_dir}/opal/{wildcards.sample}.opal_output_directory"

rule copy_out_results:
    input:
        output_dir + "/opal/{sample}.opal_report",
    output:
        output_dir + "/{sample}.finished_all",
    shell:
        "mkdir -p {final_results_directory}/{wildcards.sample}/ && cp -Rv {input} {benchmark_dir} {final_results_directory}/{wildcards.sample}/ && touch {output}"
