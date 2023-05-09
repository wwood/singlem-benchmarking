
import extern
import pandas as pd

output_dir = 'output_115_bracken'
benchmark_dir = 'benchmarks_115_bracken'

datasets = extern.run('ls ~/m/msingle/mess/115_camisim_ish_benchmarking/simulated_reads/*.1.fq.gz').split()
datasets = [os.path.basename(x).split('.')[0] for x in datasets]

original_reads_dir = '/home/woodcrob/m/msingle/mess/115_camisim_ish_benchmarking/simulated_reads'
fastq_dir = output_dir + '/reads'

kraken_db = output_dir + "/kraken/data/GTDB_release207"
# metapackage = '/home/woodcrob/git/singlem/db/S3.0.5.metapackage20220806.smpkg'
# singlem_metapackage = output_dir + '/singlem/data/S3.0.5.singlem20220806.smpkg'

truth_dir = '~/m/msingle/mess/115_camisim_ish_benchmarking/truths'

final_results_directory = "~/m/msingle/mess/115_camisim_ish_benchmarking/snakemake_results/"+output_dir

rule all:
    input:
        expand(output_dir + "/{sample}.finished_all", sample=datasets) 

rule copy_reads:
    input:
        in1=original_reads_dir + "/{sample}.1.fq.gz",
        in2=original_reads_dir + "/{sample}.2.fq.gz"
    output:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
    shell:
        "pwd && mkdir -pv {fastq_dir} && cp -vL {input} {fastq_dir}/"


rule copy_braken_data:
    input:
        "/home/woodcrob/m/db/struo2/GTDB_release207"
    output:
        db=directory(kraken_db),
        done=output_dir + "/kraken/data/done"
    shell:
        "cp -r {input} {output.db}/ && touch {output.done}"

rule kraken_run:
    input:
        reads_copied1 = fastq_dir + "/{sample}.1.fq.gz",
        reads_copied2 = fastq_dir + "/{sample}.2.fq.gz",
        db=kraken_db,
        copy_braken_data_done=output_dir + "/kraken/data/done"
    benchmark:
        benchmark_dir + "/kraken/{sample}.kraken.benchmark"
    output:
        report=output_dir + "/kraken/{sample}.kraken",
        done=output_dir + "/kraken/{sample}.kraken.done"
    shell:
        "export PATH=~/bioinfo/kraken2/install/:$PATH && " \
        "kraken2 --db {input.db} --threads 1 --output /dev/null --report {output.report} --paired {input.reads_copied1} {input.reads_copied2} && touch {output.done}"

rule braken_run:
    input:
        db=kraken_db,
        kraken_report=output_dir + "/kraken/{sample}.kraken",
        kraken_done=output_dir + "/kraken/{sample}.kraken.done"
    output:
        output_dir + "/braken/{sample}.report.S",
        output_dir + "/braken/{sample}.report.G",
        output_dir + "/braken/{sample}.report.F",
        output_dir + "/braken/{sample}.report.O",
        output_dir + "/braken/{sample}.report.C",
        output_dir + "/braken/{sample}.report.P",
        output_dir + "/braken/{sample}.report.D",
        done=output_dir + "/braken/{sample}.done"
    shell:
        "export PATH=~/bioinfo/Bracken/:$PATH && " \
        "bracken -d {input.db} -r 150 -l S -t 10 -o {output_dir}/braken/{wildcards.sample}.report.S -i {input.kraken_report} && " \
        "bracken -d {input.db} -r 150 -l G -t 10 -o {output_dir}/braken/{wildcards.sample}.report.G -i {input.kraken_report} && " \
        "bracken -d {input.db} -r 150 -l F -t 10 -o {output_dir}/braken/{wildcards.sample}.report.F -i {input.kraken_report} && " \
        "bracken -d {input.db} -r 150 -l O -t 10 -o {output_dir}/braken/{wildcards.sample}.report.O -i {input.kraken_report} && " \
        "bracken -d {input.db} -r 150 -l C -t 10 -o {output_dir}/braken/{wildcards.sample}.report.C -i {input.kraken_report} && " \
        "bracken -d {input.db} -r 150 -l P -t 10 -o {output_dir}/braken/{wildcards.sample}.report.P -i {input.kraken_report} && " \
        "bracken -d {input.db} -r 150 -l D -t 10 -o {output_dir}/braken/{wildcards.sample}.report.D -i {input.kraken_report} && " \
        "touch {output.done}"


rule bracken_to_biobox:
    input:
        output_dir + "/braken/{sample}.report.S",
        output_dir + "/braken/{sample}.report.G",
        output_dir + "/braken/{sample}.report.F",
        output_dir + "/braken/{sample}.report.O",
        output_dir + "/braken/{sample}.report.C",
        output_dir + "/braken/{sample}.report.P",
        output_dir + "/braken/{sample}.report.D",
        output_dir + "/braken/{sample}.done"
    params:
        truth = truth_dir + "/{sample}.condensed.biobox",
    output:
        biobox = output_dir + "/biobox/{sample}.biobox",
    conda:
        "singlem-dev"
    shell:
        # Convert reports to singlem condense format
        "mkdir -p {output_dir}/biobox && " \
        "~/m/msingle/mess/105_novelty_testing/kraken_to_biobox.py --report-prefix {output_dir}/braken/{wildcards.sample}.report " \
        "--bacterial-taxonomy ~/m/db/gtdb/gtdb_release207/bac120_taxonomy_r207.tsv " \
        "--archaeal-taxonomy ~/m/db/gtdb/gtdb_release207/ar53_taxonomy_r207.tsv > {output_dir}/biobox/{wildcards.sample}.braken_as_condense.csv " \ 
        "&& " \
        "~/git/singlem/extras/condensed_profile_to_biobox.py --input-condensed-table {output_dir}/biobox/{wildcards.sample}.braken_as_condense.csv " \
        "--output-biobox {output.biobox} --template-biobox {params.truth} " \
        "--no-fill > {output_dir}/biobox/{wildcards.sample}.biobox"

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