
import extern
import pandas as pd

output_dir = 'output_115_motus'
benchmark_dir = 'benchmarks_115_motus'

datasets = extern.run('ls ~/m/msingle/mess/115_camisim_ish_benchmarking/simulated_reads/*.1.fq.gz').split()[:1]
datasets = [os.path.basename(x).split('.')[0] for x in datasets]

original_reads_dir = '/home/woodcrob/m/msingle/mess/115_camisim_ish_benchmarking/simulated_reads'
fastq_dir = output_dir + '/reads'

# kraken_db = output_dir + "/kraken/data/GTDB_release207"
# metapackage = '/home/woodcrob/git/singlem/db/S3.0.5.metapackage20220806.smpkg'
# singlem_metapackage = output_dir + '/singlem/data/S3.0.5.singlem20220806.smpkg'
# sourmash_db1 = output_dir + '/sourmash/data/gtdb-rs207.taxonomy.sqldb'
# sourmash_db2 = output_dir + '/sourmash/data/gtdb-rs207.genomic-reps.dna.k31.zip'
motus_db_path_original = '/home/woodcrob/e/motus-v3.0.3/lib/python3.9/site-packages/motus/db_mOTU'
motus_db_path_local = output_dir + '/motus/data/db_mOTU'

truth_dir = '~/m/msingle/mess/115_camisim_ish_benchmarking/truths'

final_results_directory = "~/m/msingle/mess/115_camisim_ish_benchmarking/snakemake_results/"+output_dir

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
        db=directory(motus_db_path_original),
    output:
        db=directory(motus_db_path_local),
        done=output_dir + "/motus/data/done"
    shell:
        "mkdir -p {output.db} && rmdir {output.db} && cp -rvL {input.db} {output.db} && touch {output.done}"

rule motus_run:
    input:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        db=directory(motus_db_path_local),
        done=output_dir + "/motus/data/done"
    benchmark:
        benchmark_dir + "/motus/{sample}.motus.benchmark"
    output:
        report=output_dir + "/motus/{sample}.motus",
        done=output_dir + "/motus/{sample}.profile.done"
    conda:
        "motus-v3.0.3"
    shell:
        "motus profile -db {motus_db_path_local} -f {input.r1} -r {input.r2} -o {output.report} && touch {output.done}"

rule profile_to_biobox:
    input:
        report=output_dir + "/motus/{sample}.motus"
    params:
        truth = truth_dir + "/{sample}.condensed.biobox",
    output:
        biobox = output_dir + "/biobox/{sample}.biobox",
    conda:
        "singlem-dev"
    shell:
        "~/m/msingle/mess/105_novelty_testing/motus_to_condensed.py --motus {input} " \
        "--gtdb ~/m/msingle/mess/115_camisim_ish_benchmarking/motus/mOTUs_3.0.0_GTDB_tax.tsv " \
        " > {output_dir}/biobox/{wildcards.sample}.motus_as_condense.csv " \ 
        "&& ~/git/singlem/extras/condensed_profile_to_biobox.py --input-condensed-table {output_dir}/biobox/{wildcards.sample}.motus_as_condense.csv " \
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
