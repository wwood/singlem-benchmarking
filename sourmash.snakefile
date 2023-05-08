
import extern
import pandas as pd

output_dir = 'output_115_sourmash'
benchmark_dir = 'benchmarks_115_sourmash'

datasets = extern.run('ls ~/m/msingle/mess/115_camisim_ish_benchmarking/simulated_reads/*.1.fq.gz').split()[:1]
datasets = [os.path.basename(x).split('.')[0] for x in datasets]

original_reads_dir = '/home/woodcrob/m/msingle/mess/115_camisim_ish_benchmarking/simulated_reads'
fastq_dir = output_dir + '/reads'

# kraken_db = output_dir + "/kraken/data/GTDB_release207"
# metapackage = '/home/woodcrob/git/singlem/db/S3.0.5.metapackage20220806.smpkg'
# singlem_metapackage = output_dir + '/singlem/data/S3.0.5.singlem20220806.smpkg'
sourmash_db1 = output_dir + '/sourmash/data/gtdb-rs207.taxonomy.sqldb'
sourmash_db2 = output_dir + '/sourmash/data/gtdb-rs207.genomic-reps.dna.k31.zip'

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
    group: "sample-wise-group"
    input:
        db1='/home/woodcrob/m/msingle/mess/105_novelty_testing/gtdb-rs207.taxonomy.sqldb',
        db2='/home/woodcrob/m/msingle/mess/105_novelty_testing/gtdb-rs207.genomic-reps.dna.k31.zip',
    output:
        db1=sourmash_db1,
        db2=sourmash_db2,
        done=output_dir + "/sourmash/data/done"
    shell:
        "mkdir -p {output_dir}/sourmash/data && cp -rvL {input.db1} {output.db1} && cp -rvL {input.db2} {output.db2} && touch {output.done}"

rule sourmash_run:
    group: "sample-wise-group"
    input:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        db1=sourmash_db1,
        db2=sourmash_db2,
        done=output_dir + "/sourmash/data/done"
    benchmark:
        benchmark_dir + "/sourmash/{sample}.sourmash.benchmark"
    output:
        report=output_dir + "/sourmash/{sample}.gather_gtdbrs207_reps.with-lineages.csv",
        done=output_dir + "/sourmash/{sample}.profile.done"
    conda:
        "sourmash-v4.4.3"
    shell:
        # sourmash tax annotate creates a file with-lineages in the CWD, so we need to cd into the output dir before running it
        "sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge {output_dir}/sourmash/{wildcards.sample} -o {output_dir}/sourmash/{wildcards.sample}.sig {input.r1} {input.r2} && echo 'running gather..' && sourmash gather {output_dir}/sourmash/{wildcards.sample}.sig {input.db2} -o {output_dir}/sourmash/{wildcards.sample}.gather_gtdbrs207_reps.csv && echo 'running tax ..' && cd {output_dir}/sourmash && sourmash tax annotate -g {wildcards.sample}.gather_gtdbrs207_reps.csv -t ../../{input.db1} && cd - && touch {output.done}"

rule sourmash_to_condensed:
    group: "sample-wise-group"
    input:
        output_dir + "/sourmash/{sample}.gather_gtdbrs207_reps.with-lineages.csv",
    output:
        output_dir + "/biobox/{sample}.sourmash_as_condense.csv",
    shell:
        "~/m/msingle/mess/105_novelty_testing/sourmash_to_condensed.py --with-lineages {input} " \
        "--sample {wildcards.sample} " \
        " > {output} "


rule profile_to_biobox:
    group: "sample-wise-group"
    input:
        profile = output_dir + "/biobox/{sample}.sourmash_as_condense.csv",
    params:
        truth = truth_dir + "/{sample}.condensed.biobox",
    output:
        biobox = output_dir + "/biobox/{sample}.biobox",
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
    group: "sample-wise-group"
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
    group: "sample-wise-group"
    input:
        output_dir + "/opal/{sample}.opal_report",
    output:
        output_dir + "/{sample}.finished_all",
    shell:
        "mkdir -p {final_results_directory}/{wildcards.sample}/ && cp -Rv {input} {benchmark_dir} {final_results_directory}/{wildcards.sample}/ && touch {output}"
