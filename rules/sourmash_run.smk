# sourmash profiling: sketch reads -> gather vs GTDB reps -> tax metagenome
# -> condensed profile. DB staging is in rules/staging.smk (rule sourmash_copy_db).
#
# Required names (includer): num_threads, fastq_dir, benchmark_dir,
#   output_dirs_dict, sourmash_db_taxonomy_local, sourmash_db_dna_local,
#   sourmash_db_dna, sourmash_db_taxonomy.

rule sourmash_run:
    input:
        r1=fastq_dir + "/{sample}.1.fq.gz",
        r2=fastq_dir + "/{sample}.2.fq.gz",
        db1=sourmash_db_taxonomy_local,
        db2=sourmash_db_dna_local,
        done=output_dirs_dict['sourmash'] + "/sourmash/data/done",
    benchmark:
        benchmark_dir + "/sourmash/{sample}-" + str(num_threads) + "threads.benchmark"
    threads: num_threads
    output:
        report=output_dirs_dict['sourmash'] + "/sourmash/{sample}.csv",
        done=touch(output_dirs_dict['sourmash'] + "/sourmash/{sample}.profile.done"),
    params:
        output_dir=output_dirs_dict['sourmash'],
        sourmash_prefix=lambda wildcards: output_dirs_dict['sourmash'] + "/sourmash/" + wildcards.sample,
    log:
        output_dirs_dict['sourmash'] + "/logs/sourmash/{sample}.log",
    shell:
        # sourmash has no --threads option.
        'eval "$(pixi shell-hook -e sourmash)" && '
        "bash -c 'sourmash sketch dna -p k=31,abund --merge blah -o {params.sourmash_prefix}.sig {input.r1} {input.r2} "
        '&& echo "running gather.." '
        "&& sourmash gather -k 31 --dna {params.sourmash_prefix}.sig {sourmash_db_dna} -o {params.sourmash_prefix}.gather_gtdbrs207_reps.csv "
        '&& echo "running tax .." '
        "&& sourmash tax metagenome -g {params.sourmash_prefix}.gather_gtdbrs207_reps.csv -t {sourmash_db_taxonomy} >{output.report}' "
        "&> {log}"

rule sourmash_to_condensed:
    input:
        output_dirs_dict['sourmash'] + "/sourmash/{sample}.csv",
    output:
        profile=output_dirs_dict['sourmash'] + "/sourmash/{sample}.profile",
    resources:
        runtime=60,  # trivial conversion
    shell:
        "{workflow.basedir}/../bin/sourmash_to_condensed.py --summary-csv {input} "
        "--sample {wildcards.sample} > {output}"
