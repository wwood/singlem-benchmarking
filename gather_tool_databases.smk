import os.path
from os.path import join, abspath
import os

output_directory = 'tool_reference_data'
# mkdir
os.makedirs(output_directory, exist_ok=True)

singlem_metapackage = join(output_directory, 'S3.1.1.GTDB_r207.metapackage_20240302.smpkg')
singlem_metapackage_tgz = singlem_metapackage + '.zb.tar.gz'

metaphlan_db = join(output_directory, 'metaphlan_bowtiedb')
metaphlan_index = 'mpa_vOct22_CHOCOPhlAnSGB_202212'

motus_db = join(output_directory, 'motus', 'db_mOTU')
motus_gtdb_tsv = join(output_directory, 'motus', 'mOTUs_3.0.0_GTDB_tax.tsv')

kraken_db = join(output_directory, "struo2-kraken-GTDB_release207")

# Must be at genome level otherwise sourmash tax annotate fails
sourmash_db_taxonomy = join(output_directory, 'sourmash', 'gtdb-rs207.taxonomy.sqldb')
sourmash_db_dna = join(output_directory, 'sourmash', 'gtdb-rs207.genomic-reps.dna.k31.zip')

kaiju_db_original_dir = join(output_directory, 'kaiju')
kaiju_db_original_progenomes_fmi = kaiju_db_original_dir + '/kaiju_db_progenomes.fmi'
kaiju_db_original_progenomes_nodes = kaiju_db_original_dir + '/nodes.dmp'
kaiju_db_original_progenomes_names = kaiju_db_original_dir + '/names.dmp'

map2b_checkout_dir = join(output_directory, 'MAP2B-checkout')
map2b_db = os.path.join(map2b_checkout_dir, 'database/GTDB')

metabuli_db_dir = join(output_directory, 'metabuli') #'metabuli-gtdb207'

tools = ['singlem', 'metaphlan', 'motus', 'kraken', 'sourmash', 'kaiju', 'map2b', 'metabuli']

# debug
tools = ['metaphlan', 'motus', 'kraken', 'sourmash', 'kaiju', 'map2b']

rule all:
    input:
        [join(output_directory, f'{tool}.done') for tool in tools],
        join(output_directory, 'gtdb.done'),
        join(output_directory, 'shadow-genomes.done'),

rule metaphlan:
    output:
        done=touch(join(output_directory, 'metaphlan.done')),
        metaphlan_db=directory(metaphlan_db)
    conda:
        '1_novel_strains/envs/metaphlan.yml'
    log:
        join(output_directory, 'metaphlan.log')
    shell:
        'metaphlan --install --bowtie2db {metaphlan_db} --index {metaphlan_index} &> {log}'

rule kraken_download:
    output:
        done=touch(join(output_directory, 'kraken.done')),
    log:
        abspath(join(output_directory, 'kraken-download.log'))
    shell:
        "mkdir -pv {kraken_db} && cd {kraken_db} && wget --no-directories -r -np -e robots=off http://ftp.tue.mpg.de/ebio/projects/struo2/GTDB_release207/kraken2/ &> {log}"

rule sourmash:
    input:
        join(output_directory, 'sourmash-dna.done'),
        join(output_directory, 'sourmash-tax.done'),
    output:
        done=touch(join(output_directory, 'sourmash.done')),

rule sourmash_sigs:
    output:
        done=touch(join(output_directory, 'sourmash-dna.done')),
        sourmash_dna=sourmash_db_dna,
    log:
        join(output_directory, 'sourmash-sigs-download.log')
    params:
        db = os.path.basename(sourmash_db_dna)
    shell:
        'wget -O {output.sourmash_dna} https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs207/gtdb-rs207.genomic-reps.dna.k31.zip &> {log}'

rule sourmash_tax_download:
    output:
        done=touch(join(output_directory, 'sourmash-tax-download.done')),
        sourmash_taxonomy=join(output_directory, 'sourmash', 'gtdb-rs207.taxonomy.reps.csv.gz'),
    log:
        join(output_directory, 'sourmash-tax-download.log')
    shell:
        'wget -O {output.sourmash_taxonomy} https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs207/gtdb-rs207.taxonomy.reps.csv.gz &> {log}'

rule sourmash_tax_convert:
    input:
        done=join(output_directory, 'sourmash-tax-download.done'),
        sourmash_taxonomy=join(output_directory, 'sourmash', 'gtdb-rs207.taxonomy.reps.csv.gz'),
    output:
        done=touch(join(output_directory, 'sourmash-tax.done')),
        sourmash_taxonomy=sourmash_db_taxonomy,
    conda:
        '1_novel_strains/envs/sourmash.yml'
    log:
        join(output_directory, 'sourmash-tax-convert.log')
    shell:
        'sourmash tax prepare -t {input.sourmash_taxonomy} -o {output.sourmash_taxonomy} -F sql &> {log}'

rule kaiju_download:
    output:
        done=touch(join(output_directory, 'kaiju-download.done')),
    log:
        abspath(join(output_directory, 'kaiju-download.log'))
    shell:
        'mkdir -p {kaiju_db_original_dir} && cd {kaiju_db_original_dir} && wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2023/kaiju_db_progenomes_2023-05-25.tgz -O kaiju_db_progenomes_2023-05-25.tgz &> {log}'

rule kaiju_extract:
    input:
        done=join(output_directory, 'kaiju-download.done'),
    output:
        done=touch(join(output_directory, 'kaiju.done')),
        fmi=kaiju_db_original_progenomes_fmi,
        nodes=kaiju_db_original_progenomes_nodes,
        names=kaiju_db_original_progenomes_names
    log:
        abspath(join(output_directory, 'kaiju-extract.log'))
    params:
        outdir = join(output_directory, 'kaiju')
    shell:
        'cd {params.outdir} && tar -xzf ../kaiju/kaiju_db_progenomes_2023-05-25.tgz &> {log}'

rule map2b_checkout:
    output:
        done=touch(join(output_directory, 'map2b-checkout.done')),
    log:
        join(output_directory, 'map2b-checkout.log')
    shell:
        'rm -rf {map2b_checkout_dir} && git clone --branch MAP2Bv1.5 https://github.com/sunzhengCDNM/MAP2B {map2b_checkout_dir} &> {log}'

rule map2b:
    input:
        join(output_directory, 'map2b-checkout.done'),
    output:
        done=touch(join(output_directory, 'map2b.done')),
        map2b_db=directory(map2b_db)
    conda:
        "1_novel_strains/envs/MAP2B-20230420-conda.yml"
    log:
        join(output_directory, 'map2b.log')
    shell:
        'python3 {map2b_checkout_dir}/scripts/DownloadDB.py -l {map2b_checkout_dir}/config/GTDB.CjePI.database.list -d {map2b_checkout_dir}/database/GTDB &> {log}'

# rule metabuli_download:
#     output:
#         done=touch(join(output_directory, 'metabuli-download.done')),
#     log:
#         join(output_directory, 'metabuli-download.log')
#     shell:
#         'mkdir {metabuli_db_dir} && wget https://metabuli.steineggerlab.workers.dev/gtdb.tar.gz -O {metabuli_db_dir}/gtdb.tar.gz &> {log}'

# rule metabuli_extract:
#     input:
#         done=join(output_directory, 'metabuli-download.done'),
#     output:
#         done=touch(join(output_directory, 'metabuli.done')),
#         db=directory(metabuli_db)
#     log:
#         abspath(join(output_directory, 'metabuli-extract.log'))
#     shell:
#         'cd {output_directory} && tar -xzf metabuli-gtdb.tar.gz &> {log}'

rule motus_db:
    output:
        done=touch(join(output_directory, 'motus-download.done')),
    conda:
        "1_novel_strains/envs/motus.yml"
    log:
        abspath(join(output_directory, 'motus-download.log'))
    shell:
        'motus downloadDB &> {log}'

rule motus_db_move:
    input:
        done=join(output_directory, 'motus-download.done'),
    output:
        motus_db=directory(motus_db),
        done=touch(join(output_directory, 'motus-db.done')),
    conda:
        "1_novel_strains/envs/motus.yml"
    log:
        abspath(join(output_directory, 'motus-db-move.log'))
    shell:
        # ls /mnt/hpccs01/work/microbiome/msingle/mess/124_singlem-benchmarking/.snakemake/conda/b06521a4ea0bdbb4dd2eabbe19701683_/lib/python3.9/site-packages/motus/db_mOTU/
        'bin/migrate_motus_db.py {output.motus_db}'

rule motus_gtdb:
    output:
        done=touch(join(output_directory, 'motus-gtdb.done')),
        motus_gtdb_tsv=motus_gtdb_tsv,
    log:
        join(output_directory, 'motus-gtdb.log')
    shell:
        'wget https://sunagawalab.ethz.ch/share/MOTU_GTDB/mOTUs_3.0.0_GTDB_tax.tsv -O {motus_gtdb_tsv} &> {log}'

rule motus:
    input:
        join(output_directory, 'motus-db.done'),
        join(output_directory, 'motus-gtdb.done'),
    output:
        done=touch(join(output_directory, 'motus.done')),

rule singlem_download:
    output:
        done=touch(join(output_directory, 'singlem.done')),
        singlem_metapackage_tgz=singlem_metapackage_tgz
    conda:
        '1_novel_strains/envs/singlem.yml'
    log:
        join(output_directory, 'singlem-download.log')
    shell:
        "wget 'https://zenodo.org/records/11107165/files/S4.1.0.GTDB_r207.metapackage_20240502.smpkg.zb.tar.gz?download=1' -O {output.singlem_metapackage}.tar.gz &> {log}"

rule singlem_extract:
    input:
        done=join(output_directory, 'singlem.done'),
        singlem_metapackage_tgz=singlem_metapackage_tgz,
    output:
        done=touch(join(output_directory, 'singlem-extract.done')),
        singlem_metapackage=singlem_metapackage
    log:
        abspath(join(output_directory, 'singlem-extract.log'))
    shell:
        'cd {output_directory} && tar -xzf {input.done} &> {log} && mv it {output.singlem_metapackage}'

rule gtdb_download_bac120:
    output:
        done=touch(join(output_directory, 'gtdb_download.done')),
        tar = 'bac120_metadata_r207.tar.gz'
    log:
        join(output_directory, 'gtdb.log')
    shell:
        'wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/bac120_metadata_r207.tar.gz -O {output.tar} &> {log}'

rule gtdb_extract_bac120:
    input:
        done=join(output_directory, 'gtdb_download.done'),
        tar = 'bac120_metadata_r207.tar.gz'
    output:
        done=touch(join(output_directory, 'gtdb.done')),
    log:
        join(output_directory, 'gtdb-extract.log')
    shell:
        'tar -xzf {input.tar} &> {log}'

rule gtdb_download_ar53:
    output:
        done=touch(join(output_directory, 'gtdb_download_ar53.done')),
        tar = 'ar53_metadata_r207.tar.gz'
    log:
        join(output_directory, 'gtdb.log')
    shell:
        'wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/ar53_metadata_r207.tar.gz -O {output.tar} &> {log}'

rule gtdb_extract_ar53:
    input:
        done=join(output_directory, 'gtdb_download_ar53.done'),
        tar = 'ar53_metadata_r207.tar.gz'
    output:
        done=touch(join(output_directory, 'gtdb-ar.done')),
    log:
        join(output_directory, 'gtdb-extract.log')
    shell:
        'tar -xzf {input.tar} &> {log}'

rule shadow_genomes_download:
    output:
        done=touch(join(output_directory, 'shadow-genomes-download.done')),
    log:
        join(output_directory, 'shadow-genomes-download.log')
    shell:
        "rm -f reference_genomes.tar.gz && bash -c 'cd {output_directory} && wget 'https://zenodo.org/records/11002998/files/reference_genomes.tar.gz?download=1' -O reference_genomes.tar.gz' &> {log}"
        
rule shadow_genomes_extract:
    input:
        done=join(output_directory, 'shadow-genomes-download.done'),
    output:
        done=touch(join(output_directory, 'shadow-extract.done')),
    log:
        join(output_directory, 'shadow-genomes.log')
    shell:
        "bash -c 'cd {output_directory} && tar xf reference_genomes.tar.gz' &> {log}"

rule reference_genomes_link:
    input:
        done=join(output_directory, 'shadow-extract.done'),
    output:
        done=touch(join(output_directory, 'shadow-genomes.done')),
    log:
        join(output_directory, 'reference-genomes-link.log')
    shell:
        "bash -c 'cd {output_directory} && ln -s reference_genomes/shadow . && ln -s reference_genomes/new_in_r214 .' &> {log}"

rule genome_pairs_download:
    output:
        done=touch(join(output_directory, 'genome-pairs-download.done')),
    log:
        join(output_directory, 'genome-pairs-download.log')
    shell:
        """
        cd 2_phylogenetic_novelty
        cd genomes
        datasets download genome accession --inputfile ../genome_accessions.txt
        unzip ncbi_dataset.zip

        # Rename files to simple names (e.g. GCA_000508305.1_genomic.fna)
        parallel --col-sep "\t" cp {1} {2} :::: ../genome_ncbi_names.tsv

        cd ../genome_pairs
        datasets download genome accession --inputfile ../genome_pairs_accessions.txt
        unzip ncbi_dataset.zip

        # Rename files to simple names (e.g. GCA_000508305.1_genomic.fna)
        parallel --col-sep "\t" cp {1} {2} :::: ../genome_pairs_ncbi_names.tsv
        """