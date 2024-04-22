#!/bin/bash -e

# Run this in the queue using

# ~/m/msingle/mess/124_singlem-benchmarking/1_novel_strains/full_run9$ rm -rf checkout/; mqsub --run-tmp-dir -t1 --days 3 -m 32 -- bash ~/m/msingle/mess/124_singlem-benchmarking/1_novel_strains/run_benchmark_on_lyra.sh

git clone /home/woodcrob/m/msingle/mess/124_singlem-benchmarking/ checkout

cd checkout
ln -s ~/m/msingle/mess/124_singlem-benchmarking/bac120_metadata_r207.tsv
ln -s /home/woodcrob/m/msingle/mess/124_singlem-benchmarking/reference_genomes
ln -s /home/woodcrob/m/msingle/mess/124_singlem-benchmarking/tool_reference_data

cd 1_novel_strains
ln -s /home/woodcrob/m/msingle/mess/124_singlem-benchmarking/1_novel_strains/generated_reads
ln -s /home/woodcrob/m/msingle/mess/124_singlem-benchmarking/1_novel_strains/truths

git log |head -1

./run_benchmark.sh

# Remove copied reference DBs to save disk
rm -rf output_singlem/singlem/data
rm -rf output_kraken/kraken/data/
rm -rf output_sourmash/sourmash/data/
rm -rf output_metabuli/metabuli/data/
rm -rf output_metaphlan/metaphlan/data/
