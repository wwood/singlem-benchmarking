#!/bin/bash -e

# Run this in the queue using

# mqsub -m 250 --run-tmp-dir -t 1 -- run_benchmark_on_lyra.sh

git clone /home/woodcrob/m/msingle/mess/124_singlem-benchmarking/ checkout

cd checkout
ln -s ~/m/msingle/mess/124_singlem-benchmarking/bac120_metadata_r207.tsv
ln -s /home/woodcrob/m/msingle/mess/124_singlem-benchmarking/reference_genomes
ln -s /home/woodcrob/m/msingle/mess/124_singlem-benchmarking/tool_reference_data
ln -s /home/woodcrob/m/msingle/mess/124_singlem-benchmarking/generated_reads
ln -s /home/woodcrob/m/msingle/mess/124_singlem-benchmarking/truths

cd 1_novel_strains
git log |head -1

./run_benchmark.sh

# Remove copied reference DBs to save disk
rm -rf output_singlem/singlem/data
rm -rf output_kraken/kraken/data/
rm -rf output_sourmash/sourmash/data/
rm -rf output_metabuli/metabuli/data/
rm -rf output_metaphlan/metaphlan/data/
