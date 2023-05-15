#!/bin/bash

git clone /home/woodcrob/m/msingle/mess/124_singlem-benchmarking/ checkout
cd checkout/1_novel_strains
git log |head -1

./run_benchmark.sh
