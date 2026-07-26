#!/bin/bash
# Deinterleave CAMI2 marine anonymous_reads.fq.gz (interleaved R1/R2) into
# split_reads/marineN.{1,2}.fq.gz. sample_0 is a 0-byte placeholder, so only
# marine1..marine9 are staged. Reads live on /work (on-disk, ~110 MB/s).
# All 9 run in parallel; pigz -1 (fast) since these are transient inputs read
# once by the tools -- speed matters more than a few % extra size.
set -uo pipefail
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking/3_cami2_marine
eval "$(pixi shell-hook -e art 2>/dev/null)"   # for pigz
BASE=/work/microbiome/msingle/mess/114_cami2_benchmarking/data/marine/short_read/simulation_short_read
mkdir -p split_reads
stage_one() {
  local n=$1
  local src="$BASE/2018.08.15_09.49.32_sample_${n}/reads/anonymous_reads.fq.gz"
  local r1="split_reads/marine${n}.1.fq.gz"
  local r2="split_reads/marine${n}.2.fq.gz"
  if [[ -s "$r1" && -s "$r2" ]]; then echo "marine${n}: already staged"; return; fi
  echo "marine${n}: deinterleaving from $src"
  zcat "$src" \
    | paste - - - - - - - - \
    | tee >(cut -f1-4 | tr '\t' '\n' | pigz -1 -p 2 > "${r1}.tmp") \
    | cut -f5-8 | tr '\t' '\n' | pigz -1 -p 2 > "${r2}.tmp"
  mv "${r1}.tmp" "$r1"; mv "${r2}.tmp" "$r2"
  echo "marine${n}: done ($(stat -c %s "$r1") + $(stat -c %s "$r2") bytes)"
}
# all 9 in flight at once (login node has 40 cores, mostly idle)
for n in 1 2 3 4 5 6 7 8 9; do stage_one "$n" & done
wait
echo "ALL_READS_STAGED $(date)"
ls -la split_reads/