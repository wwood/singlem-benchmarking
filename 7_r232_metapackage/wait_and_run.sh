#!/bin/bash -e
# Wait for the r232 metapackage copy to finish, then launch the full benchmark 7
# run on the aqua queue. The metapackage cp writes to a .tmp path and mv's it to
# the final name on success, so "final dir exists AND cp process gone" means done.
cd /mnt/weka/scratch/microbiome/woodcrob/non_sensitive/singlem-benchmarking
FINAL=tool_reference_data/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb

echo "$(date) waiting for metapackage copy to complete .."
while true; do
  if [ -d "$FINAL" ] && ! pgrep -f "cp -rL /container_home/m/db/singlem/S6.5" >/dev/null; then
    break
  fi
  sleep 30
done
echo "$(date) metapackage present: $(du -sh $FINAL | cut -f1). Launching aqua run .."

# Sanity: metapackage must be loadable (v6) before we submit ~dozens of jobs.
PYTHONPATH=.. pixi run --environment singlem python3 -c \
  "from singlem.metapackage import Metapackage; m=Metapackage.acquire('$FINAL'); print('metapackage OK, version', m.version)"

exec 7_r232_metapackage/run.sh
