#!/bin/bash
# Split-half limma-voom DE: 15 types x 2 halves per seed. Resumable (fit_half.py skips done fits).
# Usage: run_split_half.sh "<seeds>" [parallel jobs]. Run in tmux session split_half_de.
set -euo pipefail
D=/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/subtype_vs_pooled/split_half
W=/space/scratch/rschwartz_psyrerun_work/split_half
SEEDS=${1:-"1 2 3 4 5 6 7 8 9 10"}; P=${2:-12}
mkdir -p $W/logs
cd $D
[ -s splits.tsv ] || /usr/bin/python3 make_splits.py
for s in $SEEDS; do
  t0=$(date +%s)
  /usr/bin/python3 -c "import common; [print(sch, ct) for _, sch, ct in common.jobs()]" | while read sch ct; do
    for h in A B; do echo "$sch $ct $s $h"; done
  done | xargs -P $P -L 1 bash -c '/usr/bin/python3 '$D'/fit_half.py "$0" "$1" "$2" "$3" > '$W'/logs/$2_$3_$0__$1.log 2>&1 || echo "FAILED $0 $1 $2 $3"'
  echo "seed $s done in $(( $(date +%s) - t0 )) s" | tee -a $W/logs/timing.log
done
echo SPLIT HALF DONE
