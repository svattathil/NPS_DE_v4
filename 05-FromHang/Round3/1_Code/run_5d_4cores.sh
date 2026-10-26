#!/bin/bash
###############################################################################
## run_5d_4cores.sh      (Mac or Linux)   Round 3, design 5d on 4 cores
## Run from the project directory:   bash 1_Code/run_5d_4cores.sh
## It starts 4 R processes, each doing a quarter of the replicates 0..B, and waits until all are finished.
## Windows: do not use this file; use the four-terminal method of README.md (Option A).
## Settings can be changed from the command line, e.g.   B=50 bash 1_Code/run_5d_4cores.sh
###############################################################################
B=${B:-100}
NPROT=${NPROT:-500}
TAGARG=${TAGARG:-}               # e.g. TAGARG="--outtag R3smoke"
SLEEP=${SLEEP:-60}               # seconds the first process gets to write the manifest before the others start
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

mkdir -p logs_5d
n=$((B + 1)); size=$(( (n + 3) / 4 ))
for k in 0 1 2 3; do
    from=$((k * size)); to=$((from + size - 1))
    [ $from -gt $B ] && continue
    [ $to -gt $B ] && to=$B
    echo "starting replicates $from to $to (log: logs_5d/block$((k + 1)).log)"
    Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --protect yes --B $B --nprot $NPROT \
        --rep_from $from --rep_to $to $TAGARG > logs_5d/block$((k + 1)).log 2>&1 &
    [ $k -eq 0 ] && sleep $SLEEP
done
wait
echo "All blocks finished. Check the last lines of logs_5d/block*.log, then run summarize_perm_r3.r (README.md, Step 3)."
