#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --nodes=1
#SBATCH --mem=4G 
#SBATCH --cpus-per-task=1
#SBATCH --array=1 
#SBATCH --time=0-1:00:00 
#SBATCH --output=Slurmout/%x-%A_%a.o
#SBATCH --error=Slurmout/%x-%A_%a.e
#SBATCH --chdir=.
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=svattathil@health.ucdavis.edu

## Define number of threads
nthreads=$SLURM_JOB_CPUS_PER_NODE


# Activate conda env1
source ~/.bashrc
conda activate env1

scriptdir=1_Code/03-Run_regressions/


## Functions
run_regressions01() {
    ## Run the 01a script to run per-domain regressions for each cohort in turn
    ## Bail out if run fails on any cohort
    local run=$1
    local acohort

    for acohort in OHSU rush emory
    do
        cmd=(Rscript 1_Code/03-Run_regressions/01a-limma_nps_DE.r)
        cmd+=(--cohort $acohort)
        cmd+=(--run $run)
        echo ${cmd[*]}
        eval ${cmd[*]} 2> 2_Pipeline/03-Run_regressions/Prot_perDomain/${acohort}_${run}.log
	if [[ $? -ne 0 ]]; then
            echo "Failed on $acohort (run: $run) — see log at 2_Pipeline/03-Run_regressions/Prot_perDomain/${acohort}_${run}.log"
            return 1
        fi
    done
}


## Basic model
run_regressions01 basic
source $scriptdir/02a-run-metafor.r --run basic
source $scriptdir/02a-run-metafor.r --run cond_ind


## Sensitivity - Severity
run_regressions01 severity



