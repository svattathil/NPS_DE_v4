#!/bin/bash
#SBATCH --job-name=wrapper03
#SBATCH --nodes=1
#SBATCH --mem=4G 
#SBATCH --cpus-per-task=1
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
    shift
    local cohorts=("$@")
    local acohort

    # default to all three cohorts if none specified
    if [[ ${#cohorts[@]} -eq 0 ]]; then
        cohorts=(OHSU rush emory)
    fi

    for acohort in "${cohorts[@]}"
    do

	logfile=2_Pipeline/03-Run_regressions/Prot_perDomain/"${run^}"/${run}_${acohort}.log
        cmd=(Rscript 1_Code/03-Run_regressions/01a-limma_nps_DE.r)
	cmd+=(--run $run)
        cmd+=(--cohort $acohort)
        
        echo ${cmd[*]}
        eval ${cmd[*]} 2> $logfile
	if [[ $? -ne 0 ]]; then
            echo "Failed on $acohort (run: $run) — see log at $logfile"
            return 1
        fi
    done
}


## Basic model
run_regressions01 basic
Rscript $scriptdir/02a-run-metafor.r --run basic
Rscript $scriptdir/aggregate_stats_to_excel.r --run basic


## Sensitivity - Severity
run_regressions01 severity OHSU emory
Rscript $scriptdir/02a-run-metafor.r --run severity
Rscript $scriptdir/aggregate_stats_to_excel.r --run severity


## Sensitivity - sex-specific
run_regressions01 males
Rscript $scriptdir/aggregate_stats_to_excel.r --run males

run_regressions01 females
Rscript $scriptdir/aggregate_stats_to_excel.r --run females
