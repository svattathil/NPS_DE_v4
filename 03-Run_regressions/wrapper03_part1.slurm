#!/bin/bash
#SBATCH --job-name=wrapper03_1
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --array=2-12
#SBATCH --cpus-per-task=1
#SBATCH --time=0-2:00:00 
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
run_limma_script01() {
    ## Run the 01a script to run per-domain regressions for each cohort in turn
    ## Bail out if run fails on any cohort
    ## analysis and svchoice are required arguments
    local analysis=$1
    shift
    local svchoice=$1
    shift
    local cohorts=("$@")
    local acohort

    # default to all three cohorts if none specified
    if [[ ${#cohorts[@]} -eq 0 ]]; then
        cohorts=(OHSU rush emory)
    fi

    for acohort in "${cohorts[@]}"
    do
	runid=${analysis}_${svchoice}

	logfile=2_Pipeline/03-Run_regressions/Prot_perDomain/"${runid^}"_${acohort}.log
        cmd=(Rscript 1_Code/03-Run_regressions/01a-limma_nps_DE.r)
	cmd+=(--robust)
	cmd+=(--trend)
	cmd+=(--analysis $analysis)
	cmd+=(--svchoice $svchoice)
        cmd+=(--cohort $acohort)
	
        
        echo ${cmd[*]}
        eval ${cmd[*]} 2> $logfile
	if [[ $? -ne 0 ]]; then
            echo "Failed on $acohort (run: $runid) — see log at $logfile"
            return 1
        fi
    done
}


## The first set of tasks includes SVs as covariates
## The second set does not
if (( SLURM_ARRAY_TASK_ID >= 1 && SLURM_ARRAY_TASK_ID <= 6 )); then
    svchoice="withSV"
elif (( SLURM_ARRAY_TASK_ID >= 7 && SLURM_ARRAY_TASK_ID <= 12 )); then
    svchoice="noSV"
else
    echo "Error: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} is outside 1-12" >&2
    exit 1
fi


## 1. Run limma
case $SLURM_ARRAY_TASK_ID in
    1|7)
	analysis=basic
	run_limma_script01 $analysis $svchoice
	;;
    2|8)
	analysis=severity
	run_limma_script01 $analysis $svchoice OHSU emory
        ;;
    
    3|9)
	analysis=males
	run_limma_script01 $analysis $svchoice
        ;;
    4|10)
	analysis=females
	run_limma_script01 $analysis $svchoice
        ;;

    5|11)
	analysis=e4_adj
	run_limma_script01 $analysis $svchoice
        ;;
    6|12)
	analysis=latfactor
	run_limma_script01 $analysis $svchoice
        ;;
esac


## 2. Run meta-analysis and aggregate statistics to Excel workbook
if [[ $analysis != "latfactor" ]]; then
    ## Don't run cross-cohort meta-analysis for latent factors
    Rscript $scriptdir/02a-run-metafor.r          --analysis $analysis --svchoice $svchoice
    Rscript $scriptdir/aggregate_stats_to_excel.r --analysis $analysis --svchoice $svchoice
fi


if [[ $analysis == "basic" ]]; then
    ## The cond_ind analysis is performed inside the basic analysis -- 
    ## do the meta-analysis and stats aggregation for those results here
    Rscript $scriptdir/02a-run-metafor.r          --analysis cond_ind --svchoice $svchoice
    Rscript $scriptdir/aggregate_stats_to_excel.r --analysis cond_ind --svchoice $svchoice
fi

