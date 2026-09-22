#!/bin/bash
#SBATCH --job-name=wrapper03_2
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --array=1-9 
#SBATCH --cpus-per-task=1
#SBATCH --time=0-10:00:00 
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

case $SLURM_ARRAY_TASK_ID in
    1)  
        residset=ignoreSVs
	cohort=OHSU
        ;;
    2)
        residset=noProtect
	cohort=OHSU
        ;;
    3)
        residset=resid2_multi
	cohort=OHSU
        ;;
    4)  
        residset=ignoreSVs
	cohort=rush
        ;;
    5)
        residset=noProtect
	cohort=rush
        ;;
    6)
        residset=resid2_multi
	cohort=rush
        ;;
    7)  
        residset=ignoreSVs
	cohort=emory
        ;;
    8)
        residset=noProtect
	cohort=emory
        ;;
    9)
        residset=resid2_multi
	cohort=emory
        ;;
esac



### Multivariate
## Per-cohort
echo Running: $cohort using $residset
Rscript $scriptdir/03a-run_multivariate.r --cohort $cohort --residset $residset

## Meta-analysis
#Rscript $scriptdir/04a-run_mvmeta.r
