#!/bin/bash
#SBATCH --job-name=wrapper03_2
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --array=1-3 
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

case $SLURM_ARRAY_TASK_ID in
    1)  
        residset=ignoreSVs
        ;;
    2)
        residset=noProtect
        ;;
    3)
        residset=resid2_multi
        ;;
esac



### Multivariate
## Per-cohort
for cohort in OHSU rush emory
do
	echo Running: $cohort using $residset
	Rscript $scriptdir/03a-run_multivariate.r --cohort $cohort --residset $residset

done

## Meta-analysis
#Rscript $scriptdir/04a-run_mvmeta.r
