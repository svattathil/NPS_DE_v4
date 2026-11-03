#!/bin/bash
#SBATCH --job-name=wrap4_part3
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --array=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-10:00:00 
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


### SETUP ###

case $SLURM_ARRAY_TASK_ID in
    1)  
        cohort=OHSU
        ;;
    2)
        cohort=rush
        ;;
    3)
        cohort=emory
        ;;
esac

## Step 3
Rscript 1_Code/03a_perm_null_r3.r \
	--cohort $cohort --design fl --B 1 --nprot 0 \
	--rep_from 0 --rep_to 0 \
	--perdomain \
	--outtag Perm5d_fl_protectyes_pdom_${cohort}_realall0
