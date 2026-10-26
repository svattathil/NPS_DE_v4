#!/bin/bash
#SBATCH --job-name=wrapper2_part2
#SBATCH --nodes=1
#SBATCH --mem=4G 
#SBATCH --cpus-per-task=1
#SBATCH --array=1-4 
#SBATCH --time=0-10:00:00 
#SBATCH --output=Slurmout/%x-%A_%a.o
#SBATCH --error=Slurmout/%x-%A_%a.e
#SBATCH --chdir=.
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=svattathil@health.ucdavis.edu

## Define number of threads
nthreads=$SLURM_JOB_CPUS_PER_NODE


## Load spack packages
#spack load


# Activate conda env1
source ~/.bashrc
conda activate env1


### SETUP ###
codedir=1_Code/05-FromHang/Round2/1_Code/

## Step2
case $SLURM_ARRAY_TASK_ID in
    1)  
       	Rscript $codedir/03a_perm_null_r2.r \
		--cohort emory --design fixedx --B 20 --nprot 500
        ;;
    2)
        Rscript $codedir/03a_perm_null_r2.r \
		--cohort emory --design all --protect yes --B 20 --nprot 500
        ;;

    3)
	Rscript $codedir/03a_perm_null_r2.r \
		--cohort emory --design all --protect no --B 20 --nprot 500
	;;
    4)
	Rscript $codedir/03a_perm_null_r2.r \
		--cohort emory --design batch --protect yes --B 20 --nprot 500
esac
