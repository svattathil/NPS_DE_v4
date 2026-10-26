#!/bin/bash
#SBATCH --job-name=wrapper2_part4
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

## Optional runs for different model options
case $SLURM_ARRAY_TASK_ID in
    1)  
       	Rscript $codedir/03a_run_r2.r \
		--cohort emory --corstr exchangeable --stderr_variants none
        ;;
    2)
        Rscript $codedir/03a_run_r2.r \
		--cohort emory --corstr independence --stderr_variants none
        ;;

    3)
	Rscript $codedir/03a_run_r2.r \
		--cohort emory --scalefix --stderr_variants none
	;;
    4)
	Rscript $codedir/03a_run_r2.r \
		--cohort emory --xtransform rankint --stderr_variants none
	;;
esac
