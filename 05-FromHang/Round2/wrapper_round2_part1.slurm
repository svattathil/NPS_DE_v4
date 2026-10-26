#!/bin/bash
#SBATCH --job-name=wrap2_part1
#SBATCH --nodes=1
#SBATCH --mem=4G 
#SBATCH --cpus-per-task=1
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

codedir=1_Code/05-FromHang/Round2/1_Code/

## Step 0
#Rscript $codedir/03a_run_r2.r \
#	--cohort emory --maxprot 20 --outtag R2smoke


## Step 1
Rscript $codedir//03a_run_r2.r \
	--cohort emory --residset ignoreSVs --stderr_variants jack,j1s,fij



