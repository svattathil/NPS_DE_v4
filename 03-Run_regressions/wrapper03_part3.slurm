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

residset=ignoreSVs


### Multivariate Meta-analysis
Rscript $scriptdir/04a-run_mvmeta.r
