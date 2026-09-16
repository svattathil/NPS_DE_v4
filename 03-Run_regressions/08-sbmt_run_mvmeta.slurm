#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=4G 
#SBATCH --cpus-per-task=1
#SBATCH --time=0-8:00:00 
#SBATCH --output=Slurmout/%x-%A_%a.o
#SBATCH --error=Slurmout/%x-%A_%a.e
#SBATCH --chdir=.
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=svattathil@health.ucdavis.edu


## Command line arguments
residset=$1  ## resid1 or resid2


## Activate conda env1
source ~/.bashrc
conda activate env1


## Run r script
R --slave --args --residset $residset < 1_Code/03-Run_regressions/08a-run_mvmeta.r
