#!/bin/bash
##SBATCH --job-name=limma_DE
#SBATCH --nodes=1
#SBATCH --mem=4G 
#SBATCH --cpus-per-task=1
#SBATCH --time=0-1:00:00 
#SBATCH --output=Slurmout/%x-%A_%a.o
#SBATCH --error=Slurmout/%x-%A_%a.e
#SBATCH --chdir=.
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=svattathil@health.ucdavis.edu

 
## Build command
## Redirect stderr to log file

for acohort in OHSU rush emory
do
    cmd=(Rscript 1_Code/03-Run_regressions/01a-limma_nps_DE.r)
    cmd+=(--cohort $acohort )
    
    echo ${cmd[*]}
    eval ${cmd[*]} 2> 2_Pipeline/03-Run_regressions/Prot_perDomain/$acohort.log
done
