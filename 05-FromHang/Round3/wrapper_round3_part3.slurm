#!/bin/bash
#SBATCH --job-name=wrap3_part3
#SBATCH --nodes=1
#SBATCH --mem=4G 
#SBATCH --cpus-per-task=1
#SBATCH --time=0-1:00:00 
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

codedir=1_Code/05-FromHang/Round3/1_Code/


## Step 3

fulldir=2_Pipeline/03-Run_regressions/Prot_multivar_diag_R2/R2_IgnoreSVs/ ## for convenience
Rscript $codedir/summarize_perm_r3.r \
	--dir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes \
	--full $fulldir/prot_multivar_ignoreSVs_emory_R2_regression_stats_diag.txt
