#!/bin/bash
#SBATCH --job-name=wrap2_part3
#SBATCH --nodes=1
#SBATCH --mem=4G 
#SBATCH --cpus-per-task=1
#SBATCH --time=0-5:00:00 
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


## Step 3
Rscript $codedir/03a_run_r2.r --cohort emory --residset noProtect --stderr_variants none
Rscript $codedir/03a_run_r2.r --cohort emory --residset regressSVs --stderr_variants none

## Step 4
D=2_Pipeline/03-Run_regressions/Prot_multivar_diag_R2
Rscript $codedir/summarize_calibration_r2.r --out $D/summary/step4 \
--perm 5a=$D/Perm5a_fixedx \
--perm 5b_protect=$D/Perm5b_all_protectyes \
--perm 5b_noprot=$D/Perm5b_all_protectno \
--full base=$D/R2_IgnoreSVs/prot_multivar_ignoreSVs_emory_R2_regression_stats_diag.txt


