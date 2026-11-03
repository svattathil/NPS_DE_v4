#!/bin/bash
#SBATCH --job-name=wrap4_part4-5
#SBATCH --nodes=1
#SBATCH --mem=16G
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


### SETUP ###
dir=2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3


### MAIN ###
for cohort in OHSU rush emory
do
    ## Step 4
    Rscript 1_Code/phenotype_check_r4.r \
	    --cohort ${cohort} \
	    --outdir $dir/Perm5d_fl_protectyes_pdom_${cohort}

    ## Step 5
    Rscript 1_Code/summarize_perm_r3.r \
	    --dir $dir/Perm5d_fl_protectyes_pdom_${cohort}

    Rscript 1_Code/summarize_perdomain_r4.r \
	    --dir $dir/Perm5d_fl_protectyes_pdom_${cohort} \
	    --realall $dir/Perm5d_fl_protectyes_pdom_${cohort}_realall0

done
