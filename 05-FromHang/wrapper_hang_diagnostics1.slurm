#!/bin/bash
#SBATCH --job-name=diagn
#SBATCH --nodes=1
#SBATCH --mem=16G 
#SBATCH --cpus-per-task=1
#SBATCH --time=0-8:00:00 
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

## Step 1: quick tests (a few minutes).
##   These only check that everything is set up.
##   If both finish without an error message, go on.
echo "Step 1" 1>&2
Rscript 1_Code/03a_code1_diagnostics.r --cohort emory --maxprot 20
Rscript 1_Code/03a_code5a_permutation.r --cohort emory --B 2 --nprot 30


## Step 2: full diagnostics (about 1 hour).
##   The same analysis as 03a-run_multivariate.r, with extra checks recorded for every protein.
echo "Step 2" 1>&2
time Rscript 1_Code/03a_code1_diagnostics.r --cohort emory --residset ignoreSVs


## Step 3: cleaning check (a few minutes each).
##  Run the first command, then the second.
##  They compare how the residuals were made with a copy of that step,
##  and the second writes a separate new file, so your existing residual file is not touched.
echo "Step 3" 1>&2
Rscript 1_Code/01a_code3_cleaning_variants.r --cohort emory \
	--protect yes --na zero --bookkeeping_only
Rscript 1_Code/01a_code3_cleaning_variants.r --cohort emory --protect yes --na zero


## Step 4: shuffled-data check (about 45 minutes).
##   The outcomes are shuffled among participants on purpose,
##   so no protein should be truly associated; this shows what
##   the Q–Q plot looks like when there is nothing to find.
##   If Step 4 is stopped, for example by a closed laptop,
##   run the same command again; it continues where it stopped.
echo "Step 4" 1>&2
time Rscript 1_Code/03a_code5a_permutation.r --cohort emory --B 10 --nprot 500

