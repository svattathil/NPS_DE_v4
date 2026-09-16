#!/bin/bash
#SBATCH --job-name=wrapper04
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


# Activate conda env1
source ~/.bashrc
conda activate env1


scriptdir=1_Code/04-Build_networks/
outdir=2_Pipeline/04-Build_networks/

### 1. Build networks
### Run in interactive session
scriptfile=$scriptdir/01a-run_wgcna.rmd




### -- For OHSU ###
cohort=OHSU
softpower=4


outputfile=wgcna_report_${cohort}_softpower$softpower.html
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir="'${outdir}'",')
cmd+=('output_file="'${outputfile}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")
cmd+=('--args --cohort ' $cohort)
cmd+=('--softpower '$softpower)

## execute command
eval ${cmd[*]}

### -- For rush ### 
cohort=rush
softpower=5

outputfile=wgcna_report_${cohort}_softpower$softpower.html
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir="'${outdir}'",')
cmd+=('output_file="'${outputfile}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")
cmd+=('--args --cohort ' $cohort)
cmd+=('--softpower '$softpower)


## execute command
eval ${cmd[*]}



### -- For emory ###
cohort=emory
softpower=6

outputfile=wgcna_report_${cohort}_softpower$softpower.html
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir="'${outdir}'",')
cmd+=('output_file="'${outputfile}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")
cmd+=('--args --cohort ' $cohort)
cmd+=('--softpower '$softpower)


## execute command
eval ${cmd[*]}

