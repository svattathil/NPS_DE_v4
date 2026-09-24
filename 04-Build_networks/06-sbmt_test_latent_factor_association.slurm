#!/bin/bash
#SBATCH --job-name=latent_factor_assoc
#SBATCH --nodes=1
#SBATCH --mem=4G
##SBATCH --array=1-3
#SBATCH --cpus-per-task=1
#SBATCH --time=0-1:00:00 
#SBATCH --output=Slurmout/%x-%A_%a.o
#SBATCH --error=Slurmout/%x-%A_%a.e
#SBATCH --chdir=.
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=svattathil@health.ucdavis.edu


source ~/.bashrc
conda activate env1

## This ended up being the same across all runs
cutheight="0.30"

case $SLURM_ARRAY_TASK_ID in
    1)
        cohort=OHSU
	network_set="cohort-specific"
	softpower=6
	runid=${cohort}_softpower${softpower}_cutheight${cutheight}
        ;;
    2)
        cohort=rush
	network_set="cohort-specific"
	softpower=7
	runid=${cohort}_softpower${softpower}_cutheight${cutheight}
        ;;
    3)
        cohort=emory
	network_set="cohort-specific"
	softpower=9
	runid=${cohort}_softpower${softpower}_cutheight${cutheight}
        ;;
    4)
        cohort=OHSU
	network_set="consensus"
	softpower=8
	runid=${cohort}_consensus_auto_softpower${softpower}
        ;;
    5)
        cohort=rush
	network_set="consensus"
	softpower=8
	runid=${cohort}_consensus_auto_softpower${softpower}
        ;;
    6)
        cohort=emory
	network_set="consensus"
	softpower=8
	runid=${cohort}_consensus_auto_softpower${softpower}
        ;;
esac


## Files that exist
scriptfile=1_Code/04-Build_networks/03a-test_latent_factor_association.rmd


## Files to be created
outputfile=latent_factor_regression_${runid}.html
tmpdir=Tmp/${runid}/


## Define command
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir=paste0(getwd(), "/2_Pipeline/Reports/"),')
cmd+=('output_file="'${outputfile}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=('intermediates_dir="'${tmpdir}'",')
cmd+=(")'")
cmd+=('--args --cohort ' $cohort)
cmd+=('--softpower '$softpower)
cmd+=('--cutheight '$cutheight)
cmd+=('--network_set '$network_set)


## execute command 
eval ${cmd[*]}

