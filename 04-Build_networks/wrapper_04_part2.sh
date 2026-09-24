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


: <<EOF
 


## 7. Prepare to build consensus networks
### Run in interactive session
scriptfile=$scriptdir/01a-run_wgcna.rmd
residset=ignoreSVs

### -- For OHSU ###
cohort=OHSU
softpower=9


outputfile=wgcna_for_consensus_module_${cohort}_${residset}_softpower$softpower.html
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir="'${outdir}'",')
cmd+=('output_file="'${outputfile}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")
cmd+=('--args --cohort ' $cohort)
cmd+=('--softpower '$softpower)
cmd+=('--residset 'ignoreSVs)
cmd+=('--network_set ' consensus)

## execute command
eval ${cmd[*]}


### -- For Rush ###
cohort=rush
softpower=9

outputfile=wgcna_for_consensus_module_${cohort}_${residset}_softpower$softpower.html
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir="'${outdir}'",')
cmd+=('output_file="'${outputfile}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")
cmd+=('--args --cohort ' $cohort)
cmd+=('--softpower '$softpower)
cmd+=('--residset 'ignoreSVs)
cmd+=('--network_set ' consensus)


## execute command
eval ${cmd[*]}


### -- For Emory ###
cohort=emory
softpower=9

outputfile=wgcna_for_consensus_module_${cohort}_${residset}_softpower$softpower.html
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir="'${outdir}'",')
cmd+=('output_file="'${outputfile}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")
cmd+=('--args --cohort ' $cohort)
cmd+=('--softpower '$softpower)
cmd+=('--residset 'ignoreSVs)
cmd+=('--network_set ' consensus)


## execute command
eval ${cmd[*]}

EOF

## 2. Make list of shared genes
Rscript $scriptdir/02-make_list_of_shared_genes.r


## 3. Build consensus network using one-step function
### Run in interactive session
scriptfile3=$scriptdir/03-build_consensus_network.r
scriptfile4=$scriptdir/04-consensus_module_nps_association_meta-analysis.r
scriptfile5=$scriptdir/05-extract_significant_consensus_module_NPS_associations.r

residset=ignoreSVs
for softpower in 7 8 9
do
    ## Build consensus network and run per-cohort NPS association
    Rscript $scriptfile3 \
	    --auto_softpower $softpower --residset $residset

    ## Do meta-analysis for NPS association
    Rscript $scriptfile4 \
	    --auto_softpower $softpower

    ## Extract significant associations
    Rscript $scriptfile5 \
	    --auto_softpower $softpower
done


### 3. Test association between consensus modules and latent factors ###
### Submit batch job
sbatch --array=4-6 $scriptdir/06-sbmt_test_latent_factor_association.slurm



## 9. Do gene set enrichment for consensus modules
### Submit batch job
## Specify array tasks in the submit script as necessary
scriptfile=$scriptdir/07-sbmt_test_module_enrichment.slurm
sbatch $scriptfile




