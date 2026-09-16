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

### 1. Build networks
### Run in interactive session
scriptfile=$scriptdir/01a-run_wgcna.rmd

### -- For OHSU ###
cohort=OHSU
softpower=6

outputfile=wgcna_report_${cohort}_softpower$softpower.html
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir=paste0(getwd(), "/2_Pipeline/Reports/"),')
cmd+=('output_file="'${outputfile}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")
cmd+=('--args --cohort ' $cohort)
cmd+=('--softpower '$softpower)

## execute command
eval ${cmd[*]}

### -- For rush ### 
cohort=rush
softpower=7

outputfile=wgcna_report_${cohort}_softpower$softpower.html
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir=paste0(getwd(), "/2_Pipeline/Reports/"),')
cmd+=('output_file="'${outputfile}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")
cmd+=('--args --cohort ' $cohort)
cmd+=('--softpower '$softpower)


## execute command
eval ${cmd[*]}


### -- For emory ###
cohort=emory
softpower=9

outputfile=wgcna_report_${cohort}_softpower$softpower.html
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir=paste0(getwd(), "/2_Pipeline/Reports/"),')
cmd+=('output_file="'${outputfile}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")
cmd+=('--args --cohort ' $cohort)
cmd+=('--softpower '$softpower)


## execute command
eval ${cmd[*]}


: <<EOF
 
### There is no script 02 ###


### 3. Test association with latent factors ###
### Submit batch job
jobid_3=$(sbatch --parsable $scriptdir/03-sbmt_test_latent_factor_association.slurm)


### 4. Test module enrichment ###
### Submit batch job
sbatch --dependency=afterok:$jobid_3 $scriptdir/04-sbmt_test_module_enrichment.slurm



### 5. Aggregate enrichment and association results ###
### Run in interactive session
Rscript $scriptdir/05a-aggregate_enrichment_and_association_pvalues.r


## 6. Test eigengene associations with multivariate outcomes
### Run in interactive session
for cohort in OHSU rush emory
do Rscript 1_Code/04-Build_networks/06a-run_regressions_multivariate.r \
	   --cohort $cohort
done



## 7. Prepare to build consensus networks
### Run in interactive session
scriptfile=$scriptdir/01a-run_wgcna.rmd

### -- For OHSU ###
cohort=OHSU
softpower=9

outputfile=wgcna_for_consensus_module${cohort}_softpower$softpower.html
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir=paste0(getwd(), "/2_Pipeline/Reports/"),')
cmd+=('output_file="'${outputfile}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")
cmd+=('--args --cohort ' $cohort)
cmd+=('--softpower '$softpower)
cmd+=('--network_set ' consensus)

## execute command
eval ${cmd[*]}


### -- For Rush ###
cohort=rush
softpower=7

outputfile=wgcna_for_consensus_module${cohort}_softpower$softpower.html
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir=paste0(getwd(), "/2_Pipeline/Reports/"),')
cmd+=('output_file="'${outputfile}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")
cmd+=('--args --cohort ' $cohort)
cmd+=('--softpower '$softpower)
cmd+=('--network_set ' consensus)

## execute command
eval ${cmd[*]}


### -- For Emory ###
cohort=emory
softpower=9

outputfile=wgcna_for_consensus_module${cohort}_softpower$softpower.html
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir=paste0(getwd(), "/2_Pipeline/Reports/"),')
cmd+=('output_file="'${outputfile}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")
cmd+=('--args --cohort ' $cohort)
cmd+=('--softpower '$softpower)
cmd+=('--network_set ' consensus)

## execute command
eval ${cmd[*]}


## 8. Build consensus network using one-step function
### Run in interactive session
scriptfile9=$scriptdir/09-build_consensus_network.r
scriptfile10=$scriptdir/10-consensus_module_nps_association_meta-analysis.r
scriptfile11=$scriptdir/11-extract_significant_consensus_module_NPS_associations.r

for softpower in 7 8 9
do
    ## Build consensus network and run per-cohort NPS association
    Rscript $scriptfile9 \
	    --auto_softpower $softpower

    ## Do meta-analysis for NPS association
    Rscript $scriptfile10 \
	    --auto_softpower $softpower

    ## Extract significant associations
    Rscript $scriptfile11 \
	    --auto_softpower $softpower
done


## 9. Do gene set enrichment for consensus modules
### Submit batch job
## Specify array tasks in the submit script as necessary
scriptfile=$scriptdir/04-sbmt_test_module_enrichment.slurm
sbatch $scriptfile


### 3. Test association between consensus modules and latent factors ###
### Submit batch job
sbatch --array=4-6 $scriptdir/03-sbmt_test_latent_factor_association.slurm



EOF
