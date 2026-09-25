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


scriptdir=1_Code/04-Build_networks/
outdir=2_Pipeline/04-Build_networks/


### Test association between consensus modules and latent factors ###
### Submit batch job
sbatch --array=1-3 $scriptdir/06-sbmt_test_latent_factor_association.slurm



## Do gene set enrichment for consensus modules
### Submit batch job
## Specify array tasks in the submit script as necessary
scriptfile=$scriptdir/07-sbmt_test_module_enrichment.slurm
sbatch $scriptfile




