#!/bin/bash
#SBATCH --job-name=enrichment
#SBATCH --nodes=1
#SBATCH --mem=6G 
#SBATCH --cpus-per-task=1
#SBATCH --array=5
#SBATCH --time=0-6:00:00 
#SBATCH --output=Slurmout/%x-%A_%a.o
#SBATCH --error=Slurmout/%x-%A_%a.e
#SBATCH --chdir=.
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=svattathil@health.ucdavis.edu


source ~/.bashrc
conda activate env1


### SETUP ###

case $SLURM_ARRAY_TASK_ID in
    1)  
        runid=OHSU_softpower6_cutheight0.30
        ;;
    2)
        runid=rush_softpower7_cutheight0.30
        ;;
    3)
      	runid=emory_softpower9_cutheight0.30
        ;;
    4)
	network_set=Consensus
	network_group=Auto_softpower7
	runid=auto_softpower7_OHSU
	;;
    5)
	network_set=Consensus
	network_group=Auto_softpower8
	runid=auto_softpower8_OHSU
	;;
    6)
	network_set=Consensus
	network_group=Auto_softpower9
	runid=auto_softpower9_OHSU
	;;
esac

## Directories and files that exist
scriptfile=1_Code/04-Build_networks/04a-test_module_enrichment.r


## Run script
Rscript $scriptfile \
	--network_set $network_set \
	--network_group $network_group \
	--runid $runid 
