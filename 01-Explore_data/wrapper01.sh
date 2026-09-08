#!/bin/bash
#SBATCH --job-name=wrapper_01
#SBATCH --nodes=1
#SBATCH --mem=4G 
#SBATCH --cpus-per-task=1
#SBATCH --array=1 
#SBATCH --time=0-1:00:00 
#SBATCH --output=Slurmout/%x-%A_%a.o
#SBATCH --error=Slurmout/%x-%A_%a.e
#SBATCH --chdir=.
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=svattathil@ucdavis.edu

scriptdir=1_Code/01-Explore_data


output_dir=./2_Pipeline/01-Explore_data/
## 01
scriptfile=$scriptdir/01-explore_proteomics.rmd
cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir="'${output_dir}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")

eval ${cmd[*]}


## 02
R --slave < $scriptdir/02-regress_batch_from_proteomics.r


## 03
scriptfile=$scriptdir/03-detect_outliers.rmd

cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir="'${output_dir}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")

eval ${cmd[*]}


## 04
scriptfile=$scriptdir/04-extract_phenos_OHSU.rmd

cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir="'${output_dir}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")

eval ${cmd[*]}


## 05a and b
scriptfile=$scriptdir/05a-extract_phenos_rush.rmd

cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir="'${output_dir}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")

eval ${cmd[*]}


R --slave < $scriptdir/05b-extract_phenos_rush.r

 
## 06
scriptfile=$scriptdir/06-extract_phenos_emory.rmd

cmd=(R --slave -e )
cmd+=("'rmarkdown::render(")
cmd+=('"'$scriptfile'",')
cmd+=('output_dir="'${output_dir}'",')
cmd+=('knit_root_dir=getwd(),')
cmd+=("intermediates_dir=getwd())'")

eval ${cmd[*]}


## 07
scriptfile=07-check_NPS_correlation.r
for acohort in OHSU emory
do
    R --slave --args --cohort $acohort < $scriptdir/$scriptfile
done


## 08
scriptfile=08-check_longitudinal_data.rmd
for acohort in OHSU emory
do
    outputfile=2_Pipeline/Reports/explore_longitudinal_$acohort

    cmd=(R --slave -e )
    cmd+=("'rmarkdown::render(")
    cmd+=('"'$scriptdir/$scriptfile'",')
    cmd+=('output_dir="'${output_dir}'",')
    cmd+=('output_file="'${outputfile}'",')
    cmd+=('knit_root_dir=getwd(),')
    cmd+=("intermediates_dir=getwd())'")
    cmd+=('--args --cohort ' $acohort)

    ## execute command 
    eval ${cmd[*]}
done

