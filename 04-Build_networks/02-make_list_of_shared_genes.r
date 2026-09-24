library(data.table)
library(argparse)

rm(list = ls())



###### GOAL ######



##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")
softpowers <- c(OHSU = 9, rush = 9, emory = 9)


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######
### Command line arguments
## Create parser
parser <- ArgumentParser()
#parser$add_argument("--", type="character", help="")
#parser$add_argument("--", action="store_true", default=TRUE, help="")
#parser$add_argument("--", action="store_true", default=FALSE, help="")
#parser$add_argument("--", action="store_false", dest="", help="")
#parser$add_argument("--", type="integer", default=10, help="")
#parser$add_argument("--", type="double", default=10, help="")

## Read from parser
args <- parser$parse_args()


### Files that exist
infiles <- sapply(cohorts, function(acohort) {
    paste0("2_Pipeline/04-Build_networks/Wgcna_out/Cohort-specific/",
           acohort, "_ignoreSVs_softpower",
           softpowers[acohort], "_cutheight0.30_gene_conn_stats.txt") },
    simplify = FALSE)


### Files to be created
outfiles <- list(txt = "2_Pipeline/04-Build_networks/shared_genes_for_consensus_module.txt")


###### MAIN ######
### Read in data
dats <- sapply(infiles, fread, simplify = FALSE)


### Prepare data
genecounts <- table(unlist(sapply(dats, function(x) { x$protein })))
shared_genes <- names(genecounts)[genecounts == length(cohorts)]


###### FINISH ######
write(shared_genes, file = outfiles[["txt"]], ncol = 1)
