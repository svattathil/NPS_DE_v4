library(data.table)
library(argparse)

rm(list = ls())


###### GOAL ######
## Get list of genes that were used in WGCNA in all cohorts


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
parser$add_argument("--residset", type="character", default = "regressSVs", help="")


## Read from parser
args <- parser$parse_args()


### Files that exist
infiles <- sapply(cohorts, function(acohort) {
    paste0("2_Pipeline/04-Build_networks/Wgcna_out/Cohort-specific/",
           acohort, "_", args$residset, "_softpower",
           softpowers[acohort], "_cutheight0.30_gene_conn_stats.txt") },
    simplify = FALSE)


### Files to be created
outfiles <- list(txt = "2_Pipeline/04-Build_networks/shared_genes_for_consensus_module_",
                 args$residset, ".txt")


###### MAIN ######
### Read in data
dats <- sapply(infiles, fread, simplify = FALSE)


### Prepare data
genecounts <- table(unlist(sapply(dats, function(x) { x$protein })))
shared_genes <- names(genecounts)[genecounts == length(cohorts)]


###### FINISH ######
write(shared_genes, file = outfiles[["txt"]], ncol = 1)
