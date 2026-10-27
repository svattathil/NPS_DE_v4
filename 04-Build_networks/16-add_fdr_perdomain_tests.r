library(data.table)
library(argparse)
library(qvalue)
rm(list = ls())


###### GOAL ######



##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")


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
fileid <- "auto_softpower9"
indir <- paste0("2_Pipeline/04-Build_networks/Using_ignoreSVs/",
                "Consensus/Auto_softpower9/")
infiles <- list(cohort = paste0(indir, fileid, "_nps_linear_regression_results_per-cohort.txt"),
                meta   = paste0(indir, fileid, "_nps_regression_results_meta.txt")
                )


### Files to be created
outfiles <- list(xlsx = paste0(indir, fileid, "_perdomain_association_results_fdr.xlsx")
                 )


###### MAIN ######
### Read in data
res_list <- sapply(infiles, fread, simplify = FALSE)


### Prepare data
## Split results by cohort
res_bycohort <- sapply(cohorts, function(acohort) { res_list$cohort[cohort == acohort, ] },
                       simplify = FALSE)

### add qvalue for cohort-specific results
## calculate within cohort (across n_modules x n_domains tests)
for(acohort in cohorts) {
    res_bycohort[[acohort]][, p.BH := p.adjust(Pvalue, method = "BH")]
    res_bycohort[[acohort]][, qvalue := qvalue(Pvalue)$qvalues]
    res_bycohort[[acohort]][qvalue < 0.05, sig_qvalue05 := "sig"]

    setorder(res_bycohort[[acohort]], Pvalue)
}
rm(acohort)


### The meta-analysis already have BH and qvalues
## across n_modules x n_domains tests



###### FINISH ######
forprint <- res_bycohort
forprint[["meta_analysis"]] <- res_list$meta

Write_xlsx_tables(tables = forprint, file = outfiles$xlsx)
