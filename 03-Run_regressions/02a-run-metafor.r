library(data.table)
library(argparse)
library(metafor)
library(qvalue)

rm(list = ls())


###### GOAL ######
## Run cross-cohort meta-analysis for specified outcome


##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######
### Command line arguments
## Create parser
parser <- ArgumentParser()
parser$add_argument("--run", type="character", default="basic", help = "basic or cond_ind")
parser$add_argument("--SE_type", type="character", default="SE_mod", help = "SE_mod or SE_unmod")


## Read from parser
args <- parser$parse_args()

## Define outcomes and subset of cohorts (if necessary)
if(args$run %in% c("basic", "cond_ind",  "adj_e4_exp", "adj_e4_count", "stratify_e4_any")) {
    npsvars <- npsvars.bin
}

if(args$run %in% c("severity")) {
    npsvars <- npsvars.sev
    cohorts <- c("OHSU", "emory")
}


### Files that exist
resdir <- paste0("2_Pipeline/03-Run_regressions/Prot_perDomain/", Capwords(args$run), "/")

if(args$run %in% c("basic", "severity")) {
    Statsfile <- function(acohort) {  paste0(resdir, acohort, "_", args$run, ".tsv") }
}
if(args$run == "cond_ind") {
    Statsfile <- function(acohort) {  paste0(resdir, acohort, "_joint_perDomain.tsv") }
}


### Files to be created
outdir <- paste0(resdir, "Metafor/")
MyMkdir(outdir)
out.stats <- paste0(outdir, "metafor_", args$run, "_using_", args$SE_type, ".txt")


###### MAIN ######
### Read in data
res.list.cohort <- sapply(cohorts, function(x) { fread(Statsfile(x)) }, simplify = FALSE)


### Prepare data
## Split by domain
res.list.nps <- sapply(as.character(npsvars), function(anps) {
    sapply(res.list.cohort, function(cohortdat) {
        cohortdat[domain == anps, ]
        }, simplify = FALSE)
    }, simplify = FALSE)


### Initialize list to hold results across NPS
fe_results_list.nps <- list()


### Run meta-analysis for each NPS
for(currnps in npsvars) {
    cat("Running for ", currnps, "\n")
    ### Make long table for current NPS
    res.list <- res.list.nps[[currnps]] # list with one element per cohort
    meta_long <- rbindlist(res.list, use.names = TRUE, fill = TRUE, idcol = "cohort")

    ## Keep only the columns needed for metafor
    meta_long <- meta_long[, .(protein, cohort, yi = logFC, sei = get(args$SE_type))]

    ## Filter to only the proteins that are present for all three cohorts
    proteincounts <- meta_long[, .N, by = "protein"]
    meta_long <- meta_long[protein %in% proteincounts[N == length(cohorts), protein], ]


    ### Run metafor and make table of stats
    fe_results <- meta_long[
      , {
          fit <- rma.uni(yi = yi, sei = sei, method = "FE")   # FE = fixed‑effects
          list(
              fe_beta   = coef(fit),                # pooled estimate
              fe_se     = sqrt(vcov(fit)),          # pooled SE
              fe_z      = summary(fit)$zval,        # Wald Z
              fe_p      = summary(fit)$pval,        # two‑sided p‑value
              fe_Q      = fit$QE,                   # heterogeneity Q (should be low),
              fe_Q_pval = pchisq(fit$QE, df = length(cohorts)-1, lower.tail = FALSE), # p-value of Cochran's Q
              fe_I2     = fit$I2                    # I² (should be ≈0 for FE)
          )
      },
      by = protein]


    ### Add multiple testing correction columns
    setorder(fe_results, fe_p)
    fe_results[, qvalue := qvalue(fe_p)$qvalue]
    fe_results[, p.BH := p.adjust(fe_p, method = "BH")]
    fe_results[, p.bonf := p.adjust(fe_p, method = "bonferroni")]
    fe_results[, sig.qvalue05 := ifelse(qvalue < 0.05, "sig", "")]


    ### Save result to cross-nps list
    fe_results_list.nps[[currnps]] <- fe_results
    rm(fe_results) ## clean up
}
rm(currnps)


### Make flat table for printing
## Can reshape as desired in a later step for supplementary table
## or followup analyses
fe_results_allnps <- rbindlist(fe_results_list.nps, idcol = "nps")


###### FINISH ######
fwrite(fe_results_allnps, file = out.stats, quote = FALSE, sep = "\t", na = NA)

