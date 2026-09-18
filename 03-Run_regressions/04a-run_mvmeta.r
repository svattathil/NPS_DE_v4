library(data.table)
library(argparse)
library(mixmeta)
library(qvalue)
rm(list = ls())


###### GOAL ######
## Run meta-analysis on the multivariate regression results


### Constants
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")
vars_formulti <- npsvars.bin


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######
### Files that exist
indir <- "2_Pipeline/03-Run_regressions/Prot_multivar/Mvmeta/In/"
File.betas <- function(acohort) {
    paste0(indir, "prot_multivar_", acohort, "_betas.rds") }

File.covarmat <- function(acohort) {
    paste0(indir, "prot_multivar_", acohort, "_covariances.rds") }


### Files to be created
outdir <- "2_Pipeline/03-Run_regressions/Prot_multivar/Mvmeta/Out/"
MyMkdir(outdir)
outfiles <- list(out_tab = paste0(outdir, "wald_results.txt"))


###### MAIN ######
### Read in data
betas <- sapply(cohorts, function(x) { readRDS(File.betas(x)) }, simplify = FALSE)
covarmats <- sapply(cohorts, function(x) { readRDS(File.covarmat(x)) }, simplify = FALSE)


### Prepare data
## Make vector of proteins tested in all three cohorts
proteins_to_test <- intersect(names(betas[[1]]), intersect(names(betas[[2]]), names(betas[[3]])))

## Run meta-analysis for each protein in turn
set.seed(98873)
meta_res <- sapply(1:length(proteins_to_test), function(i) {
    if(i %% 1000 == 0) { print(paste0(pN(i), " of ", pN(length(proteins_to_test)))) }

    aprot <- proteins_to_test[i]

    ## Extract beta vector and covariance matrix for current protein for each cohort
    ## covarscols define column/row names corresponding to NPS of interest
    ## beta_mat is a matrix with one row per cohort and one column per NPS
    ## cov_list is a list with one element per cohort, each element is a n.nps x n.nps matrix
    covarcols <- paste0("protein:nps_id", vars_formulti)
    beta_mat <- t(sapply(cohorts, function(acohort) {
        betas[[acohort]][[aprot]][vars_formulti] }))  ## 3 x length(vars_formulti)
    cov_list  <- lapply(cohorts, function(acohort) {
        covarmats[[acohort]][[aprot]][covarcols, covarcols] })


    ## Fit the multivariate random‑effects meta‑analysis
    ## This gives pooled protein effect
    fit_mv <- mixmeta(beta_mat,
                      S = cov_list,
                      method = "reml")
                                        #summary(fit_mv)

    ## Joint test that all 9 p-values are 0
    res <- car::linearHypothesis(fit_mv, hypothesis.matrix = diag(length(coef(fit_mv))))

    return(res[["Pr(>Chisq)"]][2])
})


### Make object for print
forprint <- data.table(protein = proteins_to_test, p = meta_res)
setorder(forprint, p)
forprint[, qvalue := qvalue(p)$qvalues]


###### FINISH ######
write.table(forprint, file = outfiles[["out_tab"]], row.names = FALSE, quote = FALSE, sep = "\t")
