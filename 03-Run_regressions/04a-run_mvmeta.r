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
source("~/comsv/Svattathil_Library/manhattan_and_qq_plot_functions.r")


###### SETUP ######
## Create parser
parser <- ArgumentParser()
parser$add_argument("--residset", type="character",
                    help="regressSVs or ignoreSVs or resid2_multi or noProtect",
                    default = "ignoreSVs")

## Read from parser
args <- parser$parse_args()



### Files that exist
indir1 <- paste0("2_Pipeline/03-Run_regressions/Prot_multivar/Using_",
                 Capwords(args$residset), "/Mvmeta/")
indir2 <- paste0(indir1, "In/")
File.betas <- function(acohort) {
    paste0(indir2, "prot_multivar_", args$residset, "_", acohort, "_betas.rds") }

File.covarmat <- function(acohort) {
    paste0(indir2, "prot_multivar_", args$residset, "_", acohort, "_covariances.rds") }


### Files to be created
outdir <- paste0(indir1, "Out/")
MyMkdir(outdir)
outfiles <- list(tab = paste0(outdir, "wald_results.txt"),
                 qq  = paste0(outdir, "qqplot.png")
                 )


###### MAIN ######
### Read in data
betas <- sapply(cohorts, function(x) { readRDS(File.betas(x)) }, simplify = FALSE)
covarmats <- sapply(cohorts, function(x) { readRDS(File.covarmat(x)) }, simplify = FALSE)


### Prepare data
## Make vector of proteins tested in all three cohorts
proteins_to_test <- intersect(names(betas[[1]]), intersect(names(betas[[2]]), names(betas[[3]])))


## Run meta-analysis for each protein in turn
meta_res <- vector("numeric", length(proteins_to_test))
contributing_cohorts <- vector("character", length(proteins_to_test))


for(i in 1:length(proteins_to_test)) {
    if(i %% 1000 == 0) { print(paste0(pN(i), " of ", pN(length(proteins_to_test)))) }

    aprot <- proteins_to_test[i]

    ## Extract beta vector and covariance matrix for current protein for each cohort
    ## covarcols define column/row names corresponding to NPS of interest
    ## beta_mat is a matrix with one row per cohort and one column per NPS
    ## cov_list is a list with one element per cohort, each element is a n.nps x n.nps matrix
    covarcols <- paste0("nps_id", vars_formulti, ":protein")
    beta_mat <- t(sapply(cohorts, function(acohort) {
        betas[[acohort]][[aprot]][vars_formulti] }))  ## 3 x length(vars_formulti)

    cov_list  <- sapply(cohorts, function(acohort) {
        currmat <- covarmats[[acohort]][[aprot]]

        ## some matrices are NULL for whatever reason; handle that
        if(length(currmat) == 1 && is.na(currmat)) {
            return(NA)
        }
        return(currmat[covarcols, covarcols])
    }, simplify = FALSE)

    ## Which cohorts have a usable (non-NA) result for this protein?
    ok_cohorts <- cohorts[!sapply(cov_list, function(x) length(x) == 1 && is.na(x))]

    if (length(ok_cohorts) >= 2) {   # need at least 2 cohorts to meta-analyze
        ## Fit the multivariate meta‑analysis
        ## To get pooled protein effect
        fit_mv <- mixmeta(
            beta_mat[ok_cohorts, , drop = FALSE],
            S = cov_list[ok_cohorts],
            method = "fixed")

        ## Joint test that all 9 p-values are 0
        res <- car::linearHypothesis(
                        fit_mv,
                        hypothesis.matrix = diag(length(coef(fit_mv))))

        protres <- res[["Pr(>Chisq)"]][2]
    }else {
        protres <- NA_real_
    }

    ## Save results to list
    meta_res[i]             <- protres
    contributing_cohorts[i] <- paste0(ok_cohorts, collapse=",")
}
#rm(i)


### Make table for print
summstats <- data.table(protein = proteins_to_test, p = meta_res, cohorts = contributing_cohorts)
setorder(summstats, p)
summstats[, qvalue := qvalue(p)$qvalues]
summstats[qvalue < 0.05, sig.qvalue05 := "sig"]


### Draw qq plot
n.na <- sum(is.na(summstats$p))
titleqq <- paste0("residset: ", args$residset, "\n", "proteins tested = ", length(proteins_to_test))
if(n.na > 0) { titleqq <- paste0(titleqq, "; excluding ", n.na, " proteins with NA P-value") }
qqplot <- qqunif.plot(summstats[!is.na(p), p], maintitle=titleqq)


###### FINISH ######
write.table(summstats, file = outfiles[["tab"]], row.names = FALSE, quote = FALSE, sep = "\t")

### Draw qqplot
png(outfiles$qq, height=7, width=7, units="in", res=300)
print(qqplot)
dev.off()
