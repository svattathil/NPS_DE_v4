library(data.table)
library(argparse)
library(metafor)
library(qvalue)
rm(list = ls())



###### GOAL ######

##### CONSTANTS #####
options(stringsAsFactors = FALSE)
cohorts <- c("OHSU", "rush", "emory")


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


##### COMMAND LINE ARGUMENTS #####
## Create parser
parser <- ArgumentParser()
parser$add_argument("--auto_softpower", type="integer", default=9, help="")

## Read from parser
args <- parser$parse_args()


###### SETUP ######

### Files that exist
indir <- paste0("2_Pipeline/04-Build_networks/Wgcna_out/Consensus/Auto_softpower", args$auto_softpower, "/")
infiles <- list(res = paste0(indir, "nps_regression_results_per-cohort.txt"))


### Files to be created
outfiles <- list(fe_results = paste0(indir, "nps_regression_results_meta.txt"))


###### MAIN ######
### Read in data
res <- fread(infiles[["res"]])


### Prepare data
setnames(res, "Std.Err", "SE")
res[, module := gsub("ME", "", eigengene)]

## Keep only the columns needed for metafor
meta_long <- res[, .(module, NPS, cohort, yi = Estimate, sei = SE)]

## Run meta-analysis for each module-NPS pair and
## save results in a data table
fe_results <- meta_long[
  , {
      fit <- rma.uni(yi = yi, sei = sei, method = "FE")  # FE = fixed-effects
      list(
          fe_beta   = coef(fit),                         # pooled estimate
          fe_se     = sqrt(vcov(fit)),                   # pooled SE
          fe_z      = summary(fit)$zval,                 # Wald Z
          fe_p      = summary(fit)$pval,                 # two-sided p-value
          fe_Q      = fit$QE,                            # heterogeneity Q (should be low),
          fe_Q_pval = pchisq(fit$QE, df = length(cohorts)-1, lower.tail = FALSE), # p-value of Cochran's Q
          fe_I2     = fit$I2                             # I² (should be ≈0 for FE)
      )
  },
  by = c("module", "NPS")]


### Add multiple testing correction columns
setorder(fe_results, fe_p)
fe_results[, qvalue := qvalue(fe_p)$qvalue]
fe_results[, p.BH := p.adjust(fe_p, method = "BH")]
fe_results[, p.bonf := p.adjust(fe_p, method = "bonferroni")]
fe_results[, sig.qvalue05 := ifelse(qvalue < 0.05, TRUE, FALSE)]


###### FINISH ######
write.table(fe_results, file = outfiles[["fe_results"]], row.names = FALSE, quote = FALSE, sep = "\t")
