library(data.table)
library(argparse)

rm(list = ls())



###### GOAL ######
## pull together association statistics to include in supplementary tables


##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")

npsmap <- data.table(code_name = npsvars.bin, text_name = names(npsvars.bin))


##### COMMAND LINE ARGUMENTS #####
## Create parser
parser <- ArgumentParser()
parser$add_argument("--auto_softpower", type="integer", default=9, help="")
parser$add_argument("--residset", type="character", default="regressSVs", help="")

## Read from parser
args <- parser$parse_args()



###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######
### Files that exist
### module eigengene ~ NPS association statistics

dir1 <- paste0("2_Pipeline/04-Build_networks/Wgcna_out/Consensus/",
                stringr::str_to_title(paste0(args$residset, "_Auto_softpower",
                                             args$auto_softpower)), "/")

dir_npsassoc <- dir1
fileid <- paste0(args$residset, "_auto_softpower", args$auto_softpower)
infiles_npsassoc <- list(per_cohort = paste0(dir1, fileid, "_nps_regression_results_per-cohort.txt"),
                        meta       = paste0(dir1, fileid, "_nps_regression_results_meta.txt"))


### module eigengene ~ latent factor association statistics
dir_lfassoc <- paste0(dir1, "DEA_latent_factors/")
infiles_lfassoc <- sapply(cohorts, function(acohort) {
    paste0(dir_lfassoc,
           fileid, "_", acohort, "_eigengene_regress_stats_latent_factors.txt")
}, simplify = FALSE)


### module eigengene ~ multivariate association statistics
dir_multassoc <- paste0(dir1, "Multivariate/")
infiles_multassoc <- sapply(cohorts, function(acohort) {
    paste0(dir_multassoc, "eigengene_multivar_", fileid, "_", acohort, "_regression_stats.txt")
}, simplify = FALSE)


### Files to be created
outpref <- paste0("aggregated_association_results_for_supp_", fileid)
outfiles <- list(nps     = paste0(dir1, outpref, "_nps.txt"),
                 lf_mult = paste0(dir1, outpref, "_latfac_multivariate.txt"))


###### MAIN ######
### Read in data
## single file with all per-cohort results
res_npsassoc_percohort <- fread(infiles_npsassoc[["per_cohort"]])

## single file with meta-analysis results
res_npsassoc_meta <- fread(infiles_npsassoc[["meta"]])

## one file per cohort
reslist_lfassoc <- sapply(cohorts, function(acohort) {
    fread(infiles_lfassoc[[acohort]])
}, simplify = FALSE)

## one file per cohort
reslist_multassoc <- sapply(cohorts, function(acohort) {
    fread(infiles_multassoc[[acohort]])
}, simplify = FALSE)



#### 1. Prepare table for individual NPS associations ####
####    One row per module-NPS                        ####

### Initiate list to hold objects to merge
formerge <- vector(mode = "list")

### Prepare per-cohort module ~ NPS association statistics
## harmonize module name
res_npsassoc_percohort[, module_color := sapply(eigengene, function(x) { gsub("ME", "", x) })]

## extract columns of interest
setnames(res_npsassoc_percohort, c("Std.Err"), c("SE"))
idcols <- c("module_color", "NPS")
namecols <- c("cohort")
valcols <- c("Estimate", "SE", "Pvalue")
tmp <- res_npsassoc_percohort[, c(idcols, namecols, valcols), with = FALSE]

## Reformat to wide table
## and group columns by cohort
formerge[["nps_assoc_percohort"]] <- data.table::dcast(tmp, module_color + NPS ~ cohort,
                                          value.var = valcols,
                                          sep = "_")

ordered_valcols <- paste0(rep(valcols, times = length(cohorts)), "_",
                          rep(cohorts, each = length(valcols)))
setcolorder(formerge[["nps_assoc_percohort"]], c(idcols, ordered_valcols))


### Prepare meta-analysis statistics for module ~ NPS association
setnames(res_npsassoc_meta, c("module", "fe_beta", "fe_se", "fe_p"),
         c("module_color", "Estimate", "SE", "Pvalue"), skip_absent = TRUE)

formerge[["nps_assoc_meta"]] <- res_npsassoc_meta[, c(idcols, valcols, "fe_I2"), with = FALSE]

### Merge into final table
forprint_nps <- merge(formerge[["nps_assoc_percohort"]], formerge[["nps_assoc_meta"]],
                      by = idcols, all = TRUE)
forprint_nps <- merge(forprint_nps, npsmap, by.x = "NPS", by.y = "code_name")
forprint_nps[, NPS := text_name]
forprint_nps[, text_name := NULL]


#### 2. Prepare table for latent factor and multivariate associations ####
####    One row per module                                            ####
### Initiate list to hold objects to merge
formerge <- vector(mode = "list")


### Prepare module ~ latent factor association statistics
## Flatten into one table for all cohorts
res_lfassoc <- rbindlist(reslist_lfassoc, idcol = "cohort", fill = TRUE)

## harmonize module name
res_lfassoc[, module := sapply(module, function(x) { gsub("ME", "", x) })]

## extract columns of interest
idcols <- c("module")
namecols <- c("cohort")
valcols <- paste0(rep(c("f1"), each = 2), "_", rep(c("coeff", "pval"), 1))
tmp <- res_lfassoc[, c(idcols, namecols, valcols),  with = FALSE]


## Reformat to wide table
## and group columns by cohort
formerge[["lf_assoc"]] <- data.table::dcast(tmp, module ~ cohort,
                                          value.var = valcols,
                                          sep = "_")
ordered_valcols <- paste0(rep(valcols, times = length(cohorts)), "_",
                          rep(cohorts, each = length(valcols)))
setcolorder(formerge[["lf_assoc"]], c(idcols, ordered_valcols))


### Prepare module ~ multivariate association statistics
## Flatten into one table for all cohorts
res_multassoc <- rbindlist(reslist_multassoc, idcol = "cohort", fill = TRUE)

## harmonize module name
res_multassoc[, module := sapply(eigengene, function(x) { gsub("ME", "", x) })]

## extract columns of interest
valcols <- c("Fstat", "P")
tmp <- res_multassoc[, c(idcols, namecols, valcols), with = FALSE]

## Reformat to wide table
## and group columns by cohort
formerge[["mult_assoc"]] <- data.table::dcast(tmp, module ~ cohort,
                                              value.var = valcols,
                                              sep = "_")
ordered_valcols <- paste0(rep(valcols, times = length(cohorts)), "_",
                          rep(cohorts, each = length(valcols)))
setcolorder(formerge[["mult_assoc"]], c(idcols, ordered_valcols))


### Merge the tables
forprint_lf_mult <- merge(formerge[["lf_assoc"]], formerge[["mult_assoc"]],
                          by = "module", all = TRUE)



###### FINISH ######
write.table(forprint_nps, file = outfiles[["nps"]],
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(forprint_lf_mult, file = outfiles[["lf_mult"]],
            row.names = FALSE, quote = FALSE, sep = "\t")
