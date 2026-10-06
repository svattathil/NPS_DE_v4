library(data.table)
library(openxlsx)
rm(list=ls())


### Constants ###
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######
### Files that exist
indir <- "2_Pipeline/02-Prepare_analysis_data/"
infiles <- sapply(cohorts, function(cohort) {
    paste0(indir, "phenos_cleaned_SVs_", cohort, ".txt")
    }, simplify = FALSE)



### Files to be created
outdir <- "2_Pipeline/Manuscript/"
outfiles <- list(xl =  paste0(outdir, "batch_correlation_tests.xlsx")
                 )


###### MAIN ######
### Read in data
phenos_list <- sapply(infiles, fread, simplify = FALSE)


### Make objects to hold outputs
## reslist and tablist are tables with one element per cohort
## each of which is a list with one element per nps
## reslist contains the fisher.test output objects
## tablist contains observed and expected counts and standardized residuals
reslist <- tablist <- sapply(cohorts, function(x) {
    alist <- vector(mode = "list", length = length(npsvars.bin))
    names(alist) <- npsvars.bin
    return(alist)
}, simplify = FALSE)


pvals <- data.table(cohort = rep(cohorts,     each =  length(npsvars.bin)),
                    nps    = rep(npsvars.bin, times = length(cohorts)))
pvals[, pval := NA_real_]


### Calculate outputs for each cohort*nps
for(currcohort in cohorts) {
    for(currnps in npsvars.bin) {
        print(paste0(currcohort, "  ", currnps))
        currtab <- table(phenos_list[[currcohort]][, c("Batch", currnps), with = FALSE])

        ## hybrid = TRUE uses chi-squared approximation if
        ##        conditions are met -- all cell counts > 5 and proportions < 0.8
        ## increase workspace as necessary
        fisher.out <- fisher.test(currtab,
                                  hybrid = TRUE,
                                  workspace = 400000)

        ## Get expected counts and standardized residuals using chisq function
        chisq.out <- suppressWarnings(chisq.test(currtab))
        E <- chisq.out$expected
        SR <- chisq.out$stdres

        ## Make tables for output
        obs <- data.table(as.data.frame.matrix(currtab), keep.rownames = "Batch")
        exp <- data.table(as.data.frame.matrix(E), keep.rownames = "Batch")
        standresid <- data.table(as.data.frame.matrix(SR), keep.rownames = "Batch")

        ## make column names
        setnames(obs, c("0", "1"), c("controls.obs", "cases.obs"))
        setnames(exp, c("0", "1"), c("controls.exp", "cases.exp"))
        setnames(standresid, c("0", "1"), c("controls.standardresid", "cases.standardresid"))
        summtab <- Reduce(function(x, y) { merge(x, y, by = "Batch") }, list(obs, exp, standresid))

        ## save to lists
        reslist[[currcohort]][[currnps]] <- fisher.out
        tablist[[currcohort]][[currnps]] <- summtab
        pvals[cohort == currcohort & nps == currnps, pval := fisher.out$p.value]
    }
}


### Add column indicating nominal sig
pvals[pval < 0.05, sig.nominal05 := "sig"]



###### FINISH ######
## Write Excel workbook
## Make one big table per cohort for printing
Write_xlsx_tables(tables = list(Fisher_pvalues = pvals,
                       OHSU = rbindlist(tablist[["OHSU"]], idcol = "nps"),
                       RADC = rbindlist(tablist[["rush"]], idcol = "nps"),
                       emory = rbindlist(tablist[["emory"]], idcol = "nps")),
                  file = outfiles$xl)



