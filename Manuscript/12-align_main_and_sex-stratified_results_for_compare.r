library(data.table)
library(argparse)
library(tidyr)

rm(list = ls())



###### GOAL ######



##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")
npsmap <- data.table(nps = npsvars.bin, nps_name = names(npsvars.bin))
setkey(npsmap, nps)

sets <- factor(c("main", "males", "females"), ordered = TRUE)

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
indir <- "/nfs/wingofss2/data/commons/svattat/Projects/NPS_protein_DE/2_Pipeline/03-Run_regressions/Prot_perDomain/"

files1 <- sapply(cohorts, function(acohort) {
  paste0(indir, "Basic_withSV/basic_withSV_", acohort, ".tsv")
})
names(files1) <- paste0("main_", cohorts)
files1[["main_meta"]] <- paste0(indir, "Basic_withSV/Metafor/metafor_basic_withSV_using_SE_mod.txt")


files2 <- sapply(cohorts, function(acohort) {
  paste0(indir, "Males_withSV/males_withSV_", acohort, ".tsv")
})
names(files2) <- paste0("males_", cohorts)
files2[["males_meta"]] <- paste0(indir, "Males_withSV/Metafor/metafor_males_withSV_using_SE_mod.txt")


files3 <- sapply(cohorts, function(acohort) {
  paste0(indir, "Females_withSV/females_withSV_", acohort, ".tsv")
})
names(files3) <- paste0("females_", cohorts)
files3[["females_meta"]] <- paste0(indir, "Females_withSV/Metafor/metafor_females_withSV_using_SE_mod.txt")

infiles <- unlist(list(files1, files2, files3), recursive = FALSE)


### Files to be created
outfiles <- list(
    xlsx = paste0("2_Pipeline/Manuscript/results_aligned_main_males_females_forcompare.xlsx")
)


###### MAIN ######
### Read in data
dat_list <- sapply(infiles, fread, simplify = FALSE)


### Set up list to hold sheets
sheetlist <- vector(mode = "list", length = (length(cohorts)+1)*length(npsvars.bin))
names(sheetlist) <- paste0(rep(c(cohorts, "meta"), each = length(npsvars.bin)),
                             "_",
                             rep(names(npsvars.bin), times = length(cohorts)+1))

### For each cohort, for each BPSD, gather results from main, males, and females analyses
### One table per cohort*BPSD
vals <- c("logFC", "CI.L", "CI.R", "t", "P.Value", "fdr") ## stats columns to keep from each results set

for(cohort in cohorts) {
    currdat_list <- sapply(sets, function(x) {
        adat <- copy(dat_list[[paste0(x, "_", cohort)]]) ## data for this cohort*set
        return(adat)
    }, simplify = FALSE)
    names(currdat_list) <- sets
    currdat <- rbindlist(currdat_list, idcol = "set")

    for(currnps in npsvars.bin) {
        ## Make wide table
        dat_nps <- currdat[domain == currnps, ] ## extract rows for this NPS domain
        wide_dat <- dcast(dat_nps, domain + protein ~ set, value.var = vals, sep = "_")
        setcolorder(wide_dat, c("domain", "protein", as.vector(outer(vals, sets, paste, sep = "_"))))
        setorder(wide_dat, P.Value_main, na.last = TRUE)

        ## Save to sheetlist
        datlabel <- paste0(cohort, "_", npsmap[currnps, nps_name])
        sheetlist[[datlabel]] <- wide_dat
        rm(list = c("datlabel", "wide_dat"))
    }
    rm(list = c("cohort", "currnps", "currdat_list", "currdat"))
}


### Do similar for the meta-analysis results
### One table per NPS
vals_meta <- c("fe_beta", "fe_se", "fe_p", "qvalue", "p.BH")

currdat_list <- sapply(sets, function(x) {
        adat <- copy(dat_list[[paste0(x, "_meta")]]) ## data for this cohort*set
        return(adat)
}, simplify = FALSE)
names(currdat_list) <- sets
currdat <- rbindlist(currdat_list, idcol = "set")

for(currnps in npsvars.bin) {
    ## Make wide table
    dat_nps <- currdat[nps == currnps, ] ## extract rows for this NPS domain
    wide_dat <- dcast(dat_nps, nps + protein ~ set, value.var = vals_meta, sep = "_")
    setcolorder(wide_dat, c("nps", "protein", as.vector(outer(vals_meta, sets, paste, sep = "_"))))
    setorder(wide_dat, fe_p_main, na.last = TRUE)

    ## Save to sheetlist
    datlabel <- paste0("meta_", npsmap[currnps, nps_name])
    sheetlist[[datlabel]] <- wide_dat
    rm(list = c("datlabel", "wide_dat"))
}



###### FINISH ######
Write_xlsx_tables(tables = sheetlist, file = outfiles$xlsx)
