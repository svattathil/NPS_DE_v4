rm(list=ls())
library(knitr)
library(data.table)
library(ggplot2)
options(stringsAsFactors = FALSE)
source("~/comsv/Svattathil_Library/svattathil_functions.r")


### CONSTANTS ###
source("1_Code/project_constants.r")



### SETUP ###
### Define input
indir1 <- "2_Pipeline/01-Explore_data/"
indir2 <- "2_Pipeline//02-Prepare_analysis_data/"
files_patho <- sapply(cohorts, function(x) { paste0(indir1, "patho_", x, ".txt") })
files_main <- sapply(cohorts, function(x) { paste0(indir2, "phenos_cleaned_SVs_", x, ".txt") })

### Define output
outfiles <- list(
    txt = paste0(indir2, "AD_pathologic_diagnosis_allcohorts.txt")
    )

### MAIN ###
### Read data
main <- sapply(files_main, fread)
patho <- sapply(files_patho, fread)

### Prepare data
## Harmonize ID columns
setnames(patho$OHSU, "bfonum", "ID")
setnames(patho$rush, "projid", "ID")
main$emory[, ID := NULL]
setnames(main$emory, "GUID", "ID")
setnames(patho$emory, "GUID", "ID")



## Merge main phenos and new patho/severity phenos
merged <- sapply(cohorts, function(acohort) {
    merge(main[[acohort]], patho[[acohort]], by = "ID", all.x = TRUE)
})

## For OHSU, prepare the path_dx variable to a more readable form
merged$OHSU[, path_dx_merged := gsub("1.0\\)", "", path_dx_merged)]
merged$OHSU[, path_dx_merged := gsub("1\\)", "", path_dx_merged)]
merged$OHSU[, path_dx_merged := gsub("2.0\\)", "+", path_dx_merged)]
merged$OHSU[, path_dx_merged := gsub("2.0", "+", path_dx_merged)]
merged$OHSU[, path_dx_merged := gsub("3.0\\)", "+", path_dx_merged)]
merged$OHSU[, path_dx_merged := gsub("2\\)", "+", path_dx_merged)]
merged$OHSU[, path_dx_merged := gsub("3\\)", "+", path_dx_merged)]
merged$OHSU[, path_dx_merged := gsub("  ", " ", path_dx_merged)]
merged$OHSU[, path_dx_merged := trimws(path_dx_merged, which = "left")]


### Initiate list to hold table of AD status per cohort
adstatus_list <- list()

### Get OHSU status- parse the pathologic diagnosis column to find those with primary
### diagnosis of AD or AP (Alzheimer's pathology)
adstatus_list[["OHSU"]] <- merged$OHSU[, .(ID, protsample, path_dx_merged)]
adstatus_list[["OHSU"]][, AD_patho_diagnosis := 0]
adstatus_list[["OHSU"]][grep("^AD|^AP", path_dx_merged), AD_patho_diagnosis := 1]
adstatus_list[["OHSU"]][, path_dx_merged := NULL]


### Get Rush status - use binary status based on Reagan criteria
adstatus_list[["rush"]] <- merged$rush[, .(ID, protsample, ad_reagan)]
setnames(adstatus_list[["rush"]], "ad_reagan", "AD_patho_diagnosis")


### Get Emory status -  all participants have neuropathologic diagnosis of AD
adstatus_list[["emory"]] <- data.table(
                     ID = merged$emory[, ID],
                     protsample = merged$emory[, protsample],
                     AD_patho_diagnosis = 1)


### Flatten to one table across cohorts
adstatus <- rbindlist(adstatus_list, idcol = "cohort")


### FINISH ###
## Write to file
write.table(adstatus, file = outfiles$txt, row.names = FALSE, quote = FALSE, sep = "\t")
