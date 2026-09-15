library(data.table)
library(argparse)
library(readr)
library(openxlsx)

rm(list = ls())


###### GOAL ######



### CONSTANTS ###
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")

PrintValidValues <- function(an_arg) {
    paste0(valid_values[[an_arg]], collapse = ", ")
}


###### SETUP ######
### Command line arguments

## Create parser
parser <- ArgumentParser()

## Active arguments
parser$add_argument("--run",    type="character", default="basic")

## Read from parser
args <- parser$parse_args()


### Set some values based on the run specifics
if(args$run == "severity") {
    npsvars <- npsvars.sev
    cohorts <- c("OHSU", "emory")
} else {
    npsvars <- npsvars.bin
}


for(acohort in cohorts) {
    ### Files that exist
    indir <- paste0("2_Pipeline/03-Run_regressions/Prot_perDomain/", Capwords(args$run), "/")
    fileid <- paste0(acohort, "_", args$run)

    infiles <- list(
        limma_per_domain = paste0(indir, fileid, ".tsv"),
        summary = paste0(indir, fileid, "_summary.tsv")
    )


    ### Files to be created
    output_file <- paste0(indir, paste(c("stats", args$run, acohort), collapse="_"), ".xlsx")


    ### Read in data
    res <- fread(infiles$limma_per_domain)
    summ_contents <- fread(infiles$summary)


    ### Split results by domain
    nps_list <- vector("list", length(npsvars))
    names(nps_list) <- npsvars
    for(adomain in npsvars) { nps_list[[adomain]] <- res[domain == adomain, ] }
    names(nps_list) <- names(npsvars)


    ### Write to Excel workbook –
    ## first sheet is logs,
    ## then one worksheet per NPS
    wb <- createWorkbook()

    addWorksheet(wb, sheetName = "summary")
    writeData(wb, sheet = "summary", summ_contents)

    for (sheet_name in names(nps_list)) {
        addWorksheet(wb, sheetName = sheet_name)
        writeData(wb, sheet = sheet_name, nps_list[[sheet_name]])
    }
    saveWorkbook(wb, file = output_file, overwrite = TRUE)

    cat("✅ Workbook saved as:", output_file, "\n")
}

