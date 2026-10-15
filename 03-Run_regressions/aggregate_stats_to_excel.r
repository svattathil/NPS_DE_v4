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
parser$add_argument("--analysis", default = "basic",
                    help = paste0("Specify the analysis (e.g., basic, severity, ",
                                  "males, females, e4_adj, latfactor)"))
parser$add_argument("--svchoice", default = "noSV",
                    help = paste0("withSV or noSV"))
## Read from parser
args <- parser$parse_args()

## Define runid
runid <- paste0(args$analysis, "_", args$svchoice)


### Set some values based on the run specifics
if(args$analysis == "severity") {
    npsvars <- npsvars.sev
    cohorts <- c("OHSU", "emory")
} else {
    npsvars <- npsvars.bin
}


for(acohort in cohorts) {
    ### Files that exist
    if(args$analysis == "cond_ind") {
        ## the cond_ind analysis is done inside the Basic analysis
        indir <- paste0("2_Pipeline/03-Run_regressions/Prot_perDomain/", paste0("Basic_", args$svchoice), "/")
    }
    else {
        indir <- paste0("2_Pipeline/03-Run_regressions/Prot_perDomain/", Capwords(runid), "/")
    }

    infiles <- list(
        limma_per_domain = paste0(indir, runid, "_", acohort, ".tsv"),
        summary          = paste0(indir, runid, "_summary_", acohort, ".tsv")
    )


    ### Files to be created
    output_file <- paste0(indir, paste(c("stats", runid, acohort), collapse="_"), ".xlsx")


    ### Read in data
    res <- fread(infiles$limma_per_domain)
    if(args$analysis != "cond_ind") {
        summ_contents <- fread(infiles$summary)
    }

    ### Split results by domain
    nps_list <- vector("list", length(npsvars))
    names(nps_list) <- npsvars
    for(adomain in npsvars) { nps_list[[adomain]] <- res[domain == adomain, ] }
    names(nps_list) <- names(npsvars)


    ### Write to Excel workbook –
    ## first sheet is logs,
    ## then one worksheet per NPS
    wb <- createWorkbook()

    if(args$analysis != "cond_ind") {
        addWorksheet(wb, sheetName = "summary")
        writeData(wb, sheet = "summary", summ_contents)
    }

    for (sheet_name in names(nps_list)) {
        addWorksheet(wb, sheetName = sheet_name)
        writeData(wb, sheet = sheet_name, nps_list[[sheet_name]])
    }
    saveWorkbook(wb, file = output_file, overwrite = TRUE)

    cat("✅ Workbook saved as:", output_file, "\n")
}

