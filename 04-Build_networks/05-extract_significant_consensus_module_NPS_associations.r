library(data.table)
library(argparse)

rm(list = ls())



###### GOAL ######


##### CONSTANTS #####
options(stringsAsFactors = FALSE)



###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")



###### SETUP ######
### Command line arguments
## Create parser
parser <- ArgumentParser()
parser$add_argument("--residset", type="character", default="regressSVs", help="")
parser$add_argument("--auto_softpower", type = "integer", help = "", default = 9)

## Read from parser
args <- parser$parse_args()


### Files that exist
indir <- paste0("2_Pipeline/04-Build_networks/Wgcna_out/Consensus/",
                stringr::str_to_title(paste0(args$residset, "_Auto_softpower",
                                             args$auto_softpower)), "/")

fileid <- paste0(args$residset, "_auto_softpower", args$auto_softpower)
infiles <- list(specific = paste0(indir, fileid, "_nps_regression_results_per-cohort.txt"),
                meta     = paste0(indir, fileid, "_nps_regression_results_meta.txt"))


### Files to be created
outfiles <- list(sig = paste0(indir, "nps_regression_sig_", fileid, ".txt"))


###### MAIN ######
### Read in data
dats <- sapply(infiles, fread, simplify = FALSE)


### Prepare data
## Add module column to per-cohort results
dats[["specific"]][, module := gsub("ME", "", eigengene)]
setnames(dats[["specific"]], "Std.Err", "SE")


## Add cohort column to meta results
dats[["meta"]][, cohort := "meta-analysis"]

### Combine significant results into one table
merge_cols <- c("module", "NPS", "cohort")
spec_cols  <- c("Estimate", "SE", "Pvalue")
meta_cols  <- c("fe_p", "fe_beta", "fe_se")

merged_sig <- merge(dats[["specific"]][Pvalue < 0.05, c(merge_cols, spec_cols), with = FALSE],
                    dats[["meta"]][fe_p < 0.05, c(merge_cols, meta_cols), with = FALSE],
                    by = merge_cols, all = TRUE)

## Reorder rows and columns
setorder(merged_sig, module, NPS, cohort)
setcolorder(merged_sig, spec_cols)


###### FINISH ######
fwrite(merged_sig, file = outfiles[["sig"]], row.names = FALSE, quote = FALSE, sep = "\t")
