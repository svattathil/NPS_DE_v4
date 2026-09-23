library(data.table)
library(argparse)
library(variancePartition)
library(gridExtra)

rm(list=ls())

### Constants
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")
covars <- c("msex", "Batch", "pmi", "age_death")

residsets <- c("resid2_multi", "noProtect", "ignoreSVs", "regressSVs")

###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######
### Command line arguments
## Create parser
parser <- ArgumentParser()
parser$add_argument("--cohort", type="character", default = "OHSU", help="OHSU, rush, or emory")

## Read from parser
args <- parser$parse_args()


### Files that exist
infiles <- list(
    log2norm = "2_Pipeline/01-Explore_data/prot.pcafiltered.log2norm.txt",
    phenos   = paste0("2_Pipeline/01-Explore_data/phenos_cleaned_", args$cohort, ".txt")
)
for(x in residsets) {
        infiles[[x]]    = Residfile(x, args$cohort)
}


### Files to be created
outdir <- paste0("2_Pipeline/02-Prepare_analysis_data/Plots_varpart/")
MyMkdir(outdir)
fileid <- paste0(args$cohort)
out.varpart <- paste0(outdir, "varpart_", fileid, ".pdf")


###### MAIN ######
### Read in data
log2norm.orig <- fread(infiles$log2norm)
phenos.orig <- fread(infiles$phenos)
resid_list <- sapply(residsets, function(x) { fread(infiles[[x]]) }, simplify = FALSE)


### Prepare data
### Subset resids to complete cases
for(residset in residsets) {
resid_list[[residset]] <-
    resid_list[[residset]][complete.cases(resid_list[[residset]]), ]
}
rm(residset)


### Subset all the other tables to the samples/proteins in first resid
### I am assuming that all the resids use the same proteins
phenos <- phenos.orig[protsample %in% colnames(resid_list[[1]])]
log2norm <- log2norm.orig[protein %in% resid_list[[1]]$protein,
                          c("protein", phenos$protsample), with = FALSE]

## Reorder columns in resids to match phenos
## The proteins are already filtered
for(residset in residsets) {
    resid_list[[residset]] <-
        resid_list[[residset]][, c("protein", phenos$protsample), with = FALSE]
}
rm(residset)

## Prepare phenos for varpart
varpartphenos <- phenos[, c(npsvars.bin, covars), with=FALSE]
varpartphenos[, Batch := paste0("Batch", Batch)]
setDF(varpartphenos)
rownames(varpartphenos) <- phenos$protsample

varpart.formula1 <- paste0("~", paste0(c(npsvars.bin, covars), collapse="+"))
vpres.log2norm <- fitExtractVarPartModel(log2norm[, !c("protein")], varpart.formula1, varpartphenos)


vpres.resid.list <- sapply(residsets, function(x) {
    fitExtractVarPartModel(resid_list[[x]][, !c("protein")], varpart.formula1, varpartphenos)
}, simplify = FALSE)


p0 <- plotVarPart(vpres.log2norm, main = "Log2norm (Data before regressing anything)")
p1 <- plotVarPart(vpres.resid.list[[1]], main = paste0("Residset = ", residsets[[1]]))
p2 <- plotVarPart(vpres.resid.list[[2]], main = paste0("Residset = ", residsets[[2]]))
p3 <- plotVarPart(vpres.resid.list[[3]], main = paste0("Residset = ", residsets[[3]]))
p4 <- plotVarPart(vpres.resid.list[[4]], main = paste0("Residset = ", residsets[[4]]))

## Draw variance partition plots
## for some reason quote/marrangeGrob complained when I put
## args$cohort and args$nps directly in the paste
acohort <- args$cohort


forprint <- marrangeGrob(list(p0, NA, p1, p3, p2, p4), ncol=3, nrow=2, top=quote(paste(acohort)))
ggsave(plot=forprint, filename=out.varpart, width=20, height=16)


###### FINISH ######
