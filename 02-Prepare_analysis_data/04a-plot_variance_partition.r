library(data.table)
library(argparse)
library(variancePartition)
library(gridExtra)

rm(list=ls())

### Constants
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")
residsets <- c("resid2_multi", "noProtect", "ignoreSVs", "regressSVs")
covars <- c("msex", "Batch", "pmi", "age_death")


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######
### Command line arguments
## Create parser
parser <- ArgumentParser()
parser$add_argument("--cohort", type="character", default = "emory", help="OHSU, rush, or emory")
parser$add_argument("--useSVs", action="store_true", default=FALSE, help="")


## Read from parser
args <- parser$parse_args()


### Files that exist
## general and cohort-specific files
infiles <- list(
    log2norm = "2_Pipeline/01-Explore_data/prot.pcafiltered.log2norm.txt",
    phenos   = paste0("2_Pipeline/02-Prepare_analysis_data/phenos_cleaned_SVs_",
                      args$cohort, ".txt")
)

## resid-specific and cohort-specific files
for(x in residsets) {
    infiles[[x]] = Residfile(x, args$cohort)
}


### Files to be created
outdir <- paste0("2_Pipeline/02-Prepare_analysis_data/Plots_varpart/")
MyMkdir(outdir)
fileid <- ifelse(args$useSVs, paste0(args$cohort, "_SVs"), args$cohort)
outfiles <- list(pdf = paste0(outdir, "varpart_", fileid, ".pdf"))


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
svcols <- grep("SV", names(phenos), value = TRUE)
varpartphenos <- phenos[, c(npsvars.bin, covars, svcols), with=FALSE]
varpartphenos[, Batch := paste0("Batch", Batch)]
setDF(varpartphenos)
rownames(varpartphenos) <- phenos$protsample


### Define formula
varpart.formula1 <- paste0("~", paste0(c(npsvars.bin, covars), collapse="+"))
if(args$useSVs) {
    varpart.formula1 <- paste0(varpart.formula1, "+", paste0(svcols, collapse="+"))
}


### Extract variance partition
vpres.log2norm <- fitExtractVarPartModel(log2norm[, !c("protein")],
                                         varpart.formula1, varpartphenos)

vpres.resid.list <- sapply(residsets, function(x) {
    fitExtractVarPartModel(resid_list[[x]][, !c("protein")], varpart.formula1, varpartphenos)
}, simplify = FALSE)


### Make plots
p0 <- plotVarPart(vpres.log2norm, main = "Log2norm (Data before regressing anything)")
p1 <- plotVarPart(vpres.resid.list[[1]], main = paste0("Residset = ", residsets[[1]]))
p2 <- plotVarPart(vpres.resid.list[[2]], main = paste0("Residset = ", residsets[[2]]))
p3 <- plotVarPart(vpres.resid.list[[3]], main = paste0("Residset = ", residsets[[3]]))
p4 <- plotVarPart(vpres.resid.list[[4]], main = paste0("Residset = ", residsets[[4]]))


### Draw plots to file
## for some reason quote/marrangeGrob complained when I put
## args$cohort and args$nps directly in the paste
width <- ifelse(args$useSVs, 26, 20)
height <- 16
acohort <- args$cohort

forprint <- marrangeGrob(grobs = list(p0, p1, p2, p3, p4),
                         layout_matrix = matrix(c(1:3, NA, 4:5), byrow = TRUE, nrow = 2),
                         top=quote(paste(acohort)))

ggsave(plot=forprint, filename=outfiles$pdf, width=width, height=height)


###### FINISH ######
