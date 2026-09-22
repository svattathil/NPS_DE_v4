library(data.table)
library(argparse)
library(variancePartition)
library(gridExtra)

rm(list=ls())

### Constants
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")
covars <- c("msex", "Batch", "pmi", "age_death")

###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######
### Command line arguments
## Create parser
parser <- ArgumentParser()
parser$add_argument("--cohort", type="character", help="OHSU, rush, or emory")
parser$add_argument("--residset", type="character",
                    help="regressSVs or ignoreSVs or resid2_multi or noProtect", default = "ignoreSVs")

## Read from parser
args <- parser$parse_args()


### Files that exist
infiles <- list(
    log2norm = "2_Pipeline/01-Explore_data/prot.pcafiltered.log2norm.txt",
    phenos   = paste0("2_Pipeline/01-Explore_data/phenos_cleaned_", args$cohort, ".txt"),
    resid    = Residfile(args$residset, args$cohort)
)


### Files to be created
outdir <- paste0("2_Pipeline/02-Prepare_analysis_data/Plots_varpart/")
MyMkdir(outdir)
fileid <- paste0(args$residset, "_", args$cohort)
out.varpart <- paste0(outdir, "varpart_", fileid, ".pdf")


###### MAIN ######
### Read in data
log2norm.orig <- fread(infiles$log2norm)
phenos.orig <- fread(infiles$phenos)
resid.orig <- fread(infiles$resid)


### Prepare data
### Subset resid to complete cases
resid <- resid.orig[complete.cases(resid.orig), ]


### Subset all the other tables to the samples/proteins in resid
phenos <- phenos.orig[protsample %in% colnames(resid)]
log2norm <- log2norm.orig[protein %in% resid$protein, c("protein", phenos$protsample), with = FALSE]

## Reorder columns in resid to match phenos
## The proteins are already filtered
resid <- copy(resid[, c("protein", phenos$protsample), with = FALSE])

## Prepare phenos for varpart
varpartphenos <- phenos[, c(npsvars.bin, covars), with=FALSE]
varpartphenos[, Batch := paste0("Batch", Batch)]
setDF(varpartphenos)
rownames(varpartphenos) <- phenos$protsample

varpart.formula1 <- paste0("~", paste0(c(npsvars.bin, covars), collapse="+"))
vpres.log2norm <- fitExtractVarPartModel(log2norm[, !c("protein")], varpart.formula1, varpartphenos)
vpres.resid <- fitExtractVarPartModel(resid[, !c("protein")], varpart.formula1, varpartphenos)


p1 <- plotVarPart(vpres.log2norm, main = "Log2norm (Data before regressing anything)")
p2 <- plotVarPart(vpres.resid, main = paste0("Residset = ", args$residset))


## Draw variance partition plots
## for some reason quote/marrangeGrob complained when I put args$cohort and args$nps directly in the paste
acohort <- args$cohort
aresidset <- args$residset


forprint <- marrangeGrob(list(p1, p2), ncol=1, nrow=2, top=quote(paste(acohort)))
ggsave(plot=forprint, filename=out.varpart, width=10, height=16)


###### FINISH ######
