library(data.table)
library(argparse)
library(sva)
library(jaffelab)

rm(list = ls())
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### GOAL ######
## Regress out covariates
## Optionally, also regress out SVs
## Optionally, protect the NPS variables


##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")

nuisance_vars <- c("Batch", "pmi", "age_death", "msex")
covar_labels <- unlist(tstrsplit(nuisance_vars, split = "_", keep = 1))
covarstring <- paste0(Capwords(covar_labels), collapse = "")
npsvars <- c(npsvars.bin)


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######
### Command line arguments
## Create parser
parser <- ArgumentParser()
parser$add_argument("--cohort", type="character", default = "OHSU", help="OHSU, rush, or emory")
parser$add_argument("--option", type = "character", default = "regressSVs",
                    help = "'regressSVs' or 'ignoreSVs' or 'noProtect'")

## Read from parser
args <- parser$parse_args()


### Set fileid
fileid <- switch(args$option,
        regressSVs =  paste0(covarstring, "SVs"),
        ignoreSVs = covarstring,
        noProtect = paste0(covarstring, "_noProtect")
        )


### Files that exist
indir <- "2_Pipeline/01-Explore_data/"
infiles <- list(prot = paste0(indir, "prot.pcafiltered.log2norm.txt"),
                phenos = paste0(indir, "phenos_cleaned_", args$cohort, ".txt")
                )


### Files to be created
outdir <- "2_Pipeline/02-Prepare_analysis_data/"
MyMkdir(outdir)

outfiles <- list()
outfiles$resid <- paste0(outdir, "/resid_regress_", fileid, "_", args$cohort, ".txt")
outfiles$log   <- paste0(outdir, "/resid_regress_", fileid, "_", args$cohort, ".log")

## Currently this is written to the same file no matter the option
## It should be fine as long as the nuisance and nps variables are constant
outfiles$phenos.svs <- paste0(outdir, "/phenos_cleaned_SVs_", args$cohort, ".txt")


###### MAIN ######
### Read in data
phenos.orig <- fread(infiles$phenos)
prot.orig <- fread(infiles$prot)


### Prepare data
## convert batch to character to ensure it is treated as a categorical variable
phenos.orig[, Batch := paste0("Batch", Batch)]


### Subset proteomics data to samples in phenos file,
### which has gone through some filtering
### This also subsets to the current cohort
prot.sampfiltered <- prot.orig[, c("protein", phenos.orig$protsample), with = FALSE]


### Get number of non-missing observations per protein
obscounts <- apply(prot.sampfiltered[, !c("protein")], 1, function(x) { sum(!is.na(x)) })

### Filter the proteomics data to the proteins with observations for at least 50 samples
prot.countfiltered <- prot.sampfiltered[obscounts >= 50, ]


### Estimate SVs
## Subset phenos to samples with complete data for phenotype(s) of interest
nomissing <- phenos.orig[, apply(.SD, 1, function(x) { !any(is.na(x)) }),
                         .SDcol = c(nuisance_vars, npsvars)]
phenos <- phenos.orig[nomissing, ]

## filter phenos to samples that are also in the prot data
phenos <- phenos[protsample %in% colnames(prot.countfiltered), ]


## Again specify samples for proteomics data
## To account for this last bit of phenos filtering
## And to make sure sample order matches between columns of edata and rows of phenos
testsamples <- as.character(phenos$protsample)
edata.countfiltered <- as.matrix(prot.countfiltered[, ..testsamples])
rownames(edata.countfiltered) <- prot.countfiltered$protein

if(!(all(phenos$protsample == colnames(edata.countfiltered)))) {
    rm(edata.countfiltered)
} else {
    print("samples match between edata.countfiltered and phenos")
}

## Extract proteins with complete data, since that is what SVA requires
edata.complete <- edata.countfiltered[complete.cases(edata.countfiltered), ]

## Specify null and full models
nullmodel <- as.formula(paste("~ ", paste0(c(nuisance_vars), collapse = "+")))
fullmodel <- as.formula(paste("~ ", paste0(c(npsvars, nuisance_vars), collapse = "+")))


mod0 <- model.matrix(nullmodel, data = phenos)
mod  <- model.matrix(fullmodel, data = phenos)


## Apply SVA
## Step 1 - determine number of surrogate variables to estimate
                                        #set.seed(68302)
n.sv <- num.sv(edata.complete, mod, method="be", seed=68302)

## Step 2 - estimate the surrogate variables
## the component svobj$sv is a matrix whose columns correspond to
## the estimated surrogate variables
svobj <- sva(edata.complete, mod = mod, mod0 = mod0, n.sv = n.sv)
svs <- data.table(svobj$sv)
sv_names <- paste0("SV", 1:n.sv)
setnames(svs, sv_names)
svs[, protsample := colnames(edata.complete)]
setcolorder(svs, "protsample")


### Add SVs to phenos
phenos_withsvs <- merge(phenos, svs, by = "protsample", all.x = TRUE, sort = FALSE)

if(!(all(phenos_withsvs$protsample == colnames(edata.countfiltered)))) {
    rm(edata.countfiltered)
} else {
    print("samples match between edata.countfiltered and phenos_withsvs")
}


### Do cleaning
modelvars <- switch(args$option,
      regressSVs = c(npsvars, nuisance_vars, sv_names),
      ignoreSVs  = c(npsvars, nuisance_vars),
      noProtect = c(nuisance_vars)
      )

## np is the  number of variables to protect in cleaning (including intercept)
np <- switch(args$option,
   regressSVs = length(npsvars) + 1,
   ignoreSVs  = length(npsvars) + 1,
   noProtect  = 1
   )

modelform_forcleaning <- as.formula(paste("~ ", paste0(modelvars, collapse = "+")))
modmat_forcleaning <- model.matrix(modelform_forcleaning, data = phenos_withsvs)

## Make matching expression data including all proteins that passed count filtering
## (don't need to restrict to complete data)
edata.cleaned <- cleaningY(edata.countfiltered,
                           modmat_forcleaning,
                           P = np)

cleaned_forprint <- as.data.table(edata.cleaned, keep.rownames = "protein")


###### FINISH ######
### Write phenos with SVs
fwrite(phenos_withsvs, file = outfiles$phenos.svs, row.names = FALSE, quote = FALSE, sep = "\t")


### Write residuals
write.table(cleaned_forprint, file = outfiles$resid,
            row.names = FALSE, quote = FALSE, sep = "\t")


### Write log file
Writelog <- function(s, doappend=TRUE) { write(s, outfiles$log, append=doappend, ncolumns=1)}
Writelog(paste(c(args$cohort, paste0("number of SVs: ", n.sv)), collapse = "\t"), doappend=FALSE)
Writelog(paste0("Subjects dropped due to incomplete phenotype data: ", sum(!nomissing)))
Writelog(paste0("Subjects with SV estimates: ", sum(nomissing)))
Writelog(paste0("SVA null model: ", paste0(nullmodel, collapse = " ")))
Writelog(paste0("SVA full model: ", paste0(fullmodel, collapse = " ")))
Writelog(paste0("model for cleaning: ", paste0(modelform_forcleaning, collapse = " ")))
Writelog(paste0("terms protected during cleaning (including intercept): ", np))
