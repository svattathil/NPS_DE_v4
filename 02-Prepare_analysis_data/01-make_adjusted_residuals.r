library(data.table)
library(argparse)
library(sva)
library(jaffelab)

rm(list = ls())
source("~/comsv/Svattathil_Library/svattathil_functions.r")



###### GOAL ######
## Regress out covariates and SVs while protecting NPS


##### CONSTANTS #####
options(stringsAsFactors = FALSE)
npsvars.bin <- c(agitation = "agit", anxiety = "anx", apathy = "apa", delusion = "del",
                      depression = "depd", disinhibition = "disn", hallucination = "hall",
                      irritability = "irr", sleep = "nite")

vars_to_protect <- c(npsvars.bin, "msex")
vars_to_regress <- c("Batch", "pmi", "age_death")
covarstring <- paste0(paste0(Capwords(sapply(covars function(x) {
    unlist(strsplit(x, split = "_"))[1] })), collapse=""), "SVs")


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
indir <- "2_Pipeline/01-Explore_data/"
infiles <- list(prot = paste0(indir, "prot.pcafiltered.log2norm.txt"),
                phenos = paste0(indir, "phenos_cleaned_", args$cohort, ".txt")
                )


### Files to be created
outdir <- paste0("2_Pipeline/02-Prepare_analysis_data/Resid4_jaffe_regress", covarstring)
MyMkdir(outdir)
outfiles <- list(out.resid = paste0(outdir, "/resid_regress",
                                    covarstring, "_", args$cohort, ".txt"),
                 out.svs = paste0(outdir, "/svs_", args$cohort, ".txt"),
                 out.log = paste0(outdir, "/", args$cohort, ".log"))


###### MAIN ######
### Read in data
phenos.orig <- fread(infiles$phenos)
prot.orig <- fread(infiles$prot)


### Prepare data
## convert batch to character to ensure it is treated as a categorical variable
phenos.orig[, Batch := paste0("Batch", Batch)]


### Subset proteomics data to samples in phenos file,
### which has gone through some filtering
prot.sampfiltered <- prot.orig[, c("protein", phenos.orig$protsample), with = FALSE]

### Get number of non-missing observations per protein
obscounts <- apply(prot.sampfiltered[, !c("protein")], 1, function(x) { sum(!is.na(x)) })

### Filter the proteomics data to the proteins with observations for at least 50 samples
prot.countfiltered <- prot.sampfiltered[obscounts >= 50, ]


### Estimate SVs
## Subset phenos to samples with complete data for phenotype(s) of interest
nomissing <- phenos.orig[, apply(.SD, 1, function(x) { !any(is.na(x)) }),
                         .SDcol = c(vars_to_regress, vars_to_protect)]
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
nullmodel <- as.formula(paste("~ ", paste0(c(vars_to_regress), collapse = "+")))
fullmodel <- as.formula(paste("~ ", paste0(c(vars_to_protect, vars_to_regress), collapse = "+")))

mod0 <- model.matrix(nullmodel, data = phenos)
mod <- model.matrix(fullmodel, data = phenos)


## Apply SVA
## Step 1 - determine number of surrogate variables to estimate
#set.seed(68302)
n.sv <- num.sv(edata.complete, mod, method="be", seed=68302)

## Step 2 - estimate the surrogate variables
## the component svobj$sv is a matrix whose columns correspond to
## the estimated surrogate variables
svobj <- sva(edata.complete, mod, mod0, n.sv=n.sv)
svs <- data.table(svobj$sv)
sv_names <- paste0("SV", 1:n.sv)
setnames(svs, sv_names)
svs[, protsample := colnames(edata.complete)]
setcolorder(svs, "protsample")


### Do cleaning
## Add SVs to phenos
phenos_withsvs <- merge(phenos, svs, by = "protsample", all.x = TRUE, sort = FALSE)
modelform_forcleaning <- as.formula(paste("~ ",
                                          paste0(c(vars_to_protect,
                                                   vars_to_regress ,
                                                   sv_names),
                                                 collapse = "+")))
modmat_forcleaning <- model.matrix(modelform_forcleaning, data = phenos_withsvs)


## Make matching expression data including all proteins that passed count filtering
## (don't need to restrict to complete data)
edata.cleaned <- cleaningY(edata.countfiltered,
                           modmat_forcleaning,
                           P = length(vars_to_protect + 1))


cleaned_forprint <- as.data.table(edata.cleaned, keep.rownames = "protein")


###### FINISH ######

### Write SVs
fwrite(svs, file = outfiles$out.svs, row.names = FALSE, quote = FALSE, sep = "\t")

### Write residuals
write.table(cleaned_forprint, file = outfiles$out.resid,
            row.names = FALSE, quote = FALSE, sep = "\t")


### Write log file
Writelog <- function(s, doappend=TRUE) { write(s, outfiles$out.log, append=doappend, ncolumns=1)}
Writelog(paste(c(args$cohort, paste0("number of SVs: ", n.sv)), collapse = "\t"), doappend=FALSE)
Writelog(paste0("Subjects dropped due to incomplete phenotype data: ", sum(!nomissing)))
Writelog(paste0("Subjects with SV estimates: ", sum(nomissing)))
Writelog(paste0("SVA full model: ", paste0(fullmodel, collapse = " ")))
Writelog(paste0("SVA null model: ", paste0(nullmodel, collapse = " ")))
