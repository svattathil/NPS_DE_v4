library(data.table)
library(argparse)
rm(list=ls())



###### FUNCTIONS ######




###### SETUP ######
### Constants
options(stringsAsFactors = FALSE)
regressvars <- c("Batch")


### Files that exist
outdir <- "2_Pipeline/01-Explore_data/"
file.log2norm <- paste0(outdir, "prot_before_sample_outlier_removal_log2norm.txt")
file.ids <- paste0(outdir, "ids_semi-filtered.txt")


### Files to be created
out.resid <- paste0(outdir, "prot_before_sample_outlier_removal_resid_regress", paste0(regressvars), ".txt")


###### MAIN ######
### Read in data
log2norm.orig <- fread(file.log2norm)
ids <- fread(file.ids)


### Prepare data
## Convert from data.table to matrix
log2norm <- copy(log2norm.orig)
protnames <- log2norm[, protein]
log2norm[, protein := NULL]
log2norm <- as.matrix(log2norm)
rownames(log2norm) <- protnames


### Convert variables to factor as necessary
ids[, Batch := paste0("batch", Batch)]
ids[, Batch := factor(Batch)]


### Define regression model
regressmodel <- paste0("protein ~ ", paste(regressvars, sep="+"))


### Estimate residuals for each protein and save results in list (one element per protein)

## Initialize empty matrix
resids <- matrix(NA, nrow=nrow(log2norm), ncol = ncol(log2norm))
rownames(resids) <- rownames(log2norm)
colnames(resids) <- colnames(log2norm)

## fit model for each protein in turn
for(i in (1:nrow(log2norm))) {
    ## Print intermittent status line
    statusfreq <- 500
    if((i - (i %/% statusfreq)*statusfreq) == 0) { print(paste0("protein ", i, " of ", nrow(log2norm))) }

    ## Get name of current protein
    currprot <- rownames(log2norm)[i]

    ## Extract observations for current protein
    currobs <- data.table(protsample = colnames(log2norm), protein = log2norm[currprot, ])

    ## Merge in phenos
    merged <- merge(currobs, ids, by = "protsample")

    ## restore proper sample order
    merged <- merged[match(colnames(resids), merged$protsample), ]

    ## Fit model
    fit <- lm(regressmodel, data = merged, na.action = na.exclude)

    ## Extract residuals
    resids[i, ] <- residuals(fit)

    ## Clean up
    rm(list=c("currprot", "currobs", "merged", "fit"))
}
rm(i)


### Convert to data.table for printing
forprint <- data.frame(resids)
setDT(forprint)
forprint[, protein := rownames(resids)]
setcolorder(forprint, "protein")


###### FINISH ######
write.table(forprint, file = out.resid, row.names = FALSE, quote = FALSE, sep = "\t")
