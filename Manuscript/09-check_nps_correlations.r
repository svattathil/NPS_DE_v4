library(data.table)
library(argparse)
library(corrplot)
rm(list=ls())


### Constants ###
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")

ScatterDiscrete <- function(adat, var1, var2 , coeff_bigger=2){
    ## modified from https://r-graph-gallery.com/5-correlation-of-discrete-variables.html
    ## Get complete cases
    forplot  <- adat[complete.cases(adat[, c(var1, var2), with = FALSE]),
                     c(var1, var2), with = FALSE]
    countsres <- xyTable(unlist(forplot[, var1, with = FALSE]),
                         unlist(forplot[, var2, with = FALSE]))

    ## Draw points with magnification
    plot(countsres$x , countsres$y , cex=countsres$number*coeff_bigger, pch=16,
         col="chocolate1", xlab= var1, ylab=var2,
         main = paste0(var1, " v. ", var2),
         axes = FALSE,
         sub = paste0("Using ", pN(nrow(forplot)), " samples with complete data"))

    ## add axes
    axis(side = 1, at = seq(0,3))
    axis(side = 2, at = seq(0,3))

    ## Add counts as overlaid text
    text(countsres$x, countsres$y, countsres$number)
}


###### SETUP ######
### Files that exist
indir <- "2_Pipeline/02-Prepare_analysis_data/"
infiles <- sapply(cohorts, function(cohort) {
    paste0(indir, "phenos_cleaned_SVs_", cohort, ".txt")
    }, simplify = FALSE)



### Files to be created
outdir <- "2_Pipeline/Manuscript/"
outfiles <- list(txt =  paste0(outdir, "nps_correlations.txt"),
                 pdf =  paste0(outdir, "nps_correlations.pdf"),
                 png =  paste0(outdir, "nps_correlations.png")
                 )


###### MAIN ######
### Read in data
phenos_list <- sapply(infiles, fread, simplify = FALSE)


### Prepare data
## Extract phenotype columns (exclude ID columns)
cordat <- sapply(phenos_list, function(dat) {
    dat[, c(npsvars.bin), with = FALSE]
}, simplify = FALSE)


### Calculate pairwise correlation
## Pearson correlation for pairs of binary variables is the
## same as Phi, which is the appropriate statistics for two binary variables
phires <- sapply(cordat, function(dat) { cor(dat, method = "pearson") }, simplify = FALSE)


### Plot pairwise correlation
pdf(outfiles$pdf, width = 8, height = 8)
par(oma = c(1, 1, 1, 1), mar = c(7, 7, 3, 1), mfrow = c(2,2))

for(cohort in cohorts) {
    corrplot(phires[[cohort]], method = "number")
    title(cohort, outer = FALSE)
    title("phi correlation", outer = TRUE)
}

if(0) {
### Make biplots for selected pair of traits
    for(x in 1:(length(npsvars)-1) ){
        for(y in (x+1):length(npsvars)) {
            ScatterDiscrete(phenos,
                            paste0(npsvars[x], ".sevmax"),
                            paste0(npsvars[y], ".sevmax"), coeff_bigger = 0.005)
        }
    }
}
dev.off()


###### FINISH ######
## Write correlation matrices to file
write.table(phires, file = outfiles$txt,
            row.names = TRUE, quote = FALSE, sep = "\t")


