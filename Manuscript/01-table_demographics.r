library(data.table)
library(ggplot2)


rm(list = ls())



###### GOAL ######




###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")

CountBinary <- function(dattab, x) {
    present <- length(which(dattab[, x, with = FALSE]==1))
    absent <- length(which(dattab[, x, with = FALSE]==0))
    percent <- paste0(sprintf("%.0f", (present / nrow(dattab)*100)), "%")
    missing <- length(which(is.na(dattab[, x, with=FALSE])))
    return(data.table(present, absent, percent, missing))
}


CountContinuous <- function(dattab, x, na.val=NA) {
    datvec <- unlist(dattab[, x, with = FALSE])
    datvec[datvec==na.val]  <- NA
    summvec <- summary(datvec)
    xmedian <-  sprintf("%.1f",summvec["Median"])
    xmean <- sprintf("%.1f", mean(datvec, na.rm=TRUE))
    xsd <- sprintf("%.1f", sd(datvec, na.rm=TRUE))
    mean.sd <- paste0(paste(xmean, xsd, sep=" ("), ")")
    xrange <- paste0(sprintf("%.1f", summvec["Min."]),
                    "-", #"\u2013", ## tried to use en dash but doesn't show up correctly in Excel
                    sprintf("%.1f", summvec["Max."]) )
    return(data.table(mean.sd = mean.sd,
                      median = xmedian,
                      range = xrange,
                      missing = ifelse("NA's" %in% names(summvec),
                                       as.vector(summvec["NA's"]), 0)))
}


## tabulate distribution of binary NPS variables
MakeNPSbintable <- function(dattab) {
    ## Tabulate
    disttab <- data.table(t(dattab[, sapply(.SD, function(x) {
        table(x, useNA = "always") }), .SDcol = npsvars.bin]))
    disttab[, nps:= npsvars.bin]
    setnames(disttab, c("0", "1", "V3"), c("no", "yes", "NA"))

    ## Add readable NPS names
    disttab <- merge(disttab, data.table(NPS = names(npsvars.bin),
                                         nps = npsvars.bin), by = "nps")
    disttab[, nps := NULL]

    ## Add total observations column
    ## Add percent cases column
    ## Add column for table with cases as count and percent in parentheses
    disttab[, total.obs := (no + yes)]
    disttab[, percent.cases := sapply((yes / (no + yes)*100), function(x) {
        paste0(sprintf("%.0f",x), "%") })]
    #disttab[, combined.col.cases := paste0(yes, " (", percent.cases, ")")]

    ## Reorder columns
    setnames(disttab, "yes", "cases")
    setcolorder(disttab, c("NPS", "cases", "percent.cases", "NA"))
    return(disttab[, 1:4])
}


###### SETUP ######


### Constants
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")
npsvars.sev <- paste0(npsvars.bin, "sev")
names(npsvars.sev) <- names(npsvars.bin)


### Files that exist
phenofiles <- sapply(cohorts, function(x) {
    paste0("2_Pipeline/02-Prepare_analysis_data/phenos_cleaned_SVs_", x, ".txt") })
residfiles <- sapply(cohorts, function(x) {
    paste0("2_Pipeline/02-Prepare_analysis_data/",
           "resid_regress_BatchPmiAgeMsex_", x, ".txt") })

### Files to be created
outdir <- "2_Pipeline/Manuscript/"
MyMkdir(outdir)
outfiles <- list(outfile = paste0(outdir, "table_demographics.txt"))


###### MAIN ######
### Read in data
phenos.all <- lapply(phenofiles, fread)
resids <- lapply(residfiles, fread, nrows = 1)

### Prepare tables - individual cohorts
## Get vector of samples in each resid file
resid_samps <- sapply(resids, function(x) { setdiff(names(x), "protein") } )

## Subset phenos to samples in corresponding resid file
phenos <- sapply(names(phenos.all), function(acohort) {
    phenos.all[[acohort]][protsample %in% resid_samps[[acohort]], ]
}, simplify = FALSE)


## Get total sample counts
totcounts <- sapply(phenos, nrow)

## Tabulate NPS
tabs.nps <- lapply(phenos, MakeNPSbintable)


## Tabulate binary covariates
covars.bin <- c("msex")
tabs.bin <- sapply(phenos, function(adattab) {   # iterate over cohorts
    rbindlist(sapply(covars.bin, function(avar) { # iterate over variables
        CountBinary(adattab, avar)
    }, simplify = FALSE), idcol = "var")
}, simplify = FALSE)


## Tabulate continuous covariates
covars.cont <- c("age_death", "pmi")
tabs.cont <- sapply(phenos, function(adattab) {   # iterate over cohorts
    rbindlist(sapply(covars.cont, function(avar) { # iterate over variables
        CountContinuous(adattab, avar)
    }, simplify = FALSE), idcol = "var")
}, simplify = FALSE)


### Test pairwise differences in age_death and pmi
### Add these manually to text
t.test(phenos[["OHSU"]]$age_death, phenos[["emory"]]$age_death)$p.value
t.test(phenos[["OHSU"]]$age_death, phenos[["rush"]]$age_death)$p.value
t.test(phenos[["rush"]]$age_death, phenos[["emory"]]$age_death)$p.value


t.test(phenos[["OHSU"]]$pmi, phenos[["emory"]]$pmi)$p.value
t.test(phenos[["OHSU"]]$pmi, phenos[["rush"]]$pmi)$p.value
t.test(phenos[["rush"]]$pmi, phenos[["emory"]]$pmi)$p.value


### Draw plots
## flatten data across cohorts to one table
combined_dt <- rbindlist(phenos, idcol = "cohort", fill = TRUE)
combined_dt[, cohort := factor(cohort, levels = cohorts)]

## Plot PMI histogram
p_pmi <- ggplot(combined_dt, aes(x = pmi)) +
    geom_histogram(binwidth = 0.5,   # adjust as needed
                   colour = "darkblue", fill = "darkblue") +
    facet_wrap(~ cohort, nrow = 1, scales = "free_y") +   # one panel per cohort
    scale_x_continuous(limits = c(0, NA),   # start at 0, upper limit auto
                       expand = c(0, 0)) + # no extra padding
    labs(title = "PMI",
         x = "PMI (hours)",
         y = "Count") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"),
          panel.spacing = unit(0.8, "lines"))





###### FINISH ######
## Write each table
write(paste0("Cohorts are ordered as follows: ",
             paste(cohorts, collapse = ", "), "\n"),
      file = outfiles[["outfile"]])

write("Total N:", file = outfiles[["outfile"]], append = TRUE)
fwrite(data.frame(rbind(totcounts)), file = outfiles[["outfile"]],
       col.names = TRUE, append = TRUE)

fwrite(Reduce(cbind, tabs.bin), file = outfiles[["outfile"]],
       col.names = TRUE, append = TRUE)
fwrite(Reduce(cbind, tabs.cont), file = outfiles[["outfile"]],
       col.names = TRUE, append = TRUE)
fwrite(Reduce(cbind, tabs.nps), file = outfiles[["outfile"]],
       col.names = TRUE, append = TRUE)
