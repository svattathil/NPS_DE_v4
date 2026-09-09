library(data.table)
library(argparse)
rm(list=ls())


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######
### Constants
options(stringsAsFactors = FALSE)

sourceset <- "Rush"


## Define nps variable names -- all 9 domains
npsvars.bin.9 <- c(agitation = "agit", anxiety = "anx", apathy = "apa", delusion = "del",
                 depression = "depd", disinhibition = "disn", hallucination = "hall",
                 irritability = "irr", sleep = "nite")

npsvars.sev.9 <- paste0(npsvars.bin.9, "sev")
names(npsvars.sev.9) <- names(npsvars.bin.9)

## Define nps variable names -- combine hall/del and drop disn
npsvars.bin.7 <- c(agitation = "agit", anxiety = "anx", apathy = "apa",
                 depression = "depd", irritability = "irr", psychosis = "psych",
                  sleep = "nite")

npsvars.sev.7 <- paste0(npsvars.bin.7, "sev")
names(npsvars.sev.7) <- names(npsvars.bin.7)


### Map each dataset to the variable set that it uses
### Note the rad201_v1 and _v2 use the same variables, but have different variable codings!
varsets <- data.frame(matrix(c("rad201", "rad201_v1",
                               "rad201", "rad201_v2",
                               "rad20x", "rad20x",
                               "radintr", "radintr"), ncol=2, byrow=TRUE))
setDT(varsets)
setnames(varsets, c("varset", "pheset"))
setkey(varsets, "pheset")


### Files that exist
## Original proteomics id file
file.ids <- "2_Pipeline/01-Explore_data/ids_pcafiltered.txt"

## Original phenotype files
filedir <- "0_Data/SharepointFiles/"
phesets <- c("rad201_v1", "rad201_v2", "rad20x", "radintr", "Wingo_Data")
phesets.nps <- setdiff(phesets, "Wingo_Data")
files.phenos <-paste0(filedir, "Phenotypes/Rush/", phesets, ".txt")
file.ids <- "2_Pipeline/01-Explore_data/ids_pcafiltered.txt"

## File with selected variables and assigned NPS that we assembled manualy
file.selectedvars <- "0_Data/Rush_variables/rush_variables_selected_grouped_set4.txt"


### Files to be created
outdir <- "2_Pipeline/01-Explore_data/"

out.cleanedphenos <- paste0(outdir, "phenos_cleaned_", tolower(sourceset), ".txt")
out.caserates <- paste0(outdir, "rush_caserates", ".txt")
out.hists <- paste0(outdir, "rush_composite_score_distributions.pdf")
out.check_long_del_disn <- paste0(outdir, "rush_check_long_del_disn.txt")


###### MAIN ######
### Read in data
ids.all <- fread(file.ids)
phenos <- lapply(files.phenos, fread)
names(phenos)  <- phesets

selectedvars <- fread(file.selectedvars)


### Prepare data
## Drop the one variable that is reverse coded
selectedvars <- copy(selectedvars[is.na(reverse), ])

## Extract ids for this source
ids <- ids.all[source=="Rush", ]
ids[, ID := as.integer(ID)]

## For radintr variables, replace "NULL" text string with NA and convert to integer
## Remove variable with character values first
phenos[["radintr"]][, onwht := NULL]
phenos[["radintr"]][, colnames(phenos[["radintr"]]) := lapply(.SD, function(x) {
    x[x=="NULL"] <- NA; return(as.integer(x))  })]


## Change case of selectedvars to match phenotype data
selectedvars[, label_201 := tolower(label_201)]
selectedvars[, label_20x := tolower(label_20x)]

### Extract phenotypes for selected variables for samples that passed QC
phenos.ana <- sapply(1:nrow(varsets), function(i) {
    cat("varset ", i, "\n") ## for TEST

    ## define current variable set and pheset
    varset <- varsets[i, varset]
    pheset <- varsets[i, pheset]

    ## 1. Identify selected vars from table. These were selected based on the codebook
    ## 2. Identify the ones that are actually in the dataset
    ## 3. Count number of missing variables
    ##    Make sure NA is not counted
    labelcol <- paste0("label_", gsub("rad", "", varset))
    vars.selected <- selectedvars[, get(labelcol)]
    vars.found <- intersect(vars.selected, names(phenos[[pheset]]))

    vars.missing <- setdiff(setdiff(vars.selected, names(phenos[[pheset]])), NA)
    print(paste0("missing ", length(vars.missing), " variables: ",
                 paste0(vars.missing, collapse=", ")))

    ## Extract variables and samples of interest
    toreturn <- phenos[[pheset]][projid %in% ids$ID, c("projid", "age_int", vars.found),
                                 with = FALSE]

    ## Fill in NAs for missing variables, if there are any
    if(length(vars.missing) > 0) {
        toreturn[, c(vars.missing) := NA_integer_]
    }

    ## Rename with label from rad20x dataset, so all phesets will have the same columns
    if(varset != "rad20x") {
        ## Extract the two sets of labels
        labelmap <- selectedvars[, c("label_20x", labelcol), with = FALSE]
        setnames(labelmap, labelcol, "oldlabel")

        ## Filter to variables that were found or filled in with NAs
        labelmap <- labelmap[oldlabel %in% names(toreturn), ]

        ## Replace names
        setnames(toreturn, labelmap$oldlabel, labelmap$label_20x)
    }

    ## Sort by age_int
    setorder(toreturn, age_int)
    return(toreturn)
}, simplify = FALSE)
names(phenos.ana) <- varsets$pheset


### Get vector of variables available in each pheset
varcols <- sapply(varsets$pheset, function(apheset) {
    setdiff(names(phenos.ana[[apheset]]), c("projid", "age_int")) }, simplify = FALSE)


### Recode all variables to have 0=No/1=Yes
Recode <- function(x, oldvals, newval) { x[x %in% oldvals] <- newval; return(x) }
Showtable <- function(y) { sapply(phenos.ana[[y]][, varcols[[y]], with = FALSE], table) }

## rad201_v1 already uses 0=No/1=Yes, with 9 for missing
apheset <- "rad201_v1"
#Showtable(apheset)

## rad201_v2 uses 2=No/1=Yes, with 7 and 8 for missing
apheset <- "rad201_v2"
phenos.ana[[apheset]][, c(varcols[[apheset]]) := lapply(.SD, Recode, oldvals=c(2), newval=0),
                      .SDcol = varcols[[apheset]]]

## rad20x uses 2=No/1=Yes, with 7, 8, and 9 for missing
apheset <- "rad20x"
phenos.ana[[apheset]][, c(varcols[[apheset]]) := lapply(.SD, Recode, oldvals=c(2), newval=0),
                      .SDcol = varcols[[apheset]]]

## radintr uses 2=No/1=Yes, with 9 for missing
## Also, variable apathy2 had slightly different wording here
## and should be recoded from c(1,2)=No/3=Yes
apheset <- "radintr"
if("apathy2" %in% names(phenos.ana[[apheset]])) {
    phenos.ana[[apheset]][apathy2 %in% c(1,2), apathy2 := 0]
    phenos.ana[[apheset]][apathy2 %in% c(3), apathy2 := 1]
}
phenos.ana[[apheset]][, c(varcols[[apheset]]) := lapply(.SD, Recode, oldvals=c(2), newval=0),
                      .SDcol = varcols[[apheset]]]


### Fill in NAs for all observations that are not 0 or 1
FillNA <- function(x) { x[!x %in% c(0,1)] <- NA_integer_; return(x) }

for(apheset in varsets$pheset) {
    phenos.ana[[apheset]][, c(varcols[[apheset]]) := lapply(.SD, FillNA),
                          .SDcol = varcols[[apheset]]]
}
rm(apheset)


### Collapse all phesets into one table
mergedvars <- rbindlist(phenos.ana, idcol = "pheset", fill=TRUE)


### Extract samples that passed QC from Wingo_Data table
phenos.ana[["Wingo_Data"]] <- phenos[["Wingo_Data"]][projid %in% ids$ID, ]


### Extract last visit observations
setorder(mergedvars, age_int)
lastobs1 <- mergedvars[, tail(.SD, 1), by = "projid"]


### Count the number of NAs per variable
### and drop the variables that have complete missing
tmp.nonmissing <- apply(lastobs1, 2, function(x) { sum(!is.na(x)) })
vars.complete.missing.todrop <- names(tmp.nonmissing)[tmp.nonmissing==0]
if(length(vars.complete.missing.todrop) >0) {
    lastobs1[, c(vars.complete.missing.todrop) := NULL]
}


### Split observations into one table per NPS
lastobs1.bynps <- sapply(unique(selectedvars$NPS), function(anps) {
    ## Have to account for the variables that were dropped due to complete missingness
    vars.tokeep <- setdiff(selectedvars[NPS == anps, label_20x], vars.complete.missing.todrop)
    print(paste0("keep ", length(vars.tokeep), " vars for ", anps))
    lastobs1[, c("projid", "age_int", "pheset", vars.tokeep), with = FALSE]
}, simplify = FALSE)


### Count the number of variables per NPS
nvars <- sapply(lastobs1.bynps, ncol) - 3


### Sum the observations per NPS to get composite score
sumscores.list <- sapply(lastobs1.bynps, function(x) {
    ## Calculate score sums and score means
    sums <- apply(x[, !c("projid", "age_int", "pheset"), with = FALSE], 1, sum, na.rm=TRUE)
    means <- apply(x[, !c("projid", "age_int", "pheset"), with = FALSE], 1, mean, na.rm=TRUE)

    ## Combine into one table with projids
    toreturn <- data.table(projid = x$projid, sumscore = sums, meanscore = means)

    ## Set NA for the participants that had missing value for all variables
    all.na <- apply(x[, !c("projid", "age_int", "pheset"), with = FALSE], 1, function(x) {
        all(is.na(x)) })
    toreturn[all.na, c("sumscore", "meanscore") := NA]

    ## Return table
    return(toreturn)
}, simplify = FALSE)


### Calculate case/control rates
## Define score threshold
## Participants with scores higher than this will be classified as cases
meanthresh <- 0

for(anps in names(sumscores.list)) {
    sumscores.list[[anps]][, c(anps) := ifelse(meanscore > meanthresh, 1, 0)]
}

if(0) {
    ## For apathy, instead of using a threshold, we're going to split based on a percentile
    percentilethresh <- 0.3
    p30 <- quantile(sumscores.list[["apathy"]]$meanscore, na.rm=TRUE, probs=percentilethresh)

    sumscores.list[["apathy"]][, apathy := NA]
    sumscores.list[["apathy"]][, apathy := ifelse(meanscore <= p30, 0, 1)]
}

caserates <- data.frame(t(sapply(names(sumscores.list), function(anps) {
    return(c(
        NPS = anps,
        percent.cases = MakeFrac(mean(sumscores.list[[anps]][, get(anps)], na.rm=TRUE),
                                 dec.places=1),
        no = sum(sumscores.list[[anps]][, get(anps)] == 0, na.rm = TRUE),
        yes = sum(sumscores.list[[anps]][, get(anps)] == 1, na.rm = TRUE),
        na = sum(is.na(sumscores.list[[anps]][, get(anps)]))
    ))
})))
setDT(caserates)
setnames(caserates, "na", "NA")
setorder(caserates, percent.cases)
caserates[, percent.cases := paste0(percent.cases, "%")]


### Assemble final dataset
## Assemble NPS
formerge <- lapply(names(sumscores.list), function(anps) {
    sumscores.list[[anps]][, c("projid", anps), with = FALSE] })
npsobs <- Reduce(function(x,y) { merge(x, y, by = "projid") }, formerge)

## Update NPS names
setnames(npsobs, names(npsvars.bin.9), npsvars.bin.9)

## Add Braak score
phenos.forprint <- merge(phenos.ana[["Wingo_Data"]][, .(projid, braaksc)], npsobs,
                         all.y=TRUE, by="projid")

## Add variables from id table
phenos.forprint <- merge(ids[, .(ID, protsample, Batch, age_death, msex, pmi)], phenos.forprint,
                         by.x="ID", by.y = "projid")



### Add psychosis variable
phenos.forprint[hall == 1 | del == 1, psych := 1]
phenos.forprint[hall == 0 & del == 0, psych := 0]
phenos.forprint[hall == 0 & is.na(del), psych := 0]
phenos.forprint[del == 0 & is.na(hall), psych := 0]


### Do some checks - These don't directly contribute to the case definitions ###
### But they help guide and interpret
### Visualize
pdf(out.hists, height = 7, width = 7)
par(mfrow=c(3,3), oma = c(0,0,2,0))

invisible(sapply(names(sumscores.list), function(anps) {
    toplot <- sumscores.list[[anps]]$sumscore
    hist(toplot, breaks=50, col="purple2",
         main = anps,
         xlab = paste0("Sum of ", nvars[anps], " variables"))
    mtext(paste0("nmissing = ", sum(is.na(toplot))), side=3, line=0.2, cex=0.8)
}))
title(paste0("Sum of scores per NPS for ", nrow(lastobs1), " participants"),
      outer = TRUE, line = 0)


invisible(sapply(names(sumscores.list), function(anps) {
    toplot <- sumscores.list[[anps]]$meanscore
    hist(toplot, breaks=50, col="orange4",
         main = anps,
         xlab = paste0("Mean of ", nvars[anps], " variables"))
    mtext(paste0("nmissing = ", sum(is.na(toplot))), side=3, line=0.2, cex=0.8)
}))
title(paste0("Mean score per NPS for ", nrow(lastobs1), " participants"), outer = TRUE, line = 0)
dev.off()


### Check missing rate by pheset for the variables for the two NPS with the high missing rate
## Check delusion items
sapply(setdiff(names(lastobs1.bynps[["delusion"]]), c("projid", "age_int", "pheset")),
       function(itemx) {
    lastobs1.bynps[["delusion"]][, mean(is.na(get(itemx))), by = "pheset"]
}, simplify = FALSE)


## Check disinhibition items
sapply(setdiff(names(lastobs1.bynps[["disinhibition"]]), c("projid", "age_int", "pheset")),
       function(itemx) {
    lastobs1.bynps[["disinhibition"]][, mean(is.na(get(itemx))), by = "pheset"]
}, simplify = FALSE)


### Extract long obs for the items for the two NPS with the high missing rate
## Add the value for our derived last-visit variable
forcheck_long_del_disn <- mergedvars[, c("projid", "age_int", "pheset",
                                         selectedvars[NPS %in% c("delusion", "disinhibition"), label_20x]),
                                     with = FALSE]
forcheck_long_del_disn <- merge(npsobs[, .(projid, del, disn)], forcheck_long_del_disn, by = "projid")
setorder(forcheck_long_del_disn, "projid", "age_int")

### Check distribution for each item used for apathy
apathy_items <- selectedvars[NPS == "apathy", label_20x]

tab_apathy_items <- lastobs1.bynps$apathy[, ..apathy_items] |>
apply(X = _, MARGIN = 2, FUN = table, useNA = "ifany") |>
addmargins(A = _, margin = 1)




###### FINISH ######
write.table(caserates, file = out.caserates,
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(phenos.forprint, file = out.cleanedphenos,
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(forcheck_long_del_disn, file = out.check_long_del_disn,
            row.names = FALSE, quote = FALSE, sep = "\t")
