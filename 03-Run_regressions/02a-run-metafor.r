library(data.table)
library(argparse)
library(metafor)
library(qvalue)

rm(list = ls())


###### GOAL ######
## Run cross-cohort meta-analysis for each NPS domain for specified model


##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")
source("~/comsv/Svattathil_Library/manhattan_and_qq_plot_functions.r")


###### SETUP ######
### Command line arguments
## Create parser
parser <- ArgumentParser()
parser <- ArgumentParser(description = "limma DE of log2 protein abundance on NPS domains")

## Control model and sample set
parser$add_argument("--analysis", default = "basic",
                    help = paste0("Specify the analysis (e.g., basic, severity, ",
                                  "males, females, e4_adj, latfactor)"))
parser$add_argument("--svchoice", default = "noSV",
                    help = paste0("withSV or noSV"))
parser$add_argument("--SE_type", type="character", default="SE_mod", help = "SE_mod or SE_unmod")


## Read from parser
args <- parser$parse_args()


## Check arguments
CheckOptions("analysis", c("basic", "severity", "males", "females", "e4_adj", "latfactor", "cond_ind"))
CheckOptions("svchoice", c("withSV", "noSV"))

## Define runid
runid <- paste0(args$analysis, "_", args$svchoice)

## Define outcomes and subset of cohorts (if necessary)
if(args$analysis %in% c("basic", "noSV", "cond_ind",  "e4_adj")) {
    npsvars <- npsvars.bin
}

if(args$analysis %in% c("severity")) {
    npsvars <- npsvars.sev
    cohorts <- c("OHSU", "emory")
}


### Files that exist
if(args$analysis == "cond_ind") {
    ## the cond_ind analysis is done inside the Basic analysis
    resdir <- paste0("2_Pipeline/03-Run_regressions/Prot_perDomain/", paste0("Basic_", args$svchoice), "/")
} else {
    resdir <- paste0("2_Pipeline/03-Run_regressions/Prot_perDomain/", Capwords(runid), "/")
}

Statsfile <- function(cohort) { paste0(resdir, runid, "_", cohort, ".tsv") }


### Files to be created
outdir <- paste0(resdir, "Metafor/")
MyMkdir(outdir)
outfiles <- list(stats = paste0(outdir, "metafor_", runid, "_using_", args$SE_type, ".txt"),
                 qq    = paste0(outdir, "metafor_", runid, "_using_", args$SE_type, "_qqplots.pdf")
                 )


###### MAIN ######
### Read in data
res.list.cohort <- sapply(cohorts, function(x) { fread(Statsfile(x)) }, simplify = FALSE)


### Prepare data
## Split by domain
res.list.nps <- sapply(as.character(npsvars), function(anps) {
    sapply(res.list.cohort, function(cohortdat) {
        cohortdat[domain == anps, ]
        }, simplify = FALSE)
    }, simplify = FALSE)


### Initialize list to hold results across NPS
fe_results_list.nps <- list()


### Run meta-analysis for each NPS
for(currnps in npsvars) {
    cat("Running for ", currnps, "\n")
    ### Make long table for current NPS
    res.list <- res.list.nps[[currnps]] # list with one element per cohort
    meta_long <- rbindlist(res.list, use.names = TRUE, fill = TRUE, idcol = "cohort")

    ## Keep only the columns needed for metafor
    meta_long <- meta_long[, .(protein, cohort, yi = logFC, sei = get(args$SE_type))]

    ## Filter to only the proteins that are present for all three cohorts
    proteincounts <- meta_long[, .N, by = "protein"]
    meta_long <- meta_long[protein %in% proteincounts[N == length(cohorts), protein], ]


    ### Run metafor and make table of stats
    fe_results <- meta_long[
      , {
          fit <- rma.uni(yi = yi, sei = sei, method = "FE")   # FE = fixed‑effects
          list(
              fe_beta   = coef(fit),                # pooled estimate
              fe_se     = sqrt(vcov(fit)),          # pooled SE
              fe_z      = summary(fit)$zval,        # Wald Z
              fe_p      = summary(fit)$pval,        # two‑sided p‑value
              fe_Q      = fit$QE,                   # heterogeneity Q (should be low),
              fe_Q_pval = pchisq(fit$QE, df = length(cohorts)-1, lower.tail = FALSE), # p-value of Cochran's Q
              fe_I2     = fit$I2                    # I² (should be ≈0 for FE)
          )
      },
      by = protein]


    ### Add multiple testing correction columns
    setorder(fe_results, fe_p)
    fe_results[, qvalue := qvalue(fe_p)$qvalue]
    fe_results[, p.BH := p.adjust(fe_p, method = "BH")]
    fe_results[, p.bonf := p.adjust(fe_p, method = "bonferroni")]
    fe_results[, sig.qvalue05 := ifelse(qvalue < 0.05, "sig", "")]


    ### Save result to cross-nps list
    fe_results_list.nps[[currnps]] <- fe_results
    rm(fe_results) ## clean up
}
rm(currnps)



###  Make qq plots
qq_list <- list()

for(currnps in npsvars) {
    pvals <- unlist(fe_results_list.nps[[currnps]][, fe_p])

    ## calculate inflation factor
    chisq.all <- qchisq(pvals, 1, lower.tail=FALSE)
    lambda.median.all <- median(chisq.all) / qchisq(0.5,1)
    lambdatext.all <- bquote(lambda ~ '=' ~ .(signif(lambda.median.all, 4)))

    ## draw plot
    n.na <- sum(is.na(pvals))
    titleqq <- paste("Fixed-effects cross-cohort meta-analysis \n",
                     names(npsvars)[npsvars == currnps])
    if(n.na > 0) { titleqq <- paste0(titleqq, "; excluding ", n.na, " proteins with NA P-value") }
    qq_list[[currnps]] <- qqunif.plot(pvals[!is.na(pvals)],
                                  maintitle=titleqq, subxlab=lambdatext.all,
                                  par.settings=list(par.sub.text=list(cex=1)))
}
rm(currnps)


### Make flat table for printing
## Can reshape as desired in a later step for supplementary table
## or followup analyses
fe_results_allnps <- rbindlist(fe_results_list.nps, idcol = "nps")


###### FINISH ######
## Write results table to file
fwrite(fe_results_allnps, file = outfiles$stats, quote = FALSE, sep = "\t", na = NA)

## Draw QQ plots to file
pdf(outfiles$qq)
invisible(sapply(qq_list, print))
dev.off()
