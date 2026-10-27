library(data.table)
library(argparse)
library(stringr)
library(qvalue)

rm(list = ls())


###### GOAL ######
## Test association of consensus module eigengenes with NPS domains,
## using eigengenes already written by 03-build_consensus_network.r
## Eigengene is the outcome in all models.
##   (1) each NPS separately (linear regression):
##         mod.eigen ~ NPS + covariates
##   (2) all NPS jointly (partial F-test):
##         mod.eigen ~ covariates + all NPS  vs  mod.eigen ~ covariates
##       then combine per-cohort F-test p-values across cohorts (Stouffer)
## Covariates are set with --covars; use --covars none (or --covars with no
## values) to fit models with the NPS variables as the only predictors.


##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")
source("~/comsv/Svattathil_Library/svattathil_functions.r")
source("~/comsv/Svattathil_Library/manhattan_and_qq_plot_functions.r")

nps <- npsvars.bin


##### COMMAND LINE ARGUMENTS #####
parser <- ArgumentParser()
parser$add_argument("--auto_softpower", type="integer", default=9, help="")
parser$add_argument("--use_powervec1", action="store_true", default=FALSE, help="")
parser$add_argument("--residset", type="character", default="ignoreSVs", help="")
parser$add_argument("--covars", type="character", nargs="*", default="",
                    help=paste0("Covariates for all association models. ",
                                "Use 'none' or give no values for NPS-only models."))
parser$add_argument("--stouffer_weights", type="character", default="sqrtN",
                    choices=c("sqrtN", "equal"),
                    help="Weights for Stouffer combination of per-cohort F-test p-values (default: sqrtN)")

args <- parser$parse_args()

if(args$use_powervec1) {
    network_group <- "powervec1"
} else {
    network_group <- paste0("auto_softpower", args$auto_softpower)
}

## Covariates: empty vector means NPS-only models
covars <- setdiff(as.character(unlist(args$covars)), c("none", ""))
print(paste("Covariates:", if(length(covars) > 0) paste(covars, collapse = ", ") else "none"))

outtag <- if(length(covars) > 0) paste0("covar", sapply(covars, Capwords)) else ""


###### FUNCTIONS ######

PrepareDataME <- function(modulename, aphenodata, amedata) {
    ## Extract eigengene for current module and add protsample as column
    mes.to.merge <- data.frame(protsample=rownames(amedata),
                               mod.eigen=amedata[,modulename])

    ## Add eigengene column to phenos
    merged <- merge(mes.to.merge, aphenodata, by = "protsample")

    ## Convert to data frame and add protsample as rownames
    samps <- merged$protsample
    setDF(merged)
    rownames(merged) <- samps
    return(merged)
}

MakeFormula <- function(predictors, response = "mod.eigen") {
    ## reformulate() fails on an empty predictor vector; use intercept-only instead
    if(length(predictors) == 0) predictors <- "1"
    reformulate(predictors, response = response)
}

RunPerNPSTest <- function(amodule, anps, aphenodata, amedata, covars) {
    ## Linear regression of eigengene on one NPS (plus covariates, if any)
    ## mod.eigen ~ NPS + covars
    regression_data <- PrepareDataME(amodule, aphenodata, amedata)
    keepvars <- c("mod.eigen", anps, covars)
    regression_data <- regression_data[complete.cases(regression_data[, keepvars]), keepvars, drop = FALSE]

    fit <- lm(MakeFormula(c(anps, covars)), data = regression_data)
    coefs <- summary(fit)$coefficients

    ## NPS coefficient is dropped (aliased) if the NPS has no variation in this cohort
    if(anps %in% rownames(coefs)) {
        res <- coefs[anps, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
    } else {
        warning("NPS ", anps, " not estimable for ", amodule, " (no variation?); returning NA")
        res <- rep(NA_real_, 4)
    }

    data.table(eigengene = amodule,
               Estimate  = res[1],
               Std.Err   = res[2],
               Tstat     = res[3],
               N         = nrow(regression_data),
               Pvalue    = res[4])
}

RunJointNPSTest <- function(amodule, aphenodata, amedata, nps, covars) {
    ## Partial F-test of all NPS domains jointly for one module
    ## Full:    mod.eigen ~ covars + NPS_1 + ... + NPS_k
    ## Reduced: mod.eigen ~ covars   (intercept-only if no covariates)
    ## H0: all NPS coefficients = 0 (conditional on covariates)
    regression_data <- PrepareDataME(amodule, aphenodata, amedata)

    ## Restrict to complete cases across eigengene, covariates, and ALL NPS
    ## so both models are fit to identical samples (required for the F-test)
    keepvars <- c("mod.eigen", covars, nps)
    regression_data <- regression_data[complete.cases(regression_data[, keepvars]), keepvars, drop = FALSE]

    fit_reduced <- lm(MakeFormula(covars), data = regression_data)
    fit_full    <- lm(MakeFormula(c(covars, nps)), data = regression_data)
    ftest <- anova(fit_reduced, fit_full)

    data.table(eigengene  = amodule,
               Fstat      = ftest[2, "F"],
               df_num     = ftest[2, "Df"],        ## taken from anova, so robust to aliased terms
               df_den     = ftest[2, "Res.Df"],
               N          = nrow(regression_data),
               partial_R2 = (ftest[1, "RSS"] - ftest[2, "RSS"]) / ftest[1, "RSS"],
               Pvalue     = ftest[2, "Pr(>F)"])
}

StoufferCombine <- function(pvals, N, weighting = "sqrtN") {
    ## Combine one-sided p-values (F-test p-values are upper-tail by construction)
    ## Z_i = qnorm(1 - p_i);  Z = sum(w_i * Z_i) / sqrt(sum(w_i^2))
    keep <- !is.na(pvals)
    pvals <- pvals[keep]; N <- N[keep]
    if(length(pvals) == 0) return(list(Z = NA_real_, Pvalue = NA_real_, n_cohorts = 0L, N_total = 0L))

    ## Keep p strictly inside (0, 1) so Z is finite
    pvals <- pmin(pmax(pvals, .Machine$double.xmin), 1 - .Machine$double.eps)
    z <- qnorm(pvals, lower.tail = FALSE)
    w <- if(weighting == "sqrtN") sqrt(N) else rep(1, length(z))
    Z <- sum(w * z) / sqrt(sum(w^2))

    list(Z = Z, Pvalue = pnorm(Z, lower.tail = FALSE),
         n_cohorts = length(pvals), N_total = sum(N))
}


###### SETUP ######
### Files that exist
netdir <- paste0("2_Pipeline/04-Build_networks/",
                 paste0("Using_", args$residset), "/",
                 "Consensus/Auto_softpower", args$auto_softpower, "/")

if(!dir.exists(netdir)) {
    stop("Network directory not found: ", netdir,
         "\nRun 03-build_consensus_network.r with matching arguments first.") }

Netfilepref <- paste0(netdir, network_group)

infiles <- list()
infiles[["eigengenes"]] <- sapply(cohorts, function(acohort) {
    paste0(Netfilepref, "_", acohort, "_eigengenes.txt") })

infiles[["phenos"]] <- sapply(cohorts, function(acohort) {
    paste0("2_Pipeline/02-Prepare_analysis_data/phenos_cleaned_SVs_", acohort, ".txt") })

missing_files <- unlist(infiles)[!file.exists(unlist(infiles))]
if(length(missing_files) > 0) {
    stop("Missing input files:\n", paste(missing_files, collapse = "\n")) }


### Files to be created
outfiles <- list()
outfiles[["regression"]]     <- paste0(Netfilepref, outtag, "_nps_linear_regression_results_per-cohort.txt")
outfiles[["qq_per-domain"]]  <- paste0(Netfilepref, outtag, "_nps_linear_regression_results_per-cohort_qqplots.pdf")
outfiles[["joint_Ftest"]]    <- paste0(Netfilepref, outtag, "_nps_joint_partialF_results_per-cohort.txt")
outfiles[["qq_joint_Ftest"]] <- paste0(Netfilepref, outtag, "_nps_joint_partialF_results_per-cohort_qqplots.pdf")
outfiles[["stouffer"]]       <- paste0(Netfilepref, outtag, "_nps_joint_partialF_results_stouffer.txt")
outfiles[["qq_stouffer"]]    <- paste0(Netfilepref, outtag, "_nps_joint_partialF_results_stouffer_qqplot.pdf")


###### MAIN ######
### Read in data ###
## Module eigengenes: one data frame per cohort, protsample as rownames
## (same structure as MEs_per_cohort in 03-build_consensus_network.r)
MEs_per_cohort <- sapply(cohorts, function(acohort) {
    tmp <- fread(infiles[["eigengenes"]][[acohort]], colClasses = list(character = "protsample"))
    tmp_sampnames <- tmp$protsample
    tmp[, protsample := NULL]
    setDF(tmp)
    rownames(tmp) <- tmp_sampnames
    return(tmp)
}, simplify = FALSE)

## Phenotype data
phenos_list <- sapply(cohorts, function(acohort) {
    fread(infiles[["phenos"]][[acohort]], colClasses = list(character = "protsample"))
}, simplify = FALSE)

## Check required phenotype columns, and that NPS are numeric
## (so the coefficient is named after the NPS variable in the lm output)
for(acohort in cohorts) {
    missing_vars <- setdiff(c("protsample", covars, nps), names(phenos_list[[acohort]]))
    if(length(missing_vars) > 0) {
        stop("Cohort ", acohort, " phenotype file missing: ", paste(missing_vars, collapse = ", ")) }

    for(anps in nps) {
        if(is.logical(phenos_list[[acohort]][[anps]])) {
            phenos_list[[acohort]][, (anps) := as.integer(get(anps))] }
        if(!is.numeric(phenos_list[[acohort]][[anps]])) {
            stop("Cohort ", acohort, ": NPS variable ", anps, " must be numeric (0/1)") }
    }
}


### Test association between each NPS and module eigengenes (linear regression)
## Nested list: one list per cohort, each a list of results per NPS
allsumms_list <- sapply(cohorts, function(acohort) {
    MEs_thiscohort <- MEs_per_cohort[[acohort]]
    testmodules <- setdiff(colnames(MEs_thiscohort), "MEgrey")   ## exclude grey

    res <- sapply(nps, function(anps) {
        out <- rbindlist(lapply(testmodules, RunPerNPSTest,
                                anps       = anps,
                                aphenodata = phenos_list[[acohort]],
                                amedata    = MEs_thiscohort,
                                covars     = covars))
        print(paste(acohort, "--", "model", anps, "-- done"))
        return(out)
    }, simplify = FALSE)
    return(res)
}, simplify = FALSE)

## Flatten: across NPS within cohort
allsumms_list2 <- sapply(allsumms_list, function(x) { rbindlist(x, idcol = "NPS") }, simplify = FALSE)

## Add qvalues
for(acohort in cohorts) {
    allsumms_list2[[acohort]][Pvalue < 0.05, sig.p05 := "sig_nominal05"]
    allsumms_list2[[acohort]][, qvalue := qvalue(Pvalue)$qvalue]
    allsumms_list2[[acohort]][qvalue < 0.05, sig.qvalue05 := "sig_qvalue05"]
    }
rm(acohort)


##  Draw qq plots
qq_perdomain <- sapply(cohorts, function(acohort) {
    ## calculate inflation factor
    chisq.all <- qchisq(allsumms_list2[[acohort]][, Pvalue], 1, lower.tail=FALSE)
    lambda.median.all <- median(chisq.all) / qchisq(0.5,1)
    lambdatext.all <- bquote(lambda ~ '=' ~ .(signif(lambda.median.all, 4)))

    ## Draw plot
    n.na <- sum(is.na(allsumms_list2[[acohort]]$Pvalue))
    titleqq <- paste0(acohort, " - per-domain tests \n (one point per module*BPSD)")
    qq_stouffer <- qqunif.plot(allsumms_list2[[acohort]][!is.na(Pvalue), Pvalue],
                               maintitle=titleqq, subxlab=lambdatext.all,
                               par.settings=list(par.sub.text=list(cex=1)))
}, simplify = FALSE)


## Flatten: across cohorts (long format)
allsumms_long <- rbindlist(allsumms_list2, idcol = "cohort")
setorder(allsumms_long, eigengene, NPS)
setcolorder(allsumms_long, c("eigengene", "NPS"))
setorder(allsumms_long, Pvalue)


### Test joint association between all NPS and module eigengenes (partial F-test)
jointsumms_list <- sapply(cohorts, function(acohort) {
    MEs_thiscohort <- MEs_per_cohort[[acohort]]
    testmodules <- setdiff(colnames(MEs_thiscohort), "MEgrey")
    rbindlist(lapply(testmodules, RunJointNPSTest,
                     aphenodata = phenos_list[[acohort]],
                     amedata    = MEs_thiscohort,
                     nps        = nps,
                     covars     = covars))
}, simplify = FALSE)

## Add multiple testing correction and sig columns (within each cohort)
for(acohort in cohorts) {
    jointsumms_list[[acohort]][, p.BH := p.adjust(Pvalue, method = "BH")]
    jointsumms_list[[acohort]][, p.bonf := p.adjust(Pvalue, method = "bonferroni")]

    jointsumms_list[[acohort]][Pvalue < 0.05, sig_nominal05 := "sig_nominal05"]
    jointsumms_list[[acohort]][p.BH < 0.05, sig_BH05 := "sig_BH05"]
    jointsumms_list[[acohort]][p.bonf < 0.05, sig_bonf05 := "sig_bonf05"]
}
rm(acohort)


##  Draw qq plots
qq_joint_Ftest <- sapply(cohorts, function(acohort) {

    ## calculate inflation factor
    chisq.all <- qchisq(jointsumms_list[[acohort]][, Pvalue], 1, lower.tail=FALSE)
    lambda.median.all <- median(chisq.all) / qchisq(0.5,1)
    lambdatext.all <- bquote(lambda ~ '=' ~ .(signif(lambda.median.all, 4)))

    ## Draw plot
    n.na <- sum(is.na(jointsumms_list[[acohort]]$Pvalue))
    titleqq <- paste0(acohort, " - partial F-test \n (one point per module)")
    qq_stouffer <- qqunif.plot(jointsumms_list[[acohort]][!is.na(Pvalue), Pvalue],
                               maintitle=titleqq, subxlab=lambdatext.all,
                               par.settings=list(par.sub.text=list(cex=1)))
}, simplify = FALSE)


## Combine into single table (aggregate across cohorts)
jointsumms_long <- rbindlist(jointsumms_list, idcol = "cohort")
jointsumms_long[, NPS := "all_NPS"]
setorder(jointsumms_long, Pvalue, eigengene, cohort)
setcolorder(jointsumms_long, c("eigengene", "NPS", "cohort"))


### Combine per-cohort F-test p-values across cohorts (Stouffer)
stouffer_summ <- jointsumms_long[, {
    sc <- StoufferCombine(Pvalue, N, weighting = args$stouffer_weights)
    list(n_cohorts = sc$n_cohorts,
         N_total   = sc$N_total,
         weighting = args$stouffer_weights,
         Stouffer_Z = sc$Z,
         Pvalue     = sc$Pvalue)
}, by = .(eigengene, NPS)]
setorder(stouffer_summ, Pvalue)

## Add multiple testing correction and sig columns
stouffer_summ[Pvalue < 0.05, sig_p05 := "sig_nominal05"]
stouffer_summ[, p.BH := p.adjust(Pvalue, method = "BH")]
stouffer_summ[, p.bonf := p.adjust(Pvalue, method = "bonferroni")]
stouffer_summ[p.BH < 0.05, sig_BH05 := "sig_BH05"]
stouffer_summ[p.bonf < 0.05, sig_bonf05 := "sig_bonf05"]


## calculate inflation factor
chisq.all <- qchisq(stouffer_summ[, Pvalue], 1, lower.tail=FALSE)
lambda.median.all <- median(chisq.all) / qchisq(0.5,1)
lambdatext.all <- bquote(lambda ~ '=' ~ .(signif(lambda.median.all, 4)))


##  Draw qq plot
n.na <- sum(is.na(stouffer_summ$Pvalue))
titleqq <- paste0("F-test cross-cohort combined p-value (Stouffer)")
qq_stouffer <- qqunif.plot(stouffer_summ[!is.na(Pvalue), Pvalue],
                   maintitle=titleqq, subxlab=lambdatext.all,
                   par.settings=list(par.sub.text=list(cex=1)))


### FINISH ###
write.table(allsumms_long,   file = outfiles[["regression"]],
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(jointsumms_long, file = outfiles[["joint_Ftest"]],
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(stouffer_summ,   file = outfiles[["stouffer"]],
            row.names = FALSE, quote = FALSE, sep = "\t")

pdf(outfiles[["qq_per-domain"]], height = 7, width = 7)
for(x in qq_perdomain) { print(x) }
dev.off()


pdf(outfiles[["qq_joint_Ftest"]], height = 7, width = 7)
for(x in qq_joint_Ftest) { print(x) }
dev.off()

pdf(outfiles[["qq_stouffer"]], height = 7, width = 7)
print(qq_stouffer)
dev.off()

