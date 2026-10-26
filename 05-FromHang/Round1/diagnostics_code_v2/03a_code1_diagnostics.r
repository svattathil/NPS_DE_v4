###############################################################################
## 03a_code1_diagnostics.r   (Code 1; Code 2 with --residset regressSVs;
##                            Code 3 with --residfile; Code 4 with --corstr / --scalefix)
##
## Same analysis as 03a-run_multivariate.r (same filters, z-scaling, geeglm with
## unstructured working correlation, CR2 covariance, HTZ Wald test), plus
## diagnostics recorded for every protein and every cohort.
## The output goes to a NEW folder; nothing of the earlier results is overwritten.
##
## Usage (from the project directory, as for 03a):
##   Rscript 1_Code/03a_code1_diagnostics.r --cohort Emory --residset ignoreSVs
##   Rscript 1_Code/03a_code1_diagnostics.r --cohort Emory --residset regressSVs   # Code 2
##
## Code 4 (working-model sensitivity; each is a separate run):
##   Rscript 1_Code/03a_code1_diagnostics.r --cohort Emory --corstr independence
##   Rscript 1_Code/03a_code1_diagnostics.r --cohort Emory --corstr exchangeable
##   Rscript 1_Code/03a_code1_diagnostics.r --cohort Emory --scalefix
##   Rscript 1_Code/03a_code1_diagnostics.r --cohort Emory --corstr independence --scalefix
##   Rscript 1_Code/03a_code1_diagnostics.r --cohort Emory --corstr exchangeable --scalefix
##
## Code 3 (residual files made by 01a_code3_cleaning_variants.r), e.g.
##   Rscript 1_Code/03a_code1_diagnostics.r --cohort Emory --residset noProtect
##   Rscript 1_Code/03a_code1_diagnostics.r --cohort Emory --residset naAware \
##       --residfile 2_Pipeline/02-Prepare_analysis_data/resid_regress_BatchPmiAgeMsex_naAware_Emory.txt
##
## Smoke test first (20 proteins, a few minutes):
##   Rscript 1_Code/03a_code1_diagnostics.r --cohort Emory --maxprot 20
###############################################################################
library(data.table)
library(argparse)
library(qvalue)
library(car)
library(geepack)
library(clubSandwich)

rm(list=ls())

### CONSTANTS ###
options(stringsAsFactors = FALSE)
options(nwarnings = 100000000)
source("1_Code/project_constants.r")

###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")
source("~/comsv/Svattathil_Library/manhattan_and_qq_plot_functions.r")
source("1_Code/diag_functions.r")


###### SETUP ######
parser <- ArgumentParser()
parser$add_argument("--cohort", type="character", default = "OHSU")
parser$add_argument("--residset", type="character",
                    help="regressSVs or ignoreSVs or resid2_multi or noProtect", default = "ignoreSVs")
parser$add_argument("--maxprot", type="integer", default = 0,
                    help="run only the first N proteins (smoke test); 0 = all")
parser$add_argument("--residfile", type="character", default = NULL,
                    help="path of a residual file; overrides Residfile(residset, cohort). --residset is then only a label")
parser$add_argument("--corstr", type="character", default = "unstructured",
                    help="working correlation: unstructured (03a), exchangeable, independence")
parser$add_argument("--scalefix", action = "store_true",
                    help="fix the dispersion at 1 (geeglm scale.fix = TRUE)")
parser$add_argument("--outtag", type="character", default = NULL,
                    help="label of the output folder / files (default: Code1, or Code4_<corstr>[_scalefix])")
args <- parser$parse_args()
if (is.null(args$outtag)) {
    args$outtag <- if (args$corstr == "unstructured" && !args$scalefix) "Code1" else
        paste0("Code4_", args$corstr, ifelse(args$scalefix, "_scalefix", ""))
}
is_default_model <- (args$corstr == "unstructured" && !args$scalefix)

### Folders and files (new folder, so no earlier result is overwritten)
outdir <- paste0("2_Pipeline/03-Run_regressions/Prot_multivar_diag/",
                 args$outtag, "_", Capwords(args$residset), "/")
MyMkdir(outdir)
fileid <- paste0("prot_multivar_", args$residset, "_", args$cohort, "_", args$outtag)
outfiles <- list(
    summstats = paste0(outdir, fileid, "_regression_stats_diag.txt"),
    cohort    = paste0(outdir, fileid, "_diag_cohort_summary.txt"),
    once      = paste0(outdir, fileid, "_diag_once.txt"),
    session   = paste0(outdir, fileid, "_diag_sessionInfo.txt"),
    qq        = paste0(outdir, fileid, "_qq.png"),
    diagqq    = paste0(outdir, fileid, "_diag_qq.png"),
    log       = paste0(outdir, fileid, ".log")
)
## the summary table of the earlier 03a run of the same cohort and option, if present
reffile <- paste0("2_Pipeline/03-Run_regressions/Prot_multivar/Using_", Capwords(args$residset),
                  "/prot_multivar_", args$residset, "_", args$cohort, "_regression_stats.txt")


###### MAIN ######
dat <- PrepareData(args$cohort, args$residset, args$residfile)
transformed <- dat$transformed
phenos      <- dat$phenos

idx_run <- seq_len(nrow(transformed))
if (args$maxprot > 0) idx_run <- head(idx_run, args$maxprot)

set.seed(24983423)   # as in 03a

print("Starting regressions for each protein")
rows <- vector("list", length(idx_run))
for (j in seq_along(idx_run)) {
    i <- idx_run[j]
    if (j %% 100 == 0) { print(paste0(pN(j), " of ", pN(length(idx_run)))) }
    res <- RunProteinAnalysis(x_vec = transformed[i, ], sample_ids = colnames(transformed),
                              phenos = phenos, npsvars = npsvars.bin,
                              corstr = args$corstr, scale_fix = args$scalefix)
    rows[[j]] <- ResultRow(rownames(transformed)[i], res,
                           n_na_resid = sum(is.na(transformed[i, ])))
}

n.errors <- sum(!sapply(rows, function(x) is.na(x$error)))
if (n.errors > 0) {
    print(paste0(pN(n.errors), " of ", pN(length(idx_run)),
                 " proteins failed and were skipped (see 'error' column in output)."))
}

summstats <- rbindlist(lapply(rows, as.data.table), fill = TRUE)
setcolorder(summstats, "protein")


### qvalue etc. (03a L252-262)
## (qvalue() can fail with very few proteins, e.g. in a smoke test; then q-values are NA)
summstats[!is.na(P), c("qvalue", "lfdr") := {
    qobj <- tryCatch(qvalue(P), error = function(e) NULL)
    if (is.null(qobj)) list(rep(NA_real_, .N), rep(NA_real_, .N)) else list(qobj$qvalues, qobj$lfdr)
}]
summstats[, P.BH := p.adjust(P, method="BH")]
summstats[, P.bonf := p.adjust(P, method="bonferroni")]
summstats[, sig.qvalue05 := ifelse(qvalue < 0.05, TRUE, FALSE)]
summstats[, sig := ifelse(qvalue < 0.05, "sig", NA_character_)]
setorder(summstats, P)
summstats[, p.rank := rank(P, ties.method="max")]
setcolorder(summstats, "p.rank")


### Q-Q plot as in 03a (L274-281, L308-311)
n.na <- sum(is.na(summstats$P))
titleqq <- paste0(fileid, "\n", "N=", nrow(phenos))
if(n.na > 0) { titleqq <- paste0(titleqq, "; excluding ", n.na, " proteins with NA P-value") }
forqq <- summstats[!is.na(P), P]
forqq[which(forqq == 0)] <- 1e-300
qqplot <- qqunif.plot(forqq, maintitle=titleqq)
png(outfiles$qq, height=7, width=7, units="in", res=300)
print(qqplot)
dev.off()


###### DIAGNOSTICS ######
### (a) once per run: software and sample bookkeeping
cr_methods <- tryCatch(as.character(methods("vcovCR")), error = function(e) NA_character_)
once <- c(
    cohort = args$cohort, residset = args$residset, corstr = args$corstr,
    scalefix = args$scalefix,
    phenos_file = dat$infiles$phenos, resid_file = dat$infiles$resid,
    dat$counts,
    n_proteins_run = length(idx_run),
    n_failed_fits = n.errors,
    fit_class_first = na.omit(summstats$fit_class)[1],
    vcovCR_has_geeglm_method = !is.null(getS3method("vcovCR", "geeglm", optional = TRUE)),
    vcovCR_methods = paste(cr_methods, collapse = ", "),
    R = R.version.string,
    geepack = as.character(packageVersion("geepack")),
    clubSandwich = as.character(packageVersion("clubSandwich")),
    qvalue = as.character(packageVersion("qvalue")),
    data.table = as.character(packageVersion("data.table")))
fwrite(data.table(item = names(once), value = as.character(once)),
       file = outfiles$once, sep = "\t", quote = FALSE)
writeLines(capture.output(print(sessionInfo())), outfiles$session)

### the failed fits tabulated by error message
if (n.errors > 0) {
    errtab <- summstats[!is.na(error), .N, by = error][order(-N)]
    fwrite(errtab, file = sub("_diag_once.txt", "_diag_errors.txt", outfiles$once),
           sep = "\t", quote = FALSE)
}

### (f) per cohort: summaries
cs <- CohortSummary(summstats)

### reproduction check against the earlier 03a output (same cohort, same option)
if (is_default_model && is.null(args$residfile) && file.exists(reffile)) {
    ref <- fread(reffile)[, .(protein, P_ref = P)]
    mm  <- merge(summstats[, .(protein, P)], ref, by = "protein")
    lp  <- function(p) log10(pmax(p, 1e-300))
    d   <- abs(lp(mm$P) - lp(mm$P_ref))
    cs <- rbind(cs, data.table(
        metric = c("repro_n_compared", "repro_max_abs_diff_log10P", "repro_n_diff_gt_1e-6"),
        value  = c(nrow(mm), max(d, na.rm = TRUE), sum(d > 1e-6, na.rm = TRUE))))
    print(paste0("Reproduction check against ", reffile, ": max |diff log10 P| = ",
                 signif(max(d, na.rm = TRUE), 3), " over ", nrow(mm), " proteins"))
} else {
    print("Reproduction check skipped (non-default model, own residual file, or no earlier summary file)")
}
fwrite(cs, file = outfiles$cohort, sep = "\t", quote = FALSE)
print(cs[grepl("lambda|frac_d2|spearman|repro|max_rel|delta", metric)])

DiagQQPlots(summstats, outfiles$diagqq, title = fileid)


###### FINISH ######
write.table(summstats, file = outfiles$summstats, row.names = FALSE, quote = FALSE, sep = "\t")
logtab <- data.table(cohort            = args$cohort,
                     model             = paste0("outcome ~ nps_id + protein:nps_id (geeglm; corstr=", args$corstr,
                                        "; scale.fix=", args$scalefix, "; CR2; HTZ)"),
                     n.proteins.tested = length(idx_run),
                     n.sig.qvalue05    = sum(summstats$sig.qvalue05, na.rm=TRUE))
write("\nModel and significance summary:", file = outfiles$log)
fwrite(logtab, file = outfiles$log, sep = "\t", quote = FALSE, col.names = TRUE, append = TRUE)
print(paste0("Done. Output in ", outdir))
