###############################################################################
## 03a_run_r2.r                                                        (Round 2)
##
## Round 2 version of 03a_code1_diagnostics.r. With the defaults it is the same analysis as
## 03a / Round 1 Code 1 (same filters, z-scaling, geeglm with unstructured working correlation,
## CR2 covariance, HTZ Wald test; P is identical) and records, for every protein, in addition
##   * the other small-sample tests of clubSandwich::Wald_test from the same CR2 fit
##   * chi-square Wald tests with the approximate jackknife covariance of geeglm
##     (--stderr_variants jack,j1s,fij): the "san.se and jack" runs of the earlier analysis
##   * the seven single-domain tests (CR2, Satterthwaite t)
## The output goes to a NEW folder (Prot_multivar_diag_R2); nothing earlier is overwritten.
##
## Usage (from the project directory, as for 03a):
##   Rscript 1_Code/03a_run_r2.r --cohort emory --residset ignoreSVs             # baseline of Round 2
##   Rscript 1_Code/03a_run_r2.r --cohort emory --residset noProtect             # existing residual file
##   Rscript 1_Code/03a_run_r2.r --cohort emory --residset regressSVs            # existing residual file
##   Rscript 1_Code/03a_run_r2.r --cohort emory --xtransform rankint             # Code 6
##   Rscript 1_Code/03a_run_r2.r --cohort emory --corstr independence            # Code 4
##   Rscript 1_Code/03a_run_r2.r --cohort emory --corstr exchangeable
## Smoke test first (20 proteins, a minute or two):
##   Rscript 1_Code/03a_run_r2.r --cohort emory --maxprot 20 --outtag R2smoke
###############################################################################
library(data.table)
library(argparse)
library(car)
library(geepack)
library(clubSandwich)

rm(list = ls())
has_qvalue <- requireNamespace("qvalue", quietly = TRUE)
if (has_qvalue) library(qvalue)
options(stringsAsFactors = FALSE)
options(nwarnings = 100000000)
source("1_Code/project_constants.r")
source("~/comsv/Svattathil_Library/svattathil_functions.r")
source("~/comsv/Svattathil_Library/manhattan_and_qq_plot_functions.r")
source("1_Code/diag_functions_r2.r")

###### SETUP ######
parser <- ArgumentParser()
parser$add_argument("--cohort", type = "character", default = "emory")
parser$add_argument("--residset", type = "character", default = "ignoreSVs",
                    help = "ignoreSVs (analysis), noProtect, regressSVs; or any label with --residfile")
parser$add_argument("--residfile", type = "character", default = NULL,
                    help = "path of a residual file; overrides Residfile(residset, cohort); --residset is then only a label")
parser$add_argument("--maxprot", type = "integer", default = 0, help = "run only the first N proteins (smoke test); 0 = all")
parser$add_argument("--corstr", type = "character", default = "unstructured",
                    help = "working correlation: unstructured (03a), exchangeable, independence")
parser$add_argument("--scalefix", action = "store_true", help = "fix the dispersion at 1 (geeglm scale.fix = TRUE)")
parser$add_argument("--xtransform", type = "character", default = "z",
                    help = "z (03a: z-scaling) or rankint (rank-based inverse-normal transform, Code 6)")
parser$add_argument("--stderr_variants", type = "character", default = "jack",
                    help = "comma-separated subset of jack,j1s,fij (refits; chi-square test with that jackknife covariance); none = skip")
parser$add_argument("--no_alltests", action = "store_true", help = "skip the other small-sample tests of Wald_test")
parser$add_argument("--no_perdomain", action = "store_true", help = "skip the seven single-domain tests")
parser$add_argument("--outtag", type = "character", default = NULL, help = "label of the output folder / files")
args <- parser$parse_args()
stopifnot(args$xtransform %in% c("z", "rankint"))
sev <- if (args$stderr_variants %in% c("none", "")) NULL else strsplit(args$stderr_variants, ",")[[1]]
stopifnot(all(sev %in% c("jack", "j1s", "fij")))
extra <- list(alltests = !args$no_alltests, stderr_variants = sev, perdomain = !args$no_perdomain)

if (is.null(args$outtag)) {
    args$outtag <- paste0("R2",
        ifelse(args$corstr != "unstructured", paste0("_", args$corstr), ""),
        ifelse(args$scalefix, "_scalefix", ""),
        ifelse(args$xtransform != "z", paste0("_", args$xtransform), ""))
}
is_default_model <- (args$corstr == "unstructured" && !args$scalefix && args$xtransform == "z")

outdir <- paste0("2_Pipeline/03-Run_regressions/Prot_multivar_diag_R2/", args$outtag, "_", Capwords(args$residset), "/")
MyMkdir(outdir)
fileid <- paste0("prot_multivar_", args$residset, "_", args$cohort, "_", args$outtag)
outfiles <- list(
    summstats = paste0(outdir, fileid, "_regression_stats_diag.txt"),
    cohort    = paste0(outdir, fileid, "_diag_cohort_summary.txt"),
    once      = paste0(outdir, fileid, "_diag_once.txt"),
    session   = paste0(outdir, fileid, "_diag_sessionInfo.txt"),
    qq        = paste0(outdir, fileid, "_qq.png"),
    diagqq    = paste0(outdir, fileid, "_diag_qq.png"),
    refqq     = paste0(outdir, fileid, "_refs_qq.png"),
    log       = paste0(outdir, fileid, ".log"))
## earlier summary files used for the reproduction checks (if present)
reffile_03a <- paste0("2_Pipeline/03-Run_regressions/Prot_multivar/Using_", Capwords(args$residset),
                      "/prot_multivar_", args$residset, "_", args$cohort, "_regression_stats.txt")
reffile_R1  <- paste0("2_Pipeline/03-Run_regressions/Prot_multivar_diag/Code1_", Capwords(args$residset),
                      "/prot_multivar_", args$residset, "_", args$cohort, "_Code1_regression_stats_diag.txt")


###### MAIN ######
dat <- PrepareData(args$cohort, args$residset, args$residfile, xtransform = args$xtransform)
transformed <- dat$transformed
phenos      <- dat$phenos
idx_run <- seq_len(nrow(transformed))
if (args$maxprot > 0) idx_run <- head(idx_run, args$maxprot)
set.seed(24983423)   # as in 03a

print("Starting regressions for each protein")
rows <- vector("list", length(idx_run))
for (j in seq_along(idx_run)) {
    i <- idx_run[j]
    if (j %% 100 == 0) print(paste0(pN(j), " of ", pN(length(idx_run))))
    res <- RunProteinAnalysis(x_vec = transformed[i, ], sample_ids = colnames(transformed),
                              phenos = phenos, npsvars = npsvars.bin,
                              corstr = args$corstr, scale_fix = args$scalefix, extra = extra)
    rows[[j]] <- ResultRow(rownames(transformed)[i], res, n_na_resid = sum(is.na(transformed[i, ])))
}
n.errors <- sum(!sapply(rows, function(x) is.na(x$error)))
if (n.errors > 0)
    print(paste0(pN(n.errors), " of ", pN(length(idx_run)), " proteins failed and were skipped (see 'error' column)."))

summstats <- rbindlist(lapply(rows, as.data.table), fill = TRUE)
setcolorder(summstats, "protein")

### q-values etc. (as 03a); qvalue() can fail with very few proteins, then the q-values are NA
summstats[!is.na(P), c("qvalue", "lfdr") := {
    qobj <- if (has_qvalue) tryCatch(qvalue(P), error = function(e) NULL) else NULL
    if (is.null(qobj)) list(rep(NA_real_, .N), rep(NA_real_, .N)) else list(qobj$qvalues, qobj$lfdr)
}]
summstats[, P.BH := p.adjust(P, method = "BH")]
summstats[, P.bonf := p.adjust(P, method = "bonferroni")]
summstats[, sig.qvalue05 := ifelse(qvalue < 0.05, TRUE, FALSE)]
summstats[, sig := ifelse(qvalue < 0.05, "sig", NA_character_)]
setorder(summstats, P)
summstats[, p.rank := rank(P, ties.method = "max")]
setcolorder(summstats, "p.rank")

### Q-Q plot as in 03a
n.na <- sum(is.na(summstats$P))
titleqq <- paste0(fileid, "\n", "N=", nrow(phenos))
if (n.na > 0) titleqq <- paste0(titleqq, "; excluding ", n.na, " proteins with NA P-value")
forqq <- summstats[!is.na(P), P]; forqq[which(forqq == 0)] <- 1e-300
qqplot <- qqunif.plot(forqq, maintitle = titleqq)
png(outfiles$qq, height = 7, width = 7, units = "in", res = 300); print(qqplot); dev.off()


###### DIAGNOSTICS ######
cr_methods <- tryCatch(as.character(methods("vcovCR")), error = function(e) NA_character_)
once <- c(cohort = args$cohort, residset = args$residset, corstr = args$corstr, scalefix = args$scalefix,
          xtransform = args$xtransform, stderr_variants = paste(sev, collapse = ","),
          alltests = extra$alltests, perdomain = extra$perdomain,
          phenos_file = dat$infiles$phenos, resid_file = dat$infiles$resid, dat$counts,
          n_proteins_run = length(idx_run), n_failed_fits = n.errors,
          fit_class_first = na.omit(summstats$fit_class)[1],
          vcovCR_has_geeglm_method = !is.null(getS3method("vcovCR", "geeglm", optional = TRUE)),
          vcovCR_methods = paste(cr_methods, collapse = ", "),
          R = R.version.string, geepack = as.character(packageVersion("geepack")),
          clubSandwich = as.character(packageVersion("clubSandwich")),
          qvalue = if (has_qvalue) as.character(packageVersion("qvalue")) else NA,
          data.table = as.character(packageVersion("data.table")))
fwrite(data.table(item = names(once), value = as.character(once)), file = outfiles$once, sep = "\t", quote = FALSE)
writeLines(capture.output(print(sessionInfo())), outfiles$session)
if (n.errors > 0) {
    errtab <- summstats[!is.na(error), .N, by = error][order(-N)]
    fwrite(errtab, file = sub("_diag_once.txt", "_diag_errors.txt", outfiles$once), sep = "\t", quote = FALSE)
}

cs <- CohortSummary(summstats)

### internal consistency of the Round 2 columns (should be ~0)
lp <- function(p) log10(pmax(p, 1e-300))
chk <- list()
if ("p_HTZ" %in% names(summstats))
    chk$check_max_diff_log10P_vs_pHTZ <- max(abs(lp(summstats$P) - lp(summstats$p_HTZ)), na.rm = TRUE)
if ("p_chisq" %in% names(summstats))
    chk$check_max_diff_log10P_chisq_vs_cr2 <- max(abs(lp(summstats$p_chisq) - lp(summstats$p_chisq_cr2)), na.rm = TRUE)
if (length(chk)) cs <- rbind(cs, data.table(metric = names(chk), value = unlist(chk)))

### reproduction checks against the earlier analyses of the same cohort (default model only)
repro <- function(ref_p, label) {
    mm <- merge(summstats[, .(protein, P)], ref_p, by = "protein")
    d  <- abs(lp(mm$P) - lp(mm$P_ref))
    print(paste0("Reproduction check against ", label, ": max |diff log10 P| = ", signif(max(d, na.rm = TRUE), 3),
                 " over ", nrow(mm), " proteins"))
    data.table(metric = paste0(c("repro_n_compared_", "repro_max_abs_diff_log10P_", "repro_n_diff_gt_1e-6_"), label),
               value = c(nrow(mm), max(d, na.rm = TRUE), sum(d > 1e-6, na.rm = TRUE)))
}
if (is_default_model && is.null(args$residfile)) {
    if (file.exists(reffile_03a)) cs <- rbind(cs, repro(fread(reffile_03a)[, .(protein, P_ref = P)], "03a"))
    if (file.exists(reffile_R1))  cs <- rbind(cs, repro(fread(reffile_R1)[, .(protein, P_ref = P)], "Round1"))
} else print("Reproduction check skipped (non-default model, own residual file or other transform)")
fwrite(cs, file = outfiles$cohort, sep = "\t", quote = FALSE)
print(cs[grepl("lambda_GC|frac_p_lt_0.05|frac_d2|repro|check_", metric)])

DiagQQPlots(summstats, outfiles$diagqq, title = fileid)
RefQQPlot(summstats, outfiles$refqq, title = fileid)


###### FINISH ######
write.table(summstats, file = outfiles$summstats, row.names = FALSE, quote = FALSE, sep = "\t")
logtab <- data.table(cohort = args$cohort,
                     model = paste0("outcome ~ nps_id + protein:nps_id (geeglm; corstr=", args$corstr,
                                    "; scale.fix=", args$scalefix, "; x=", args$xtransform, "; CR2; HTZ)"),
                     n.proteins.tested = length(idx_run),
                     n.sig.qvalue05 = sum(summstats$sig.qvalue05, na.rm = TRUE))
write("\nModel and significance summary:", file = outfiles$log)
fwrite(logtab, file = outfiles$log, sep = "\t", quote = FALSE, col.names = TRUE, append = TRUE)
print(paste0("Done. Output in ", outdir))
