###############################################################################
## 03a_perm_null_r2.r                    (Round 2: permutation nulls 5a / 5b / 5c)
##
## Under a permutation of the NPS vectors among subjects every protein-NPS null hypothesis is true,
## so the p-values of the permuted data show the calibration of the test. Three designs:
##
##   --design fixedx  (5a, Round 1) permute the NPS vectors, keep the cleaned expression x FIXED.
##   --design all     (5b) permute the NPS vectors among all subjects AND repeat the cleaning of
##                    script V with the permuted NPS columns, then z-scale and test.
##   --design batch   (5c) as 5b, but the NPS vectors are permuted only within Batch (keeps the
##                    batch-NPS association of the real data under the null).
##
## Why 5b: with --protect yes (the analysis) the cleaning keeps the fit of the 7 NPS outcomes in x
## (8 protected columns, about 66 residual df). In real data this can make x slightly more related to
## the outcomes than chance; 5a cannot see that, 5b can. Expected: 5b null reproduces the real-data
## shift of Q, and --protect no (option noProtect) removes it.
##
##   --protect yes|no     NPS columns protected in the cleaning (yes = ignoreSVs, no = noProtect)
##   --na zero|aware      missing values set to 0 when estimating the coefficients (as jaffelab) or
##                        estimated from observed samples only
## Replicate 0 is the REAL data processed the same way (same cleaning, scaling, test).
##
## Usage (from the project directory):
##   Rscript 1_Code/03a_perm_null_r2.r --cohort emory --design fixedx --B 20 --nprot 500   # = Round 1 5a
##   Rscript 1_Code/03a_perm_null_r2.r --cohort emory --design all   --protect yes --B 20 --nprot 500
##   Rscript 1_Code/03a_perm_null_r2.r --cohort emory --design all   --protect no  --B 20 --nprot 500
##   Rscript 1_Code/03a_perm_null_r2.r --cohort emory --design batch --protect yes --B 20 --nprot 500
## Smoke test first:  --B 2 --nprot 30
##
## Each replicate is written to its own file; if the job is stopped, running the same command again
## continues. A manifest records the settings; if the settings or the protein subset differ from the
## existing replicate files, those files are MOVED to a stale_<time>/ sub-folder, never reused.
###############################################################################
library(data.table)
library(argparse)
library(car)
library(geepack)
library(clubSandwich)

rm(list = ls())
options(stringsAsFactors = FALSE)
options(nwarnings = 100000000)
source("1_Code/project_constants.r")
source("~/comsv/Svattathil_Library/svattathil_functions.r")
source("~/comsv/Svattathil_Library/manhattan_and_qq_plot_functions.r")
source("1_Code/diag_functions_r2.r")

parser <- ArgumentParser()
parser$add_argument("--cohort", type = "character", default = "emory")
parser$add_argument("--design", type = "character", default = "all", help = "fixedx (5a), all (5b), batch (5c)")
parser$add_argument("--protect", type = "character", default = "yes", help = "yes (ignoreSVs) or no (noProtect); designs all/batch")
parser$add_argument("--na", type = "character", default = "zero", help = "zero (jaffelab) or aware; designs all/batch")
parser$add_argument("--B", type = "integer", default = 20, help = "number of permutations")
parser$add_argument("--nprot", type = "integer", default = 500, help = "random subset of proteins (0 = all)")
parser$add_argument("--seed", type = "integer", default = 20261011)
parser$add_argument("--residset", type = "character", default = "ignoreSVs", help = "design fixedx only: which residual file")
parser$add_argument("--residfile", type = "character", default = NULL)
parser$add_argument("--corstr", type = "character", default = "unstructured")
parser$add_argument("--xtransform", type = "character", default = "z")
parser$add_argument("--stderr_variants", type = "character", default = "none", help = "none or comma list of jack,j1s,fij (slower)")
parser$add_argument("--no_perdomain", action = "store_true")
parser$add_argument("--outtag", type = "character", default = NULL)
args <- parser$parse_args()
stopifnot(args$design %in% c("fixedx", "all", "batch"), args$protect %in% c("yes", "no"), args$na %in% c("zero", "aware"))
sev <- if (args$stderr_variants %in% c("none", "")) NULL else strsplit(args$stderr_variants, ",")[[1]]
extra <- list(alltests = TRUE, stderr_variants = sev, perdomain = !args$no_perdomain)
stopifnot(args$B >= 1)
reclean <- args$design != "fixedx"

if (is.null(args$outtag))
    args$outtag <- if (!reclean) "Perm5a_fixedx" else
        paste0(ifelse(args$design == "all", "Perm5b_all", "Perm5c_batch"), "_protect", args$protect,
               ifelse(args$na == "aware", "_naAware", ""))
if (args$corstr != "unstructured") args$outtag <- paste0(args$outtag, "_", args$corstr)
if (args$xtransform != "z") args$outtag <- paste0(args$outtag, "_", args$xtransform)

outdir <- paste0("2_Pipeline/03-Run_regressions/Prot_multivar_diag_R2/", args$outtag, "/")
MyMkdir(outdir)
fileid <- paste0("prot_multivar_", args$cohort, "_", args$outtag)
repfile <- function(b) sprintf("%s%s_rep%03d.txt", outdir, fileid, b)
manifest_file <- paste0(outdir, fileid, "_manifest.txt")

dat <- PrepareData(args$cohort, args$residset, args$residfile, xtransform = args$xtransform)
transformed <- dat$transformed
phenos      <- dat$phenos

### the random subset of proteins (same set and seed as Round 1 Code 5a; the same in every replicate)
set.seed(args$seed)
nsel <- if (args$nprot > 0) min(args$nprot, nrow(transformed)) else nrow(transformed)
sel  <- sort(sample.int(nrow(transformed), nsel))
selnames <- rownames(transformed)[sel]
print(paste0("Proteins used: ", nsel, " of ", nrow(transformed)))

### data for the cleaning designs
if (reclean) {
    cin <- LoadCleaningInputs(args$cohort)
    ## align the GEE phenos (PrepareData) to the order of the cleaning phenos; the sample sets must agree
    stopifnot(setequal(cin$phenos$protsample, phenos$protsample),
              all(selnames %in% rownames(cin$edata)))
    phenos_g <- phenos[match(cin$phenos$protsample, phenos$protsample), ]
    stopifnot(all(phenos_g$protsample == cin$phenos$protsample))
    ## real-data NPS vs the NPS in the two tables must agree
    for (v in npsvars.bin) stopifnot(isTRUE(all.equal(as.numeric(phenos_g[[v]]), as.numeric(cin$phenos[[v]]))))
    edata_sel <- cin$edata[selnames, , drop = FALSE]
}

### manifest and stale-replicate guard
cur <- c(design = args$design, protect = ifelse(reclean, args$protect, NA), na = ifelse(reclean, args$na, NA),
         nprot = nsel, seed = args$seed, corstr = args$corstr, xtransform = args$xtransform,
         n_samples = nrow(phenos), residset = ifelse(reclean, NA, args$residset),
         proteins_md5 = as.character(tools::md5sum(local({ f <- tempfile(); writeLines(selnames, f); f }))),
         stderr_variants = paste(sev, collapse = ","), perdomain = extra$perdomain)
cur <- setNames(ifelse(is.na(cur) | cur == "", "-", as.character(cur)), names(cur))   # "-" = not applicable (survives fwrite/fread)
if (file.exists(manifest_file)) {
    old <- fread(manifest_file); oldv <- setNames(old$value, old$item)
    same <- all(names(cur) %in% names(oldv)) && all(mapply(function(a, b) identical(as.character(a), as.character(b)),
                                                           cur, oldv[names(cur)]))
    if (!same) {
        stale <- paste0(outdir, "stale_", format(Sys.time(), "%Y%m%d_%H%M%S"), "/")
        MyMkdir(stale)
        ff <- list.files(outdir, pattern = "_rep[0-9]+\\.txt$", full.names = TRUE)
        file.rename(ff, paste0(stale, basename(ff)))
        warning("Settings differ from the existing replicate files: ", length(ff), " files moved to ", stale)
    }
}
fwrite(data.table(item = names(cur), value = cur), manifest_file, sep = "\t", quote = FALSE)

RunOneRep <- function(b) {
    if (b > 0) set.seed(args$seed + b)
    if (!reclean) {
        ph <- copy(phenos)
        if (b > 0) ph <- PermuteNPS(ph, MakeDonor(ph, "all"))
        xm <- transformed[sel, , drop = FALSE]
    } else {
        donor <- if (b > 0) MakeDonor(cin$phenos, args$design) else seq_len(nrow(cin$phenos))
        ph <- PermuteNPS(phenos_g, donor)
        cp <- PermuteNPS(cin$phenos, donor)
        cl <- CleanOnce(edata_sel, cp, protect = args$protect, na = args$na)
        xm <- ScaleRows(cl, args$xtransform)
        keep <- intersect(selnames, rownames(xm)); xm <- xm[keep, , drop = FALSE]
    }
    rows <- vector("list", nrow(xm))
    for (j in seq_len(nrow(xm))) {
        res <- RunProteinAnalysis(x_vec = xm[j, ], sample_ids = colnames(xm), phenos = ph,
                                  npsvars = npsvars.bin, corstr = args$corstr, extra = extra)
        rows[[j]] <- ResultRow(rownames(xm)[j], res, n_na_resid = sum(is.na(xm[j, ])))
    }
    dtb <- rbindlist(lapply(rows, as.data.table), fill = TRUE)
    dtb[, rep := b]
    dtb
}

t_start <- proc.time()[["elapsed"]]
for (b in 0:args$B) {
    if (file.exists(repfile(b))) {
        chk <- tryCatch(nrow(fread(repfile(b), select = "protein")), error = function(e) -1)
        if (chk > 0 && chk <= nsel) { print(paste0("replicate ", b, " exists, skipped")); next }
        file.rename(repfile(b), paste0(repfile(b), ".bad"))
    }
    dtb <- RunOneRep(b)
    fwrite(dtb, file = repfile(b), sep = "\t", quote = FALSE)
    el <- proc.time()[["elapsed"]] - t_start
    print(sprintf("replicate %d done; %.1f sec per protein; elapsed %.1f min", b,
                  sum(dtb$elapsed_sec, na.rm = TRUE) / nrow(dtb), el / 60))
}

###### SUMMARY (the full comparison is in summarize_calibration_r2.r) ######
all <- rbindlist(lapply(0:args$B, function(b) fread(repfile(b))), fill = TRUE)
real <- all[rep == 0]; null <- all[rep > 0]
lam <- all[, .(n_p = sum(!is.na(P)), lambda_GC_HTZ = LambdaGC(P), lambda_GC_chisq_CR2 = LambdaGC(p_chisq_cr2),
               frac_p_lt_0.05_HTZ = mean(P < 0.05, na.rm = TRUE), frac_p_lt_0.001_HTZ = mean(P < 0.001, na.rm = TRUE),
               median_Q = median(Q_cr2, na.rm = TRUE), median_d2 = median(df_denom, na.rm = TRUE)), by = rep][order(rep)]
fwrite(lam, paste0(outdir, fileid, "_null_lambda_by_replicate.txt"), sep = "\t", quote = FALSE)
print(lam)
pooled <- rbind(data.table(set = "real (rep 0)", t(CohortSummary(real)[, setNames(value, metric)])),
                data.table(set = "null pooled", t(CohortSummary(null)[, setNames(value, metric)])), fill = TRUE)
fwrite(pooled, paste0(outdir, fileid, "_null_pooled_summary.txt"), sep = "\t", quote = FALSE)

png(paste0(outdir, fileid, "_null_qq.png"), height = 6, width = 12, units = "in", res = 200)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
for (v in c("P", "p_chisq_cr2")) {
    plist <- c(list(real[[v]], null[[v]]), lapply(1:args$B, function(b) all[rep == b][[v]]))
    EmptyQQ(plist, paste0(fileid, "\n", ifelse(v == "P", "CR2 + HTZ", "CR2 + chi-square"), "; ", nsel, " proteins, B = ", args$B))
    for (b in 1:args$B) AddQQ(all[rep == b][[v]], "grey70", cex = 0.25)
    AddQQ(null[[v]], "black", cex = 0.3); AddQQ(real[[v]], "red", cex = 0.5)
    legend("topleft", bty = "n", pch = 16, col = c("red", "black", "grey70"), cex = 0.8,
           legend = c("real data", "pooled permutations", "single permutations"))
}
dev.off()
print(paste0("Done. Output in ", outdir))
