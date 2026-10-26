###############################################################################
## summarize_calibration_r2.r                                         (Round 2)
##
## Summarizes calibration from the Round 2 outputs, for any number of runs.
##   --perm label=<folder>   permutation-null folder written by 03a_perm_null_r2.r (rep000 = real data
##                           for the same proteins, rep001... = permutations); repeatable
##   --full  label=<file>    *_regression_stats_diag.txt of a full run of 03a_run_r2.r; repeatable
##                           (used for the Q-Q curves and the BH counts of that run)
##   --out   <prefix>        prefix of the output files (a folder is created if needed)
## Example:
##   Rscript 1_Code/summarize_calibration_r2.r --out 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R2/summary \
##     --perm 5a=.../Perm5a_fixedx --perm 5b_protect=.../Perm5b_all_protectyes --perm 5b_noprot=.../Perm5b_all_protectno \
##     --full base=.../R2_IgnoreSVs/prot_multivar_ignoreSVs_emory_R2_regression_stats_diag.txt
## Tables written (all tab-separated):
##   *_lambda_table.txt    lambda_GC and fractions of p<0.05/0.01/0.001 for every p-value column, real and null
##   *_Q_ratio.txt         quantiles of Q (CR2, 7 df): real and null, with their ratios to the chi2_7 quantile
##                         and real/null ratio -- the "bulk shift" of Round 1
##   *_by_d2_tertile.txt   null rejection rates and median Q by tertile of the HTZ denominator df d2
##   *_perdomain.txt       lambda_GC and p<0.05 fraction of the seven single-domain tests, real and null
##   *_scaled_chisq.txt    recalibration p = pchisq(Q/c, 7), c = null median Q / chi2_7 median; BH counts
##                         of the real data at q<0.05/0.10/0.20 (post-hoc: NOT a valid analysis)
##   *_qq.png              Q-Q curves (real red, null black) of the main p-values, one panel per set
###############################################################################
library(data.table)
library(argparse)
source("1_Code/diag_functions_r2.r", chdir = FALSE)

parser <- ArgumentParser()
parser$add_argument("--perm", action = "append", default = character(0))
parser$add_argument("--full", action = "append", default = character(0))
parser$add_argument("--out", type = "character", default = "summary")
args <- parser$parse_args()
dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)
spl <- function(x) { i <- regexpr("=", x, fixed = TRUE); c(substr(x, 1, i - 1), substr(x, i + 1, nchar(x))) }

## --- load sets: list of list(label, real, null) ---
sets <- list()
for (a in args$perm) {
    kv <- spl(a); ff <- list.files(kv[2], pattern = "_rep[0-9]+\\.txt$", full.names = TRUE)
    stopifnot(length(ff) >= 2)
    all <- rbindlist(lapply(ff, fread), fill = TRUE)
    sets[[kv[1]]] <- list(real = all[rep == 0], null = all[rep > 0], B = uniqueN(all[rep > 0]$rep))
}
fulls <- lapply(args$full, function(a) { kv <- spl(a); list(label = kv[1], dt = fread(kv[2])) })
stopifnot(length(sets) + length(fulls) > 0)

pcols <- function(dt) intersect(c("P", "p_chisq", "p_chisq_cr2", "p_chisq_cr0", "p_chisq_san", "p_NaiveF", "p_NaiveFp",
                                  "p_HTA", "p_HTB", "p_HTZ", "p_EDF", "p_EDT", "p_chisq_jack", "p_chisq_j1s", "p_chisq_fij"), names(dt))
lamrow <- function(dt, label, which) {
    rbindlist(lapply(pcols(dt), function(v) data.table(set = label, data = which, pvalue = v, n = sum(!is.na(dt[[v]])),
        lambda_GC = LambdaGC(dt[[v]]), lambda_5pct = LambdaTail(dt[[v]]),
        frac_lt_0.05 = mean(dt[[v]] < 0.05, na.rm = TRUE), frac_lt_0.01 = mean(dt[[v]] < 0.01, na.rm = TRUE),
        frac_lt_0.001 = mean(dt[[v]] < 0.001, na.rm = TRUE))))
}
## 1. lambda table
lt <- rbindlist(c(lapply(names(sets), function(s) rbind(lamrow(sets[[s]]$real, s, "real"), lamrow(sets[[s]]$null, s, "null"))),
                  lapply(fulls, function(f) lamrow(f$dt, f$label, "full real"))))
fwrite(lt, paste0(args$out, "_lambda_table.txt"), sep = "\t", quote = FALSE)

## 2. Q ratio table (CR2 Q, 7 df)
qs <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
qr_rows <- list()
for (s in names(sets)) {
    r <- sets[[s]]$real$Q_cr2; n <- sets[[s]]$null$Q_cr2
    qr_rows[[s]] <- data.table(set = s, quantile = qs, chisq7 = qchisq(qs, 7),
        Q_real = quantile(r, qs, na.rm = TRUE), Q_null = quantile(n, qs, na.rm = TRUE))[,
        `:=`(null_over_chisq7 = Q_null / chisq7, real_over_chisq7 = Q_real / chisq7, real_over_null = Q_real / Q_null)]
}
for (f in fulls) if ("Q_cr2" %in% names(f$dt))
    qr_rows[[f$label]] <- data.table(set = f$label, quantile = qs, chisq7 = qchisq(qs, 7),
        Q_real = quantile(f$dt$Q_cr2, qs, na.rm = TRUE), Q_null = NA_real_)[, `:=`(null_over_chisq7 = NA_real_, real_over_chisq7 = Q_real / chisq7, real_over_null = NA_real_)]
qr <- rbindlist(qr_rows); fwrite(qr, paste0(args$out, "_Q_ratio.txt"), sep = "\t", quote = FALSE)

## 3. by tertile of d2 (null)
tert <- rbindlist(lapply(names(sets), function(s) {
    nl <- copy(sets[[s]]$null)[!is.na(df_denom)]
    nl[, tert := cut(df_denom, quantile(df_denom, 0:3 / 3, na.rm = TRUE), include.lowest = TRUE, labels = c("low d2", "mid d2", "high d2"))]
    nl[, .(set = s, n = .N, median_d2 = median(df_denom), median_Q = median(Q_cr2, na.rm = TRUE),
           frac_HTZ_lt_0.05 = mean(P < 0.05), frac_chisq_lt_0.05 = mean(p_chisq_cr2 < 0.05, na.rm = TRUE),
           median_kurt_x = if ("kurt_x" %in% names(nl)) median(kurt_x, na.rm = TRUE) else NA_real_), by = tert]
}))
fwrite(tert, paste0(args$out, "_by_d2_tertile.txt"), sep = "\t", quote = FALSE)

## 4. per domain
pdc <- function(dt) grep("^pdom_", names(dt), value = TRUE)
pd <- rbindlist(lapply(names(sets), function(s) rbindlist(lapply(c("real", "null"), function(w) {
    dt <- sets[[s]][[w]]; v <- pdc(dt); if (!length(v)) return(NULL)
    data.table(set = s, data = w, domain = sub("^pdom_", "", v), lambda_GC = sapply(v, function(x) LambdaGC(dt[[x]])),
               frac_lt_0.05 = sapply(v, function(x) mean(dt[[x]] < 0.05, na.rm = TRUE)))
}))), fill = TRUE)
if (nrow(pd)) fwrite(pd, paste0(args$out, "_perdomain.txt"), sep = "\t", quote = FALSE)

## 5. scaled chi-square recalibration (post-hoc; c from the null median of Q)
sc <- rbindlist(lapply(names(sets), function(s) {
    cc <- median(sets[[s]]$null$Q_cr2, na.rm = TRUE) / qchisq(0.5, 7)
    r <- sets[[s]]$real; p <- pchisq(r$Q_cr2 / cc, 7, lower.tail = FALSE); q <- p.adjust(p, "BH")
    nn <- sets[[s]]$null; pn <- pchisq(nn$Q_cr2 / cc, 7, lower.tail = FALSE)
    data.table(set = s, source = "real rep0 and null of this folder", c_scale = cc, n = sum(!is.na(p)),
               lambda_GC_real_scaled = LambdaGC(p), lambda_GC_null_scaled = LambdaGC(pn),
               frac_null_lt_0.05 = mean(pn < 0.05, na.rm = TRUE), n_BH_q05 = sum(q < 0.05, na.rm = TRUE),
               n_BH_q10 = sum(q < 0.10, na.rm = TRUE), n_BH_q20 = sum(q < 0.20, na.rm = TRUE))
}))
for (f in fulls) if ("Q_cr2" %in% names(f$dt) && nrow(sc)) for (s in names(sets)) {
    cc <- sc[set == s, c_scale]; p <- pchisq(f$dt$Q_cr2 / cc, 7, lower.tail = FALSE); q <- p.adjust(p, "BH")
    sc <- rbind(sc, data.table(set = s, source = paste0("full run ", f$label), c_scale = cc, n = sum(!is.na(p)),
        lambda_GC_real_scaled = LambdaGC(p), lambda_GC_null_scaled = NA_real_, frac_null_lt_0.05 = NA_real_,
        n_BH_q05 = sum(q < 0.05, na.rm = TRUE), n_BH_q10 = sum(q < 0.10, na.rm = TRUE), n_BH_q20 = sum(q < 0.20, na.rm = TRUE)))
}
if (nrow(sc)) fwrite(sc, paste0(args$out, "_scaled_chisq.txt"), sep = "\t", quote = FALSE)

## 6. Q-Q plot
panels <- c(names(sets), sapply(fulls, `[[`, "label"))
png(paste0(args$out, "_qq.png"), height = 5.2, width = 5.2 * min(length(panels), 3), units = "in", res = 150)
par(mfrow = c(ceiling(length(panels) / 3), min(length(panels), 3)), mar = c(4.2, 4.2, 3, 1))
for (s in panels) {
    if (s %in% names(sets)) { r <- sets[[s]]$real$P; n <- sets[[s]]$null$P } else { r <- fulls[[which(sapply(fulls, `[[`, "label") == s)]]$dt$P; n <- NULL }
    EmptyQQ(list(r, n), paste0(s, "  (CR2 + HTZ)"))
    if (!is.null(n)) AddQQ(n, "black", cex = 0.3)
    AddQQ(r, "red", cex = 0.5)
    legend("topleft", bty = "n", pch = 16, col = c("red", "black"), cex = 0.8, legend = c("real", "permutation null"))
}
dev.off()
print(lt[pvalue %in% c("P", "p_chisq_cr2", "p_chisq_jack")]); print(qr[quantile %in% c(0.5, 0.95)]); print(tert); if (nrow(sc)) print(sc)
cat("Written with prefix", args$out, "\n")
