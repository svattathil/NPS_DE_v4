###############################################################################
## summarize_perm_r3.r                                                (Round 3)
##
## Reads the replicate files of ONE permutation folder made by 03a_perm_null_r3.r
## (rep000 = real data, rep001... = permutations) and writes
##   <prefix>_r3_table.txt        one row per reference: lambda_GC and the share of p-values below 0.05 / 0.01 / 0.001
##                                (null pooled, with the range over the single permutations; real data)
##   <prefix>_r3_by_d2.txt        the same shares by tertile of the denominator df d2 of HTZ (null pooled)
##   <prefix>_r3_by_replicate.txt lambda_GC and share below 0.05 for every replicate and reference
##   <prefix>_r3_pvalues.txt      the p-values of every protein and replicate for every reference (to redraw the figures)
##   <prefix>_r3_qq_curves.txt    numbers of the Q-Q curves (about 200 points per curve): real, pooled null, and the
##                                lowest / highest curve among the single permutations
##   <prefix>_r3_qq_2x3.png       the Q-Q plots, 2 x 3 panels (red = real, black = pooled permutations, grey = single ones)
##
## The references are the SAME CR2 statistic with different reference distributions.
##
## Usage (from the project directory):
##   Rscript 1_Code/summarize_perm_r3.r --dir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes
##   optional:  --full <path>/prot_multivar_..._regression_stats_diag.txt   (the full run; then the red dots use all proteins)
###############################################################################
library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--dir", type = "character", required = TRUE, help = "folder with the *_repNNN.txt files")
parser$add_argument("--full", type = "character", default = NULL,
                    help = "optional: regression_stats_diag.txt of the full run (all proteins) for the red dots")
parser$add_argument("--prefix", type = "character", default = NULL, help = "output prefix (default: <dir>/<folder name>)")
args <- parser$parse_args()

dir <- sub("/*$", "/", args$dir)
files <- list.files(dir, pattern = "_rep[0-9]+\\.txt$", full.names = TRUE)
stopifnot(length(files) >= 2)
repno <- as.integer(sub(".*_rep([0-9]+)\\.txt$", "\\1", files))
files <- files[order(repno)]; repno <- sort(repno)
stopifnot(repno[1] == 0)
if (any(diff(repno) != 1)) warning("Replicate numbers are not consecutive: ", paste(setdiff(0:max(repno), repno), collapse = ","), " missing")
prefix <- if (is.null(args$prefix)) paste0(dir, basename(sub("/$", "", dir))) else args$prefix
all <- rbindlist(lapply(files, fread), fill = TRUE)
real <- all[rep == 0]; null <- all[rep > 0]
B <- length(unique(null$rep))
cat(sprintf("%d replicate files; real data: %d proteins; permutations: B = %d, %d p-values\n", length(files), nrow(real), B, nrow(null)))
realfull <- if (!is.null(args$full)) fread(args$full) else NULL

refs <- c("HTZ" = "P", "HTA" = "p_HTA", "HTB" = "p_HTB", "EDF" = "p_EDF", "naive F" = "p_NaiveF",
          "chi2(7), CR2" = "p_chisq_cr2", "chi2(7), san.se" = "p_chisq_san")
refs <- refs[refs %in% names(all)]

lam <- function(p) { p <- p[!is.na(p)]; if (!length(p)) return(NA_real_); p[p == 0] <- 1e-300
                     median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1) }
shr <- function(p, a) { p <- p[!is.na(p)]; if (!length(p)) NA_real_ else mean(p < a) }

###### table 1: by reference ######
tab <- rbindlist(lapply(names(refs), function(nm) {
    v <- refs[[nm]]
    per <- null[, .(l = lam(get(v)), s05 = shr(get(v), 0.05), s001 = shr(get(v), 0.001)), by = rep]
    r <- data.table(reference = nm, n_null = sum(!is.na(null[[v]])), n_real = sum(!is.na(real[[v]])),
        lambda_null = lam(null[[v]]), lambda_null_rep_min = min(per$l, na.rm = TRUE), lambda_null_rep_max = max(per$l, na.rm = TRUE),
        share05_null = shr(null[[v]], 0.05), share05_null_rep_min = min(per$s05, na.rm = TRUE), share05_null_rep_max = max(per$s05, na.rm = TRUE),
        share01_null = shr(null[[v]], 0.01),
        share001_null = shr(null[[v]], 0.001), share001_null_rep_min = min(per$s001, na.rm = TRUE), share001_null_rep_max = max(per$s001, na.rm = TRUE),
        lambda_real = lam(real[[v]]), share05_real = shr(real[[v]], 0.05), share001_real = shr(real[[v]], 0.001))
    if (!is.null(realfull) && v %in% names(realfull))
        r[, `:=`(lambda_real_full = lam(realfull[[v]]), share05_real_full = shr(realfull[[v]], 0.05), n_real_full = sum(!is.na(realfull[[v]])))]
    r
}), fill = TRUE)
fwrite(tab, paste0(prefix, "_r3_table.txt"), sep = "\t", quote = FALSE)
print(tab[, .(reference, lambda_null = round(lambda_null, 2), share05_null = round(share05_null, 4), share001_null = round(share001_null, 4),
              lambda_real = round(lambda_real, 2))])

###### table 2: by tertile of d2 (HTZ denominator df, from the null) ######
nl <- copy(null)[!is.na(df_denom)]
qs <- quantile(nl$df_denom, 0:3 / 3)
nl[, tert := cut(df_denom, unique(qs), include.lowest = TRUE, labels = c("low d2", "mid d2", "high d2")[seq_len(length(unique(qs)) - 1)])]
bytert <- rbindlist(lapply(names(refs), function(nm) {
    v <- refs[[nm]]
    nl[, .(reference = nm, n = sum(!is.na(get(v))), median_d2 = median(df_denom), share05 = shr(get(v), 0.05), share001 = shr(get(v), 0.001),
           lambda = lam(get(v))), by = tert]
}))
fwrite(bytert, paste0(prefix, "_r3_by_d2.txt"), sep = "\t", quote = FALSE)

###### the p-values themselves (small file; no individual-level data) ######
pv <- all[, c("rep", "protein", "df_denom", "Q_cr2", unname(refs)), with = FALSE]
fwrite(pv, paste0(prefix, "_r3_pvalues.txt"), sep = "\t", quote = FALSE)

###### table 3: by replicate ######
byrep <- rbindlist(lapply(names(refs), function(nm) {
    v <- refs[[nm]]; all[, .(reference = nm, n = sum(!is.na(get(v))), lambda = lam(get(v)), share05 = shr(get(v), 0.05),
                              share001 = shr(get(v), 0.001), median_Q_cr2 = median(Q_cr2, na.rm = TRUE)), by = rep]
}))
fwrite(byrep[order(rep, reference)], paste0(prefix, "_r3_by_replicate.txt"), sep = "\t", quote = FALSE)

###### Q-Q curves on a common grid: e = expected -log10(p); observed -log10 of the p-value at upper-tail fraction 10^-e ######
obs_at <- function(p, e) {
    p <- sort(p[!is.na(p)]); p[p == 0] <- 1e-300; m <- length(p)
    r <- ceiling(10^-e * m)
    ok <- r >= 1 & r <= m
    out <- rep(NA_real_, length(e)); out[ok] <- -log10(p[r[ok]]); out
}
egrid <- seq(0, 4.4, by = 0.02)
curves <- rbindlist(lapply(names(refs), function(nm) {
    v <- refs[[nm]]
    reps <- sapply(split(null[[v]], null$rep), obs_at, e = egrid)           # grid x B
    d <- data.table(reference = nm, expected = egrid, real_500 = obs_at(real[[v]], egrid), null_pooled = obs_at(null[[v]], egrid),
                    null_rep_min = apply(reps, 1, function(z) if (all(is.na(z))) NA_real_ else min(z, na.rm = TRUE)),
                    null_rep_max = apply(reps, 1, function(z) if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE)))
    if (!is.null(realfull) && v %in% names(realfull)) d[, real_full := obs_at(realfull[[v]], egrid)]
    d
}), fill = TRUE)
fwrite(curves, paste0(prefix, "_r3_qq_curves.txt"), sep = "\t", quote = FALSE)

###### figure: 2 x 3, red = real, black = pooled permutations, grey = single permutations ######
panel <- c("HTZ", "HTA", "HTB", "EDF", "naive F", "chi2(7), CR2")
lab   <- c("HTZ", "HTA", "HTB", "EDF", "naive F", expression(chi[7]^2))
qq <- function(p) { p <- sort(p[!is.na(p)]); p[p == 0] <- 1e-300; m <- length(p)
                    list(e = -log10((seq_len(m) - 0.5) / (m + 1)), o = -log10(p), m = m) }
png(paste0(prefix, "_r3_qq_2x3.png"), height = 8.4, width = 12, units = "in", res = 200)
par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.5, 1), pty = "s")
for (k in seq_along(panel)) {
    if (!panel[k] %in% names(refs)) { plot.new(); next }
    v <- refs[[panel[k]]]
    rl <- if (!is.null(realfull) && v %in% names(realfull)) realfull[[v]] else real[[v]]
    cr <- qq(rl); cn <- qq(null[[v]])
    mx <- max(cr$e, cr$o, cn$e, cn$o, 1)
    plot(NA, xlim = c(0, mx), ylim = c(0, mx), xlab = expression(Expected~-log[10](p)), ylab = expression(Observed~-log[10](p)),
         main = lab[[k]], cex.main = 1.1)
    abline(0, 1, col = "grey40")
    for (b in unique(null$rep)) { cb <- qq(null[rep == b][[v]]); points(cb$e, cb$o, col = "grey75", pch = 16, cex = 0.25) }
    points(cn$e, cn$o, col = "black", pch = 16, cex = 0.35)
    points(cr$e, cr$o, col = "red", pch = 16, cex = 0.45)
    if (k == 1) legend("topleft", bty = "n", pch = 16, col = c("red", "black", "grey75"), cex = 0.9,
                       legend = c(ifelse(is.null(realfull), "real data (500 proteins)", "real data"), "pooled permutations", "single permutations"))
}
dev.off()
cat("Written: ", prefix, "_r3_table.txt, _r3_by_d2.txt, _r3_by_replicate.txt, _r3_pvalues.txt, _r3_qq_curves.txt, _r3_qq_2x3.png\n", sep = "")
