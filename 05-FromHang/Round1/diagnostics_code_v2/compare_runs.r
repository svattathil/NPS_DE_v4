###############################################################################
## compare_runs.r
## Compares two runs of 03a_code1_diagnostics.r protein by protein
## (Code 1 vs Code 2, 3, 4, 6, 7).
##
## Usage:
##   Rscript 1_Code/compare_runs.r --a <stats_diag file of run A> --b <stats_diag file of run B> \
##       --labela Code1 --labelb Code2 --out <prefix>
## where the files are the *_regression_stats_diag.txt written by 03a_code1_diagnostics.r
###############################################################################
library(data.table)
library(argparse)
parser <- ArgumentParser()
parser$add_argument("--a", type = "character", required = TRUE)
parser$add_argument("--b", type = "character", required = TRUE)
parser$add_argument("--labela", type = "character", default = "A")
parser$add_argument("--labelb", type = "character", default = "B")
parser$add_argument("--out", type = "character", default = "compare")
args <- parser$parse_args()
source("1_Code/diag_functions.r")

A <- fread(args$a); B <- fread(args$b)
cols <- intersect(c("protein", "P", "df_denom", "p_chisq_cr2", "qvalue", "phi_hat"), intersect(names(A), names(B)))
M <- merge(A[, ..cols], B[, ..cols], by = "protein", suffixes = c(".a", ".b"))
lg <- function(p) -log10(pmax(p, 1e-300))

s <- list(n_common_proteins = nrow(M),
          lambda_GC_HTZ_a = LambdaGC(A$P), lambda_GC_HTZ_b = LambdaGC(B$P),
          lambda_5pct_HTZ_a = LambdaTail(A$P), lambda_5pct_HTZ_b = LambdaTail(B$P),
          frac_p_lt_0.05_a = mean(A$P < 0.05, na.rm = TRUE), frac_p_lt_0.05_b = mean(B$P < 0.05, na.rm = TRUE),
          spearman_P = cor(M$P.a, M$P.b, method = "spearman", use = "pairwise.complete.obs"),
          pearson_minus_log10P = cor(lg(M$P.a), lg(M$P.b), use = "pairwise.complete.obs"),
          mean_minus_log10P_diff_b_minus_a = mean(lg(M$P.b) - lg(M$P.a), na.rm = TRUE),
          n_abs_diff_minus_log10P_gt_0.5 = sum(abs(lg(M$P.b) - lg(M$P.a)) > 0.5, na.rm = TRUE))
if ("qvalue.a" %in% names(M)) {
    s$n_qvalue_lt_0.05_a <- sum(A$qvalue < 0.05, na.rm = TRUE)
    s$n_qvalue_lt_0.05_b <- sum(B$qvalue < 0.05, na.rm = TRUE)
    s$n_qvalue_lt_0.05_both <- sum(M$qvalue.a < 0.05 & M$qvalue.b < 0.05, na.rm = TRUE)
}
if ("df_denom.a" %in% names(M)) {
    s$median_d2_a <- median(M$df_denom.a, na.rm = TRUE); s$median_d2_b <- median(M$df_denom.b, na.rm = TRUE)
    s$spearman_d2 <- cor(M$df_denom.a, M$df_denom.b, method = "spearman", use = "pairwise.complete.obs")
}
if ("phi_hat.a" %in% names(M)) { s$median_phi_a <- median(M$phi_hat.a, na.rm = TRUE); s$median_phi_b <- median(M$phi_hat.b, na.rm = TRUE) }
out <- data.table(metric = names(s), value = unlist(s))
fwrite(out, paste0(args$out, "_compare_summary.txt"), sep = "\t", quote = FALSE)
print(out)

png(paste0(args$out, "_compare.png"), height = 5, width = 15, units = "in", res = 200)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1))
mx <- max(lg(M$P.a), lg(M$P.b), na.rm = TRUE)
plot(lg(M$P.a), lg(M$P.b), pch = 16, cex = 0.4, xlab = paste("-log10 p,", args$labela), ylab = paste("-log10 p,", args$labelb),
     xlim = c(0, mx), ylim = c(0, mx), main = "CR2 + HTZ, protein by protein"); abline(0, 1, col = "grey40")
if ("df_denom.a" %in% names(M)) {
    r <- range(M$df_denom.a, M$df_denom.b, na.rm = TRUE)
    plot(M$df_denom.a, M$df_denom.b, pch = 16, cex = 0.4, xlab = paste("d2,", args$labela), ylab = paste("d2,", args$labelb),
         xlim = r, ylim = r, main = "HTZ denominator df d2"); abline(0, 1, col = "grey40")
} else plot.new()
EmptyQQ(list(A$P, B$P), "Q-Q plot, CR2 + HTZ")
AddQQ(A$P, "red"); AddQQ(B$P, "blue")
legend("topleft", bty = "n", pch = 16, col = c("red", "blue"), legend = c(args$labela, args$labelb), cex = 0.9)
dev.off()
print(paste0("Written: ", args$out, "_compare_summary.txt and _compare.png"))
