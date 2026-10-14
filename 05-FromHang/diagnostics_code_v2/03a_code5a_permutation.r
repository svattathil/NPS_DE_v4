###############################################################################
## 03a_code5a_permutation.r        (Code 5a: permutation null, test only)
##
## For replicate b = 1..B the NPS vectors (all 7 domains together) are permuted
## among subjects, the cleaned expression x_ig is NOT recomputed, and the test
## of 03a is applied to a random subset of proteins. Under the permutation every
## null hypothesis is true, so the p-values show the calibration of the test.
## Replicate 0 is the real (unpermuted) data for the same proteins, run through
## the same code.
##
## Usage (from the project directory):
##   Rscript 1_Code/03a_code5a_permutation.r --cohort Emory --B 20 --nprot 1000
## Smoke test first:
##   Rscript 1_Code/03a_code5a_permutation.r --cohort Emory --B 2 --nprot 30
##
## Each replicate is written to its own file; if the job is stopped, running the
## same command again continues with the replicates that are missing.
###############################################################################
library(data.table)
library(argparse)
library(car)
library(geepack)
library(clubSandwich)

rm(list=ls())
options(stringsAsFactors = FALSE)
options(nwarnings = 100000000)
source("1_Code/project_constants.r")
source("~/comsv/Svattathil_Library/svattathil_functions.r")
source("~/comsv/Svattathil_Library/manhattan_and_qq_plot_functions.r")
source("1_Code/diag_functions.r")

parser <- ArgumentParser()
parser$add_argument("--cohort", type="character", default = "OHSU")
parser$add_argument("--residset", type="character", default = "ignoreSVs")
parser$add_argument("--residfile", type="character", default = NULL,
                    help="path of a residual file; overrides Residfile(residset, cohort)")
parser$add_argument("--B", type="integer", default = 20, help="number of permutations")
parser$add_argument("--nprot", type="integer", default = 1000,
                    help="random subset of proteins (0 = all proteins)")
parser$add_argument("--seed", type="integer", default = 20261011)
parser$add_argument("--outtag", type="character", default = "Code5a")
args <- parser$parse_args()

outdir <- paste0("2_Pipeline/03-Run_regressions/Prot_multivar_diag/",
                 args$outtag, "_", Capwords(args$residset), "/")
MyMkdir(outdir)
fileid <- paste0("prot_multivar_", args$residset, "_", args$cohort, "_", args$outtag)
repfile <- function(b) sprintf("%s%s_rep%03d.txt", outdir, fileid, b)

dat <- PrepareData(args$cohort, args$residset, args$residfile)
transformed <- dat$transformed
phenos      <- dat$phenos

### the random subset of proteins (the same in every replicate)
set.seed(args$seed)
nsel <- if (args$nprot > 0) min(args$nprot, nrow(transformed)) else nrow(transformed)
sel  <- sort(sample.int(nrow(transformed), nsel))
print(paste0("Proteins used: ", nsel, " of ", nrow(transformed)))

t_start <- proc.time()[["elapsed"]]
for (b in 0:args$B) {
    if (file.exists(repfile(b))) { print(paste0("replicate ", b, " exists, skipped")); next }

    ## b = 0: real data; b >= 1: permute the NPS vectors among subjects
    phenos_b <- copy(phenos)
    if (b > 0) {
        set.seed(args$seed + b)
        perm <- sample.int(nrow(phenos))
        for (v in npsvars.bin) set(phenos_b, j = v, value = phenos[[v]][perm])
    }

    rows <- vector("list", length(sel))
    for (j in seq_along(sel)) {
        i <- sel[j]
        res <- RunProteinAnalysis(x_vec = transformed[i, ], sample_ids = colnames(transformed),
                                  phenos = phenos_b, npsvars = npsvars.bin)
        rows[[j]] <- ResultRow(rownames(transformed)[i], res,
                               n_na_resid = sum(is.na(transformed[i, ])))
    }
    dtb <- rbindlist(lapply(rows, as.data.table), fill = TRUE)
    dtb[, rep := b]
    fwrite(dtb, file = repfile(b), sep = "\t", quote = FALSE)

    el <- proc.time()[["elapsed"]] - t_start
    print(sprintf("replicate %d done; %.1f sec per protein fit in this replicate; elapsed %.1f min",
                  b, sum(dtb$elapsed_sec, na.rm = TRUE) / nrow(dtb), el / 60))
}


###### SUMMARY ######
all <- rbindlist(lapply(0:args$B, function(b) fread(repfile(b))), fill = TRUE)
real <- all[rep == 0]
null <- all[rep > 0]

lam <- all[, .(n_p = sum(!is.na(P)),
               lambda_GC_HTZ        = LambdaGC(P),
               lambda_GC_chisq_CR2  = LambdaGC(p_chisq_cr2),
               lambda_GC_chisq_san  = LambdaGC(p_chisq_san),
               frac_p_lt_0.05_HTZ   = mean(P < 0.05, na.rm = TRUE),
               frac_p_lt_0.01_HTZ   = mean(P < 0.01, na.rm = TRUE),
               frac_p_lt_0.001_HTZ  = mean(P < 0.001, na.rm = TRUE),
               median_d2            = median(df_denom, na.rm = TRUE)), by = rep][order(rep)]
fwrite(lam, file = paste0(outdir, fileid, "_null_lambda_by_replicate.txt"), sep = "\t", quote = FALSE)
print(lam)

pooled <- rbind(data.table(set = "real (rep 0)", t(CohortSummary(real)[, setNames(value, metric)])),
                data.table(set = "null pooled",  t(CohortSummary(null)[, setNames(value, metric)])),
                fill = TRUE)
fwrite(pooled, file = paste0(outdir, fileid, "_null_pooled_summary.txt"), sep = "\t", quote = FALSE)

### empirical p-values of the real data from the pooled null (caveat: proteins are correlated)
nullp <- sort(null$P[!is.na(null$P)])
real[, P_emp := (1 + findInterval(P, nullp)) / (1 + length(nullp))]
fwrite(real, file = paste0(outdir, fileid, "_real_with_empirical_p.txt"), sep = "\t", quote = FALSE)

### Q-Q plot: null replicates (grey), pooled null (black), real data (red)
png(paste0(outdir, fileid, "_null_qq.png"), height = 6, width = 12, units = "in", res = 200)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
for (v in c("P", "p_chisq_cr2")) {
    plist <- c(list(real[[v]], null[[v]]), lapply(1:args$B, function(b) all[rep == b][[v]]))
    EmptyQQ(plist, paste0(fileid, "\n", ifelse(v == "P", "CR2 + HTZ", "CR2 + chi-square"),
                          "; ", nsel, " proteins, B = ", args$B))
    for (b in 1:args$B) AddQQ(all[rep == b][[v]], "grey70", cex = 0.25)
    AddQQ(null[[v]], "black", cex = 0.3)
    AddQQ(real[[v]], "red", cex = 0.5)
    legend("topleft", bty = "n", pch = 16, col = c("red", "black", "grey70"), cex = 0.8,
           legend = c("real data", "pooled permutations", "single permutations"))
}
dev.off()
print(paste0("Done. Output in ", outdir))
