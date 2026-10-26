###############################################################################
## tests/test_make_synthetic_data_r2.R
## Makes a SYNTHETIC "emory" cohort in a scratch project folder so that the Round 2 scripts can be
## tested without real data. Run from the scratch project directory (NOT your real project):
##   Rscript tests/test_make_synthetic_data_r2.R
## Writes: 2_Pipeline/01-Explore_data/{prot.pcafiltered.log2norm.txt, phenos_cleaned_emory.txt}
##         2_Pipeline/02-Prepare_analysis_data/{phenos_cleaned_SVs_emory.txt,
##             resid_regress_BatchPmiAgeMsex_emory.txt, resid_regress_BatchPmiAgeMsex_noProtect_emory.txt}
## The two residual files are made with the cleaning replica (cleaning_functions_r2.r), so with
## default options 03a_run_r2.r must reproduce 03a's P exactly.
###############################################################################
library(data.table)
source("1_Code/project_constants.r")
source("1_Code/cleaning_functions_r2.r")
set.seed(11)
n <- 73; G <- 300
ids <- sprintf("S%03d", 1:n)
batch <- sample(1:4, n, TRUE)
z <- rnorm(n)
ph <- data.table(protsample = ids, Batch = batch, pmi = rnorm(n, 8, 3), age_death = rnorm(n, 85, 6),
                 msex = rbinom(n, 1, 0.5), sv1 = rnorm(n))
for (k in seq_along(npsvars.bin))
    set(ph, j = npsvars.bin[k], value = rbinom(n, 1, plogis(-0.6 + 0.2 * k + 0.5 * z + 0.4 * (batch == 2))))
ph[sample(n, 3), (npsvars.bin[1]) := NA]   # 3 subjects lack one NPS value -> excluded, as in script V
X <- matrix(rnorm(G * n), G, n) + 0.5 * outer(rnorm(G), ph$pmi - mean(ph$pmi)) / 3 +
     outer(rnorm(G, 0, 0.5), (batch == 3) * 1)
X[sample(length(X), 0.04 * length(X))] <- NA
X[7, sample(n, 30)] <- NA                      # a protein with too many missing (< 50 observed): dropped
X[1:5, ] <- X[1:5, ] + 0.9 * matrix(rep(z, 5), 5, byrow = TRUE)
X[10, ] <- rt(n, 3)                            # a heavy-tailed protein
colnames(X) <- ids
expr <- data.table(protein = sprintf("P%04d", 1:G), X)
dir.create("2_Pipeline/01-Explore_data", recursive = TRUE, showWarnings = FALSE)
dir.create("2_Pipeline/02-Prepare_analysis_data", recursive = TRUE, showWarnings = FALSE)
phraw <- copy(ph); phraw[, Batch := as.integer(Batch)]
fwrite(phraw, "2_Pipeline/01-Explore_data/phenos_cleaned_emory.txt", sep = "\t")
fwrite(expr, "2_Pipeline/01-Explore_data/prot.pcafiltered.log2norm.txt", sep = "\t")
fwrite(ph[, c("protsample", npsvars.bin), with = FALSE], "2_Pipeline/02-Prepare_analysis_data/phenos_cleaned_SVs_emory.txt", sep = "\t")

## "existing" residual files, made with the cleaning replica
cin <- LoadCleaningInputs("emory")
for (pr in c("yes", "no")) {
    cl <- CleanOnce(cin$edata, cin$phenos, protect = pr, na = "zero")
    out <- data.table(protein = rownames(cl), cl)
    fwrite(out, Residfile(ifelse(pr == "yes", "ignoreSVs", "noProtect"), "emory"), sep = "\t")
}
cat("synthetic emory data written; n samples used:", nrow(cin$phenos), "proteins:", nrow(cin$edata), "\n")
