###############################################################################
## phenotype_check_r4.r                                                (Round 4)
##
## Phenotype-level table (no proteins are used in the results): how are the seven NPS domains related to the
## nuisance variables that the cleaning removes from the expression (Batch, PMI, age at death, sex), and to each other?
## Useful to judge whether a domain-specific excess of small p-values could be carried by a hidden factor
## that follows Batch / PMI / age / sex. Same subjects as the permutation run (LoadCleaningInputs).
##
## Usage (from the project directory):
##   Rscript 1_Code/phenotype_check_r4.r --cohort rush --outdir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes_pdom_rush
## Writes  <outdir>/<cohort>_r4_pheno_domains.txt, _cases_by_batch.txt, _phi.txt
###############################################################################
library(data.table)
library(argparse)
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")
source("1_Code/cleaning_functions_r2.r")

parser <- ArgumentParser()
parser$add_argument("--cohort", type = "character", required = TRUE)
parser$add_argument("--outdir", type = "character", required = TRUE)
args <- parser$parse_args()
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

cin <- LoadCleaningInputs(args$cohort)
ph <- cin$phenos
cat(sprintf("%d subjects; Batch sizes: %s\n", nrow(ph), paste(names(table(ph$Batch)), table(ph$Batch), sep = "=", collapse = ", ")))

wald_p <- function(y, x) tryCatch({ f <- suppressWarnings(glm(y ~ x, family = binomial)); summary(f)$coefficients[2, 4] }, error = function(e) NA_real_)
wald_b <- function(y, x) tryCatch({ f <- suppressWarnings(glm(y ~ x, family = binomial)); unname(coef(f)[2]) }, error = function(e) NA_real_)
set.seed(1)
res <- rbindlist(lapply(npsvars.bin, function(v) {
    y <- as.integer(ph[[v]])
    tb <- table(ph$Batch, y)
    pb <- tryCatch(fisher.test(tb, simulate.p.value = TRUE, B = 1e5)$p.value, error = function(e) NA_real_)
    ps <- tryCatch(fisher.test(table(ph$msex, y))$p.value, error = function(e) NA_real_)
    data.table(domain = v, n = length(y), cases = sum(y), p_Batch = pb,
               p_pmi = wald_p(y, ph$pmi), log_or_per_unit_pmi = wald_b(y, ph$pmi),
               p_age_death = wald_p(y, ph$age_death), log_or_per_unit_age = wald_b(y, ph$age_death),
               p_sex = ps)
}))
fwrite(res, file.path(args$outdir, paste0(args$cohort, "_r4_pheno_domains.txt")), sep = "\t", quote = FALSE)
cb <- rbindlist(lapply(npsvars.bin, function(v) ph[, .(domain = v, n = .N, cases = sum(as.integer(get(v)))), by = Batch]))
fwrite(cb, file.path(args$outdir, paste0(args$cohort, "_r4_pheno_cases_by_batch.txt")), sep = "\t", quote = FALSE)
phi <- cor(sapply(npsvars.bin, function(v) as.integer(ph[[v]])))
pl <- CJ(d1 = seq_along(npsvars.bin), d2 = seq_along(npsvars.bin))[d1 < d2]
pl[, `:=`(domain1 = npsvars.bin[d1], domain2 = npsvars.bin[d2], phi = phi[cbind(d1, d2)])]
fwrite(pl[, !c("d1", "d2")], file.path(args$outdir, paste0(args$cohort, "_r4_pheno_phi.txt")), sep = "\t", quote = FALSE)
print(res)
cat("Written to ", args$outdir, "\n")
