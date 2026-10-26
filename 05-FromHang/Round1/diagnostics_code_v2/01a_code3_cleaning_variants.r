###############################################################################
## 01a_code3_cleaning_variants.r                                        (Code 3)
##
## Variants of the covariate-regression (cleaning) step of
## 01a-estimate_SVs_make_adjusted_residuals.r (script V), for the conjectures
##   C6: the NPS outcomes are protected and the cleaning's degrees of freedom are ignored
##   C7: missing expression values are treated as 0 by jaffelab::cleaningY
## Sample and protein filters, nuisance variables and model matrix are those of
## script V (option ignoreSVs / noProtect). The surrogate variables are NOT needed
## here, so sva is not run.
##
## Two switches, giving the 2 x 2 design of the plan:
##   --protect  yes | no     protect intercept + NPS outcomes (ignoreSVs) or only the intercept (noProtect)
##   --na       zero | aware zero  = missing -> 0 when estimating coefficients (what jaffelab does)
##                           aware = coefficients estimated from the observed values of each protein
##
## Output (2_Pipeline/02-Prepare_analysis_data/), read by 03a_code1_diagnostics.r through --residfile:
##   resid_regress_BatchPmiAgeMsex_naAware_<cohort>.txt              (protect yes, na aware)
##   resid_regress_BatchPmiAgeMsex_noProtect_naAware_<cohort>.txt    (protect no,  na aware)
##   resid_regress_BatchPmiAgeMsex_Code3check_<cohort>.txt           (protect yes, na zero; replica)
##   resid_regress_BatchPmiAgeMsex_noProtect_Code3check_<cohort>.txt (protect no,  na zero; replica)
## The 'zero' files are our own replica of jaffelab::cleaningY. They are compared with the
## existing files made by script V; they should agree to ~1e-8, which also shows whether the
## jaffelab version used replaces missing values by 0. Existing files are never overwritten.
## Bookkeeping (Code 3b) is written to *_code3_bookkeeping.txt, *_batch_by_nps.txt,
## *_na_per_protein.txt (aggregate numbers only).
##
## Usage:  Rscript 1_Code/01a_code3_cleaning_variants.r --cohort emory --protect yes --na aware
##         (cohort names as in script V: OHSU, rush, emory)
###############################################################################
library(data.table)
library(argparse)

rm(list = ls())
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")
source("~/comsv/Svattathil_Library/svattathil_functions.r")

nuisance_vars <- c("Batch", "pmi", "age_death", "msex")
npsvars <- c(npsvars.bin)

parser <- ArgumentParser()
parser$add_argument("--cohort", type = "character", default = "OHSU", help = "OHSU, rush, or emory")
parser$add_argument("--protect", type = "character", default = "yes", help = "yes or no")
parser$add_argument("--na", type = "character", default = "aware", help = "zero or aware")
parser$add_argument("--bookkeeping_only", action = "store_true", help = "write only the Code 3b tables")
args <- parser$parse_args()
stopifnot(args$protect %in% c("yes", "no"), args$na %in% c("zero", "aware"))

label <- paste0("BatchPmiAgeMsex", ifelse(args$protect == "no", "_noProtect", ""),
                ifelse(args$na == "aware", "_naAware", "_Code3check"))
indir  <- "2_Pipeline/01-Explore_data/"
outdir <- "2_Pipeline/02-Prepare_analysis_data/"
infiles <- list(prot   = paste0(indir, "prot.pcafiltered.log2norm.txt"),
                phenos = paste0(indir, "phenos_cleaned_", args$cohort, ".txt"))
outfile <- paste0(outdir, "resid_regress_", label, "_", args$cohort, ".txt")
bookpre <- paste0(outdir, "resid_regress_", label, "_", args$cohort)
if (file.exists(outfile)) stop("Output file exists, refusing to overwrite: ", outfile)

###### FUNCTIONS ######
## exact replica of the current GitHub version of jaffelab::cleaningY
## (missing values are set to 0 before the coefficients are estimated)
CleaningY_zero <- function(y, mod, P) {
    stopifnot(P <= ncol(mod))
    stopifnot(qr(mod)$rank == ncol(mod))
    Hat <- solve(t(mod) %*% mod) %*% t(mod)
    ty <- t(y)
    ty[is.na(ty)] <- 0
    beta <- Hat %*% ty
    y - t(as.matrix(mod[, -c(seq_len(P))]) %*% beta[-seq_len(P), ])
}

## same cleaning, but the coefficients of each protein are estimated from its observed
## samples only; missing values stay missing. Without missing values it equals CleaningY_zero.
## Returns the matrix and the number of proteins whose observed samples left a column unidentified
## (its coefficient is then set to 0).
CleaningY_naAware <- function(y, mod, P) {
    stopifnot(P <= ncol(mod))
    nm <- ncol(mod); out <- y; n_deficient <- 0L
    patt <- apply(is.na(y), 1, function(r) paste(which(r), collapse = ","))
    for (pt in unique(patt)) {
        rows <- which(patt == pt)
        obs  <- if (pt == "") seq_len(ncol(y)) else setdiff(seq_len(ncol(y)), as.integer(strsplit(pt, ",")[[1]]))
        qrM  <- qr(mod[obs, , drop = FALSE])
        if (qrM$rank < nm) n_deficient <- n_deficient + length(rows)
        beta <- qr.coef(qrM, t(y[rows, obs, drop = FALSE]))
        beta[is.na(beta)] <- 0
        fitn <- mod[obs, -seq_len(P), drop = FALSE] %*% beta[-seq_len(P), , drop = FALSE]
        out[rows, obs] <- y[rows, obs, drop = FALSE] - t(fitn)
    }
    list(cleaned = out, n_deficient = n_deficient)
}

###### MAIN (data preparation as in script V, L71-L118 without the SVA part) ######
phenos.orig <- fread(infiles$phenos)
prot.orig   <- fread(infiles$prot)
phenos.orig[, Batch := paste0("Batch", Batch)]
prot.sampfiltered <- prot.orig[, c("protein", phenos.orig$protsample), with = FALSE]
obscounts <- apply(prot.sampfiltered[, !c("protein")], 1, function(x) { sum(!is.na(x)) })
prot.countfiltered <- prot.sampfiltered[obscounts >= 50, ]

nomissing <- phenos.orig[, apply(.SD, 1, function(x) { !any(is.na(x)) }),
                         .SDcol = c(nuisance_vars, npsvars)]
phenos <- phenos.orig[nomissing, ]
phenos <- phenos[protsample %in% colnames(prot.countfiltered), ]
testsamples <- as.character(phenos$protsample)
edata <- as.matrix(prot.countfiltered[, ..testsamples])
rownames(edata) <- prot.countfiltered$protein
stopifnot(all(phenos$protsample == colnames(edata)))

modelvars <- if (args$protect == "yes") c(npsvars, nuisance_vars) else nuisance_vars
np        <- if (args$protect == "yes") length(npsvars) + 1 else 1
mod       <- model.matrix(as.formula(paste("~ ", paste0(modelvars, collapse = "+"))), data = phenos)

###### Code 3b: bookkeeping (aggregate numbers only) ######
nbatch <- length(unique(phenos$Batch))
nuis_mat <- model.matrix(as.formula(paste("~ ", paste0(nuisance_vars, collapse = "+"))), data = phenos)
pW <- ncol(nuis_mat) - 1                                   # nuisance columns without the intercept
na_per_prot <- rowSums(is.na(edata))
book <- c(cohort = args$cohort,
          n_subjects_in_phenos_file = nrow(phenos.orig),
          n_subjects_dropped_incomplete = sum(!nomissing),
          n_subjects_used = nrow(phenos),
          n_proteins_in_file = nrow(prot.orig),
          n_proteins_ge50_obs = nrow(edata),
          n_batches = nbatch, p_W = pW,
          protected_columns_P = np, cleaning_model_columns = ncol(mod),
          n_proteins_no_NA = sum(na_per_prot == 0),
          n_proteins_with_NA = sum(na_per_prot > 0),
          frac_missing_entries = mean(is.na(edata)),
          NA_per_protein_median = median(na_per_prot),
          NA_per_protein_max = max(na_per_prot))
fwrite(data.table(item = names(book), value = as.character(book)),
       file = paste0(bookpre, "_code3_bookkeeping.txt"), sep = "\t", quote = FALSE)
fwrite(data.table(n_NA_category = c("0", "1-5", "6-10", ">10"),
                  n_proteins = c(sum(na_per_prot == 0), sum(na_per_prot %in% 1:5),
                                 sum(na_per_prot %in% 6:10), sum(na_per_prot > 10))),
       file = paste0(bookpre, "_na_per_protein.txt"), sep = "\t", quote = FALSE)

## Batch x NPS association: cases per batch (cells < 5 shown as "<5"), Monte-Carlo chi-square p, Cramer's V
set.seed(1)
bn <- lapply(npsvars, function(v) {
    tab <- table(phenos$Batch, phenos[[v]])
    st  <- suppressWarnings(chisq.test(tab)$statistic)
    pm  <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 10000)$p.value)
    cv  <- sqrt(as.numeric(st) / (sum(tab) * min(nrow(tab) - 1, ncol(tab) - 1)))
    cases <- if ("1" %in% colnames(tab)) tab[, "1"] else rep(0, nrow(tab))
    data.table(nps = v, batch = rownames(tab), n = rowSums(tab),
               cases = ifelse(cases < 5, "<5", as.character(cases)),
               mc_chisq_p = pm, cramers_V = cv)
})
fwrite(rbindlist(bn), file = paste0(bookpre, "_batch_by_nps.txt"), sep = "\t", quote = FALSE)
print(book)
if (args$bookkeeping_only) { print("bookkeeping written; stopping (--bookkeeping_only)"); quit(save = "no") }

###### Cleaning ######
if (args$na == "zero") {
    cleaned <- CleaningY_zero(edata, mod, P = np)
    extra <- NULL
} else {
    res <- CleaningY_naAware(edata, mod, P = np)
    cleaned <- res$cleaned
    extra <- c(n_proteins_with_unidentified_column = res$n_deficient)
}

## compare with the existing file of script V (zero versions only) and with the zero replica (aware version)
refname <- paste0(outdir, "resid_regress_BatchPmiAgeMsex", ifelse(args$protect == "no", "_noProtect", ""),
                  "_", args$cohort, ".txt")
cmp <- NULL
if (file.exists(refname)) {
    ref <- fread(refname)
    rm_ <- as.matrix(ref[, -1]); rownames(rm_) <- ref[[1]]
    common <- intersect(rownames(rm_), rownames(cleaned))
    sm <- intersect(colnames(rm_), colnames(cleaned))
    d <- abs(rm_[common, sm] - cleaned[common, sm])
    cmp <- c(compared_with = refname, n_proteins_compared = length(common),
             max_abs_diff = max(d, na.rm = TRUE),
             n_proteins_diff_gt_1e_6 = sum(apply(d, 1, max, na.rm = TRUE) > 1e-6))
    if (args$na == "zero")
        print(paste0("Replica of cleaningY vs existing file: max |diff| = ", signif(as.numeric(cmp["max_abs_diff"]), 3),
                     " (should be ~1e-8 if the jaffelab version used sets missing values to 0)"))
}
fwrite(data.table(item = c(names(extra), names(cmp)), value = as.character(c(extra, cmp))),
       file = paste0(bookpre, "_code3_check.txt"), sep = "\t", quote = FALSE)

### Write (same format as script V)
cleaned_forprint <- as.data.table(cleaned, keep.rownames = "protein")
write.table(cleaned_forprint, file = outfile, row.names = FALSE, quote = FALSE, sep = "\t")
print(paste0("Written: ", outfile))
