###############################################################################
## cleaning_functions_r2.r                                           (Round 2)
##
## Functions that repeat the covariate-regression ("cleaning") step of script V
## (01a-estimate_SVs_make_adjusted_residuals.r) inside R, so that the permutation
## null of 03a_perm_null_r2.r can repeat the cleaning with permuted NPS outcomes
## (designs 5b and 5c).
##
## Everything here is taken from 01a_code3_cleaning_variants.r (Round 1), whose
## replica of jaffelab::cleaningY reproduced your existing residual file to 5e-14
## in Emory. The surrogate variables are not needed (option ignoreSVs / noProtect),
## so sva is not run.
##
## Requires: data.table; project_constants.r (npsvars.bin) sourced before.
###############################################################################

nuisance_vars <- c("Batch", "pmi", "age_death", "msex")     # as in script V

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
## samples only; missing values stay missing (without missing values: identical to the above)
CleaningY_naAware <- function(y, mod, P) {
    stopifnot(P <= ncol(mod))
    nm <- ncol(mod); out <- y
    patt <- apply(is.na(y), 1, function(r) paste(which(r), collapse = ","))
    for (pt in unique(patt)) {
        rows <- which(patt == pt)
        obs  <- if (pt == "") seq_len(ncol(y)) else setdiff(seq_len(ncol(y)), as.integer(strsplit(pt, ",")[[1]]))
        qrM  <- qr(mod[obs, , drop = FALSE])
        beta <- qr.coef(qrM, t(y[rows, obs, drop = FALSE]))
        beta[is.na(beta)] <- 0
        fitn <- mod[obs, -seq_len(P), drop = FALSE] %*% beta[-seq_len(P), , drop = FALSE]
        out[rows, obs] <- y[rows, obs, drop = FALSE] - t(fitn)
    }
    out
}


## Read the raw expression and the phenotype file and prepare them exactly as script V does
## (script V L71-L118; same as 01a_code3_cleaning_variants.r)
## returns list(edata = proteins x samples (raw, log2), phenos = data.table in the column order of edata)
LoadCleaningInputs <- function(cohort, npsvars = npsvars.bin) {
    indir <- "2_Pipeline/01-Explore_data/"
    phenos.orig <- fread(paste0(indir, "phenos_cleaned_", cohort, ".txt"))
    prot.orig   <- fread(paste0(indir, "prot.pcafiltered.log2norm.txt"))
    phenos.orig[, Batch := paste0("Batch", Batch)]
    prot.sampfiltered <- prot.orig[, c("protein", phenos.orig$protsample), with = FALSE]
    obscounts <- apply(prot.sampfiltered[, !c("protein")], 1, function(x) sum(!is.na(x)))
    prot.countfiltered <- prot.sampfiltered[obscounts >= 50, ]
    nomissing <- phenos.orig[, apply(.SD, 1, function(x) !any(is.na(x))),
                             .SDcols = c(nuisance_vars, npsvars)]
    phenos <- phenos.orig[nomissing, ]
    phenos <- phenos[protsample %in% colnames(prot.countfiltered), ]
    testsamples <- as.character(phenos$protsample)
    edata <- as.matrix(prot.countfiltered[, ..testsamples])
    rownames(edata) <- prot.countfiltered$protein
    stopifnot(all(phenos$protsample == colnames(edata)))
    list(edata = edata, phenos = phenos)
}

## One cleaning of the proteins `proteins` (names; NULL = all)
##   protect = "yes": intercept + the 7 NPS outcomes protected (option ignoreSVs, used in the analysis)
##   protect = "no" : only the intercept protected (option noProtect)
##   na      = "zero" (what jaffelab does) or "aware"
## phenos must have the columns nuisance_vars and npsvars (in the order of the columns of edata)
CleanOnce <- function(edata, phenos, protect = "yes", na = "zero", proteins = NULL,
                      npsvars = npsvars.bin) {
    stopifnot(protect %in% c("yes", "no"), na %in% c("zero", "aware"))
    if (!is.null(proteins)) edata <- edata[proteins, , drop = FALSE]
    modelvars <- if (protect == "yes") c(npsvars, nuisance_vars) else nuisance_vars
    np  <- if (protect == "yes") length(npsvars) + 1 else 1
    mod <- model.matrix(as.formula(paste("~ ", paste0(modelvars, collapse = "+"))), data = phenos)
    if (na == "zero") CleaningY_zero(edata, mod, P = np) else CleaningY_naAware(edata, mod, P = np)
}


## Permutation of the NPS vectors among subjects (all 7 domains together, as in Code 5a)
##   design "all"   : among all subjects (designs 5a, 5b)
##   design "batch" : among subjects of the same Batch (design 5c); needs the column Batch
## returns the donor index: the NPS vector of subject i becomes the vector of subject donor[i]
MakeDonor <- function(phenos, design = "all") {
    n <- nrow(phenos)
    if (design == "all") return(sample.int(n))
    donor <- seq_len(n)
    for (b in unique(phenos$Batch)) {
        idx <- which(phenos$Batch == b)
        donor[idx] <- if (length(idx) > 1) idx[sample.int(length(idx))] else idx
    }
    donor
}

## apply a donor index to the NPS columns of a table (returns a copy)
PermuteNPS <- function(phenos, donor, npsvars = npsvars.bin) {
    out <- copy(phenos)
    for (v in npsvars) set(out, j = v, value = phenos[[v]][donor])
    out
}


## z-scaling or rank-based inverse-normal transform of each row (03a L94-105 for "z")
## keeps the proteins with more than 50 non-missing values, as PrepareData of Round 1 does
ScaleRows <- function(resid, xtransform = "z") {
    tokeep <- apply(resid, 1, function(x) sum(!is.na(x)) > 50)
    x <- resid[tokeep, , drop = FALSE]
    out <- t(apply(x, 1, function(v) {
        if (xtransform == "rankint") {
            ok <- !is.na(v); r <- rank(v[ok], ties.method = "average")
            w <- rep(NA_real_, length(v)); w[ok] <- qnorm((r - 0.5) / sum(ok))
            as.numeric(scale(w))
        } else as.numeric(scale(v))
    }))
    rownames(out) <- rownames(x); colnames(out) <- colnames(x)
    out
}
