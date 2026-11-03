###############################################################################
## fl_functions_r3.r                                             (Round 3, design 5d)
##
## Freedman-Lane permutation on the PROTEIN side.
##
## Idea (in words):  write the raw expression of one protein as
##        x_i  =  f(Z_i)  +  e_i ,        Z = Batch, PMI, age, sex
## where f(Z_i) is the part explained by the nuisance variables and e_i is the rest.
## If the protein is unrelated to the NPS outcomes, e_i carries no information about the NPS vector
## of subject i. So we may shuffle the e_i among subjects (inside the same Batch) and rebuild
##        x*_i  =  f(Z_i)  +  e_donor(i) .
## The NPS outcomes are NOT touched: their relation to Batch, PMI, age and sex stays as in the real data.
## Then the whole pipeline is repeated on x*: cleaning of script V (NPS protected), z-scaling, GEE, tests.
##
## Missing values: the missing value stays together with its residual (an NA moves with e_i).
## The same donor index is used for all proteins of a replicate, so the correlation between proteins is kept.
##
## Requires cleaning_functions_r2.r (nuisance_vars, MakeDonor) to be sourced first.
###############################################################################

## fits f(Z) of every protein by least squares on its OBSERVED samples (Z = nuisance_vars),
## and returns the fitted values and the residuals (NA where the expression is missing)
##   edata  : proteins x samples, raw log2 expression (as LoadCleaningInputs returns it)
##   phenos : table in the column order of edata, with the columns nuisance_vars
FLSetup <- function(edata, phenos, nuisance = nuisance_vars) {
    stopifnot(ncol(edata) == nrow(phenos))
    modZ <- model.matrix(as.formula(paste("~", paste(nuisance, collapse = "+"))), data = phenos)
    fitted <- matrix(NA_real_, nrow(edata), ncol(edata), dimnames = dimnames(edata))
    for (j in seq_len(nrow(edata))) {
        y  <- edata[j, ]
        ok <- !is.na(y)
        cf <- rep(0, ncol(modZ))
        if (sum(ok) > ncol(modZ)) {
            f  <- lm.fit(modZ[ok, , drop = FALSE], y[ok])
            cf <- f$coefficients
            cf[is.na(cf)] <- 0                      # aliased column (a Batch with too few samples): coefficient 0
        } else if (sum(ok) > 0) cf[1] <- mean(y[ok])
        fitted[j, ] <- as.numeric(modZ %*% cf)
    }
    list(fitted = fitted, resid = edata - fitted, rank_Z = qr(modZ)$rank, n_par_Z = ncol(modZ))
}

## rebuilds the expression with the residuals of subject donor[i] given to subject i
## (donor = identity gives back the real expression)
FLBuild <- function(fl, donor) {
    stopifnot(length(donor) == ncol(fl$fitted))
    x <- fl$fitted + fl$resid[, donor, drop = FALSE]
    dimnames(x) <- dimnames(fl$fitted)              # sample names stay with the SUBJECT, not with the residual
    x
}
