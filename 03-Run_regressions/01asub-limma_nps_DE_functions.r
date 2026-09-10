



Adjust_p <- function(p) {
    ## Multiple-testing adjustment. Storey q-values by default, with an
    ## explicit Benjamini-Hochberg fallback if pi0 estimation fails.
    tryCatch(qvalue(p)$qvalues,
             error = function(e) {
                 warning("qvalue() failed (", conditionMessage(e),
                         "); using Benjamini-Hochberg instead.")
                 p.adjust(p, method = "BH")
             })
}


Annotate_result <- function(dt, ...) {
    ## Attach run metadata and the FDR column to a topTable result.
    meta <- list(...)
    for (nm in names(meta)) set(dt, j = nm, value = meta[[nm]])
    set(dt, j = "fdr", value = Adjust_p(dt[["P.Value"]]))
    dt[]
}


Build_design <- function(predictors, covariates, pheno_dt, sample_ids,
                         svs = NULL) {
    ## Build a design matrix from predictors + covariates (+ optional SVs),
    ## drop constant/aliased columns, keep the intercept.
    f <- as.formula(paste("~", paste(c(predictors, covariates),
                                     collapse = " + ")))
    design_dt <- data.table(model.matrix(f, data = pheno_dt))
    if (!is.null(svs)) design_dt <- cbind(design_dt, as.data.table(svs))
  design <- as.matrix(design_dt)
  rownames(design) <- sample_ids
  colnames(design) <- make.names(colnames(design))
  informative <- apply(design, 2, function(x) length(unique(x)) > 1L)
  informative[1] <- TRUE                                # keep intercept
  design[, informative, drop = FALSE]
}

Fit_limma <- function(edat, design) {
  lmFit(edat, design) |>
    eBayes(robust = args$robust, trend = args$trend)
}

Estimate_svs <- function(protected_vars, pheno_dt, sample_ids) {
  nullmodel <- paste("~", paste(covars, collapse = " + "))
  fullmodel <- paste("~", paste(c(protected_vars, covars), collapse = " + "))
  message("  SVA null model: ", nullmodel)
  message("  SVA full model: ", fullmodel)

  mod0 <- model.matrix(as.formula(nullmodel), data = pheno_dt)
  mod  <- model.matrix(as.formula(fullmodel), data = pheno_dt)

  edat <- expr[complete_rows, sample_ids, drop = FALSE]
  if (nrow(edat) < 100L) {
    warning("Fewer than 100 complete proteins; skipping SVA.")
    return(NULL)
  }

  n_sv <- num.sv(edat, mod, method = "be")
  if (n_sv < 1L) {
    message("  no surrogate variables estimated")
    return(NULL)
  }
  message("  estimating ", n_sv, " surrogate variable(s)")
  sv <- sva(edat, mod, mod0, n.sv = n_sv)$sv
  colnames(sv) <- paste0("SV", seq_len(ncol(sv)))
  sv
}


Split_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

## Storey pi0, returned as NA rather than an error when estimation fails.
Estimate_pi0 <- function(p) {
  tryCatch(qvalue(p)$pi0, error = function(e) NA_real_)
}


Pvalue_diagnostics <- function(fit, p, domain, n_samples, n_cases, n_svs) {
    ## Reported numerically per domain: the proportion of p < 0.05, the
    ## observed/expected ratio at p < 0.05 (1 under the null), the Storey pi0
    ## estimate (proportion of true nulls), and limma's prior degrees of
    ## freedom, which is Inf when variance moderation collapses to a single
    ## common variance.

    p <- p[is.finite(p)]
    data.table(
        cohort          = args$cohort,
        domain          = domain,
        model           = "per_domain",
        n_samples       = n_samples,
    n_cases         = n_cases,
    n_svs           = n_svs,
    n_proteins      = length(p),
    prop_p_lessthan_05    = mean(p < 0.05),
    obs_exp_p_lessthan_05 = mean(p < 0.05) / 0.05,
    prop_p_lesstthan_001   = mean(p < 0.001),
    median_p        = median(p),
    ks_unif_p       = tryCatch(suppressWarnings(ks.test(p, "punif")$p.value),
                               error = function(e) NA_real_),
    pi0             = Estimate_pi0(p),
    df_prior        = if (length(fit$df.prior) == 1L) fit$df.prior
                      else median(fit$df.prior),
    s2_prior        = if (length(fit$s2.prior) == 1L) fit$s2.prior
                      else median(fit$s2.prior)
  )
}



Write_diagnostic_plot <- function(fit, p, domain, n_samples, n_cases,
                                  n_svs, path) {
    ## Write a single two-panel PNG per domain: p-value histogram (left) and
## plotSA (right).
## plotSA =  limma's residual standard deviation versus average
##     log2 abundance, with the empirical Bayes prior overlaid. A clear
##     downward or upward trend argues for eBayes(trend = TRUE); points
##     flagged far above the prior are the hypervariable proteins that
##     eBayes(robust = TRUE) protects against.
##
  png(path, width = args$plot_width, height = args$plot_height,
      units = "in", res = 200)
  on.exit(dev.off(), add = TRUE)

  op <- par(mfrow = c(1, 2), mar = c(4.6, 4.6, 4.2, 1.4), mgp = c(2.6, 0.8, 0),
            cex.main = 1.0, cex.lab = 0.95)
  on.exit(par(op), add = TRUE)

  p_ok <- p[is.finite(p)]
  hist(p_ok, breaks = seq(0, 1, by = 0.02), col = "grey80", border = "white",
       xlab = "unadjusted p-value", ylab = "number of proteins",
       main = paste0(args$cohort, ": ", domain, "\n",
                     "n = ", n_samples, " (cases = ", n_cases,
                     "), SVs = ", n_svs))
  abline(h = length(p_ok) * 0.02, col = "red", lty = 2, lwd = 1.5)
  legend("topright", bty = "n", lty = 2, col = "red", lwd = 1.5,
         legend = "uniform (null) expectation", cex = 0.85)

  plotSA(fit, main = paste0("plotSA: ", domain,
                            "\ndf.prior = ",
                            signif(median(fit$df.prior), 3)),
         xlab = "average log2 abundance",
         ylab = expression(sqrt(sigma)))
  invisible(path)
}
