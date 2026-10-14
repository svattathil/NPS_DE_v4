###############################################################################
## diag_functions.r
##
## Helper functions shared by
##   03a_code1_diagnostics.r   (Code 1; with --residset regressSVs also Code 2)
##   03a_code5a_permutation.r  (Code 5a)
##
## Purpose: run EXACTLY the per-protein analysis of 03a-run_multivariate.r
## (z-scaled expression -> long format -> geeglm -> CR2 -> Wald_test "HTZ")
## and, in addition, record diagnostics. The diagnostics only READ the fitted
## objects; they never change the p-value that 03a would give.
##
## Place this file in 1_Code/ next to 03a-run_multivariate.r.
## (v2: PrepareData() accepts a residual file by path, needed for the Code 3 variants.)
## Requires: data.table, geepack, clubSandwich (as 03a does)
###############################################################################


###### 1. Reading and preparing the data (03a, lines 40-104, unchanged) ######

PrepareData <- function(cohort, residset, residfile = NULL) {
    indir <- "2_Pipeline/02-Prepare_analysis_data/"
    infiles <- list(phenos = paste0(indir, "phenos_cleaned_SVs_", cohort, ".txt"),
                    resid  = if (is.null(residfile)) Residfile(residset, cohort) else residfile)

    resid.all  <- fread(infiles$resid)
    phenos.all <- fread(infiles$phenos)
    n_phenos_file   <- nrow(phenos.all)
    n_resid_samples <- ncol(resid.all) - 1L        # first column = protein id
    n_proteins_file <- nrow(resid.all)

    ### Filter samples (03a L80)
    phenos.all <- phenos.all[protsample %in% colnames(resid.all), ]

    ### Restrict phenos columns (03a L83-86); modelform("agit") = "agit ~ protein"
    representative_formula <- as.formula("agit ~ protein")
    allvars <- setdiff(unique(c(npsvars.bin, all.vars(representative_formula))), "protein")
    phenos  <- phenos.all[, c("protsample", allvars), with = FALSE]

    ### Restrict resid to samples in phenos and convert to matrix (03a L89-91)
    resid <- as.matrix(resid.all[, phenos.all$protsample, with = FALSE])
    rownames(resid) <- resid.all$protein

    ### Proteins with data for more than 50 participants (03a L94-97)
    tokeep <- apply(resid, 1, function(x) { sum(!is.na(x)) > 50 })
    for_transformation <- resid[tokeep, ]

    ### Z-scale (03a L100-104)
    print("Z-scaling data for each protein")
    transformed <- t(apply(for_transformation, 1, scale))
    rownames(transformed) <- rownames(for_transformation)
    colnames(transformed) <- colnames(for_transformation)

    list(transformed = transformed,
         phenos      = phenos,
         infiles     = infiles,
         counts      = c(n_phenos_file          = n_phenos_file,
                         n_resid_samples        = n_resid_samples,
                         n_phenos_after_filter  = nrow(phenos.all),
                         n_proteins_file        = n_proteins_file,
                         n_proteins_kept        = nrow(transformed),
                         protsample_unique      = as.integer(anyDuplicated(phenos$protsample) == 0),
                         n_proteins_with_any_NA = sum(apply(for_transformation, 1,
                                                            function(x) any(is.na(x))))))
}


###### 2. The per-protein analysis of 03a (unchanged), plus diagnostics ######

DomainCounts <- function(longtest) {
    tab <- longtest[, .(n = .N, cases = sum(outcome, na.rm = TRUE)), by = nps_id]
    paste(sprintf("%s:n=%d,cases=%d", tab$nps_id, tab$n, tab$cases), collapse = "; ")
}

## x_vec      : named numeric vector, the z-scaled expression of ONE protein
##              (a row of 'transformed'), names = sample ids
## sample_ids : colnames(transformed)
## phenos     : data.table with protsample + the NPS outcome columns
##              (Code 5a passes a phenos table whose NPS vectors are permuted)
## corstr, scale_fix : defaults reproduce 03a; Code 4 will change them
RunProteinAnalysis <- function(x_vec, sample_ids, phenos, npsvars,
                               corstr = "unstructured", scale_fix = FALSE,
                               diagnostics = TRUE) {
    t0 <- proc.time()[["elapsed"]]
    domain_counts <- NA_character_
    fit_class     <- NA_character_
    warn_msgs     <- character(0)
    diag_row      <- NULL

    fit_result <- withCallingHandlers(
        tryCatch({
            ### 1. Set up data (03a L127-153)
            test <- data.frame(protein    = unlist(x_vec),
                               protsample = sample_ids)
            test <- merge(test, phenos, by = "protsample")
            test <- test[complete.cases(test), ]
            setDT(test)

            idvars <- all.vars(as.formula("protsample ~ protein"))
            longtest <- melt(test, id.vars = idvars, measure.vars = npsvars,
                             variable.name = "nps_id", value.name = "outcome")
            setDT(longtest)
            longtest[, protsample := factor(protsample)]
            longtest[, nps_id := factor(nps_id)]
            longtest[, wave := as.integer(nps_id)]
            setorder(longtest, protsample, wave)

            domain_counts <- DomainCounts(longtest)

            ### 2. GEE model (03a L166-173)
            ## (scale.fix = TRUE fixes the dispersion at its default value 1; the argument
            ##  is only passed when requested, so the default call is identical to 03a)
            if (scale_fix) {
                gee_fit_int <- geeglm(outcome ~ nps_id + protein:nps_id,
                                      id      = protsample,
                                      waves   = wave,
                                      data    = longtest,
                                      family  = binomial(link = "logit"),
                                      corstr  = corstr,
                                      scale.fix = TRUE,
                                      std.err = "san.se")
            } else {
                gee_fit_int <- geeglm(outcome ~ nps_id + protein:nps_id,
                                      id      = protsample,
                                      waves   = wave,
                                      data    = longtest,
                                      family  = binomial(link = "logit"),
                                      corstr  = corstr,
                                      std.err = "san.se")
            }
            fit_class <- paste(class(gee_fit_int), collapse = "/")

            ## Some geepack versions return linear.predictors as an N x 1 matrix,
            ## which makes clubSandwich stop with "non-conformable arrays".
            ## Flattening it changes no number; for versions where it is already a
            ## vector this line does nothing.
            gee_fit_int$linear.predictors <- as.vector(gee_fit_int$linear.predictors)

            if (anyNA(coef(gee_fit_int))) {
                stop("Aliased coefficient(s) in GEE fit (likely a constant/near-constant ",
                     "NPS domain for this protein): ",
                     paste(names(coef(gee_fit_int))[is.na(coef(gee_fit_int))], collapse = ", "))
            }

            ### 3. Interaction coefficients (03a L184-187)
            coef_vec <- coef(gee_fit_int)
            prot_terms <- grep("protein", names(coef_vec), value = TRUE)
            beta_prot  <- coef_vec[prot_terms]
            names(beta_prot) <- sub(":protein", "", sub("nps_id", "", prot_terms))

            ### 4. CR2 covariance (03a L192-201)
            V_full <- clubSandwich::vcovCR(gee_fit_int,
                                           cluster = longtest$protsample,
                                           type    = "CR2")
            V_prot <- V_full[prot_terms, prot_terms]

            ### 5. Joint Wald test, HTZ (03a L208-212)
            wald_i <- clubSandwich::Wald_test(gee_fit_int,
                         constraints = constrain_zero(prot_terms),
                         vcov        = "CR2",
                         cluster     = longtest$protsample,
                         test        = "HTZ")

            n <- length(unique(gee_fit_int$id))

            ### 6. Diagnostics (read-only; a failure here never affects the result)
            if (diagnostics) {
                diag_row <- tryCatch(
                    ProteinDiagnostics(gee_fit_int, longtest, prot_terms, V_prot,
                                       wald_i, corstr),
                    error = function(e) data.frame(diag_error = conditionMessage(e)))
            }

            list(beta = beta_prot, V = V_prot, stats = wald_i, error = NA_character_, n_g = n)

        }, error = function(e) {
            list(beta = NA, V = NA,
                 stats = data.frame(test = NA_character_, Fstat = NA_real_, df_num = NA_integer_,
                                    df_denom = NA_real_, p_val = NA_real_, sig = NA_real_),
                 error = conditionMessage(e), n_g = NA_integer_)
        }),
        warning = function(w) {
            warn_msgs <<- c(warn_msgs, conditionMessage(w))
            invokeRestart("muffleWarning")
        })

    fit_result$domain_counts <- domain_counts
    fit_result$fit_class     <- fit_class
    fit_result$diag          <- diag_row
    fit_result$warnings      <- warn_msgs
    fit_result$elapsed       <- proc.time()[["elapsed"]] - t0
    fit_result
}


###### 3. Diagnostics of one fit ######

WaldQ <- function(beta, V) {
    tryCatch(as.numeric(t(beta) %*% solve(V, beta)), error = function(e) NA_real_)
}

## working correlation matrix R-hat as used by geeglm (from geese$alpha)
WorkingCorMatrix <- function(alpha, K, corstr) {
    R <- diag(K)
    if (corstr == "independence" || length(alpha) == 0) return(R)
    if (corstr == "exchangeable") { R[upper.tri(R)] <- alpha[1]; R[lower.tri(R)] <- alpha[1]; return(R) }
    if (corstr == "unstructured" && length(alpha) == K * (K - 1) / 2) {
        ## names are "alpha.j:k", j<k
        jk <- do.call(rbind, strsplit(sub("^alpha\\.", "", names(alpha)), ":"))
        for (r in seq_len(nrow(jk))) {
            j <- as.integer(jk[r, 1]); k <- as.integer(jk[r, 2])
            R[j, k] <- R[k, j] <- alpha[r]
        }
        return(R)
    }
    matrix(NA_real_, K, K)
}

EigStats <- function(R) {
    tryCatch({
        ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
        list(min_eig = min(ev), cond = ifelse(min(ev) > 0, max(ev) / min(ev), Inf))
    }, error = function(e) list(min_eig = NA_real_, cond = NA_real_))
}

ProteinDiagnostics <- function(fit, longtest, prot_terms, V_cr2, wald_i, corstr) {
    d <- list()
    coef_vec <- coef(fit)
    beta <- coef_vec[prot_terms]
    q    <- length(prot_terms)
    K    <- nlevels(longtest$nps_id)
    n_g  <- length(unique(longtest$protsample))
    N    <- nrow(longtest)
    p    <- length(coef_vec)

    ## ---- (b) data of the protein ----
    cases <- longtest[, .(cases = sum(outcome)), by = nps_id]$cases
    d$n_g          <- n_g
    d$min_cases    <- min(cases)
    d$min_noncases <- min(n_g - cases)
    d$n_rare       <- sum(cases < 5)               # threshold 5 is a guess
    x <- longtest[wave == 1, protein]
    d$kurt_x       <- mean((x - mean(x))^4) / mean((x - mean(x))^2)^2   # normal = 3
    d$max_abs_x    <- max(abs(x))
    d$n_distinct_x <- length(unique(x))

    ## ---- (c) fit ----
    d$phi_hat    <- tryCatch(as.numeric(fit$geese$gamma)[1], error = function(e) NA_real_)
    d$geese_error <- tryCatch(as.numeric(fit$geese$error)[1], error = function(e) NA_real_)
    alpha <- fit$geese$alpha
    d$alpha_min  <- if (length(alpha)) min(alpha) else NA_real_
    d$alpha_max  <- if (length(alpha)) max(alpha) else NA_real_

    Rw <- WorkingCorMatrix(alpha, K, corstr)
    es <- EigStats(Rw)
    d$min_eig_R <- es$min_eig                       # working correlation used by geeglm
    d$cond_R    <- es$cond

    ## fitted values and Pearson residuals recomputed from the coefficients
    mm  <- model.matrix(outcome ~ nps_id + protein:nps_id, data = longtest)
    mu  <- plogis(drop(mm %*% coef_vec[colnames(mm)]))
    u   <- (longtest$outcome - mu) / sqrt(mu * (1 - mu))
    d$phi_pearson <- sum(u^2) / (N - p)             # moment estimate of the dispersion
    d$min_mu <- min(mu); d$max_mu <- max(mu)
    U <- matrix(u, ncol = K, byrow = TRUE)          # rows = subjects (longtest is sorted)
    Remp <- suppressWarnings(cor(U))
    ee <- EigStats(Remp)
    d$min_eig_Remp  <- ee$min_eig                   # empirical correlation of Pearson residuals
    d$cond_Remp     <- ee$cond
    d$max_abs_rho_emp <- tryCatch(max(abs(Remp[lower.tri(Remp)])), error = function(e) NA_real_)

    d$max_abs_beta <- max(abs(beta))
    d$max_se_cr2   <- sqrt(max(diag(V_cr2)))

    ## ---- (d) the test ----
    Fstat    <- as.numeric(wald_i$Fstat)[1]
    df_num   <- as.numeric(wald_i$df_num)[1]
    df_denom <- as.numeric(wald_i$df_denom)[1]
    d$eta_from_df <- df_denom + df_num - 1          # eta = d2 + q - 1
    ## 'delta' is a column of the Wald_test output; it is NOT eta (checked on test data), kept for reference
    d$delta_wald  <- if ("delta" %in% names(wald_i)) as.numeric(wald_i$delta)[1] else NA_real_
    d$Q_from_F    <- Fstat * df_num * d$eta_from_df / df_denom
    d$Q_cr2       <- WaldQ(beta, V_cr2)             # beta' V^-1 beta with V = CR2 block
    d$rel_diff_Q  <- abs(d$Q_from_F - d$Q_cr2) / d$Q_cr2

    ## ---- (e) alternative p-values, chi-square(q) reference ----
    d$p_chisq_cr2 <- pchisq(d$Q_cr2, df = q, lower.tail = FALSE)

    V_cr0 <- clubSandwich::vcovCR(fit, cluster = longtest$protsample, type = "CR0")[prot_terms, prot_terms]
    d$Q_cr0       <- WaldQ(beta, V_cr0)
    d$p_chisq_cr0 <- pchisq(d$Q_cr0, df = q, lower.tail = FALSE)

    idx   <- match(prot_terms, names(coef_vec))
    V_san <- fit$geese$vbeta[idx, idx]              # robust (san.se) covariance of geeglm
    d$Q_san       <- WaldQ(beta, V_san)
    d$p_chisq_san <- pchisq(d$Q_san, df = q, lower.tail = FALSE)
    d$max_se_san  <- sqrt(max(diag(V_san)))

    as.data.frame(d)
}

## one flat list per protein: the columns of 03a plus the diagnostics
ResultRow <- function(protein_name, res, n_na_resid = NA_integer_) {
    g1 <- function(v) if (is.null(v) || length(v) == 0) NA_real_ else as.numeric(v)[1]
    st <- res$stats
    row <- list(protein       = protein_name,
                Fstat         = g1(st$Fstat),
                df_num        = g1(st$df_num),
                df_denom      = g1(st$df_denom),
                P             = g1(st$p_val),
                nobs          = g1(res$n_g),
                error         = res$error,
                domain_counts = res$domain_counts,
                n_na_resid    = n_na_resid,
                n_warn        = length(res$warnings),
                warn_first    = if (length(res$warnings)) substr(res$warnings[1], 1, 200) else NA_character_,
                fit_class     = res$fit_class,
                elapsed_sec   = res$elapsed)
    if (!is.null(res$diag)) row <- c(row, as.list(res$diag))
    row
}


###### 4. Calibration summaries ######

LambdaGC <- function(p) {
    p <- p[!is.na(p)]; p[p == 0] <- 1e-300
    median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1)
}

## median over the 5% smallest p-values of (observed / expected) chi-square(1) quantile
LambdaTail <- function(p, frac = 0.05) {
    p <- sort(p[!is.na(p)]); p[p == 0] <- 1e-300
    m <- length(p); if (m < 20) return(NA_real_)
    r <- seq_len(max(1, floor(frac * m)))
    obs <- qchisq(p[r], df = 1, lower.tail = FALSE)
    exp <- qchisq((r - 0.5) / (m + 1), df = 1, lower.tail = FALSE)   # as in qqunif.plot
    median(obs / exp)
}

CohortSummary <- function(dt) {
    out <- list()
    pv <- list(HTZ = dt$P, chisq_CR2 = dt$p_chisq_cr2,
               chisq_CR0 = dt$p_chisq_cr0, chisq_san = dt$p_chisq_san)
    for (nm in names(pv)) {
        p <- pv[[nm]]
        out[[paste0("n_p_", nm)]]            <- sum(!is.na(p))
        out[[paste0("lambda_GC_", nm)]]      <- LambdaGC(p)
        out[[paste0("lambda_5pct_", nm)]]    <- LambdaTail(p)
        out[[paste0("frac_p_lt_0.05_", nm)]] <- mean(p < 0.05, na.rm = TRUE)   # 0.05 if calibrated
        out[[paste0("frac_p_lt_0.001_", nm)]] <- mean(p < 0.001, na.rm = TRUE)
    }
    qs <- c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)
    for (v in intersect(c("n_g", "df_denom", "phi_hat", "phi_pearson", "min_cases",
                          "min_eig_R", "min_eig_Remp", "kurt_x", "max_abs_x", "n_na_resid",
                          "max_abs_beta", "max_se_cr2", "max_se_san", "elapsed_sec"), names(dt))) {
        qv <- quantile(dt[[v]], qs, na.rm = TRUE)
        for (j in seq_along(qs)) out[[sprintf("%s_q%03d", v, round(100 * qs[j]))]] <- unname(qv[j])
    }
    if ("df_denom" %in% names(dt)) {
        out$frac_d2_lt_10 <- mean(dt$df_denom < 10, na.rm = TRUE)
        out$frac_d2_lt_20 <- mean(dt$df_denom < 20, na.rm = TRUE)
        for (v in intersect(c("n_g", "min_cases", "min_eig_R", "kurt_x", "max_abs_x", "n_na_resid"), names(dt)))
            out[[paste0("spearman_d2_", v)]] <-
                suppressWarnings(cor(dt$df_denom, dt[[v]], method = "spearman", use = "pairwise.complete.obs"))
        out$spearman_d2_minus_log10P <-
            suppressWarnings(cor(dt$df_denom, -log10(dt$P), method = "spearman", use = "pairwise.complete.obs"))
    }
    if ("rel_diff_Q" %in% names(dt)) out$max_rel_diff_Q <- max(dt$rel_diff_Q, na.rm = TRUE)
    if ("nobs" %in% names(dt)) out$n_failed_fits <- sum(is.na(dt$P))
    data.table(metric = names(out), value = unlist(out))
}


###### 5. Plots (base graphics) ######

QQcurve <- function(p) {
    p <- sort(p[!is.na(p)]); p[p == 0] <- 1e-300
    m <- length(p)
    list(e = -log10((seq_len(m) - 0.5) / (m + 1)), o = -log10(p), m = m)
}

AddQQ <- function(p, col, lwd = 1.5, pch = 16, cex = 0.4, type = "p") {
    cv <- QQcurve(p); if (cv$m == 0) return(invisible(NULL))
    points(cv$e, cv$o, col = col, pch = pch, cex = cex, type = type, lwd = lwd)
}

EmptyQQ <- function(plist, main) {
    mx <- max(unlist(lapply(plist, function(p) { cv <- QQcurve(p); if (cv$m) c(cv$e, cv$o) else 0 })), 1)
    plot(NA, xlim = c(0, mx), ylim = c(0, mx), xlab = expression(Expected~-log[10](p)),
         ylab = expression(Observed~-log[10](p)), main = main, cex.main = 0.9)
    abline(0, 1, col = "grey40")
}

DiagQQPlots <- function(dt, file, title = "") {
    png(file, height = 5, width = 15, units = "in", res = 200)
    par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1))

    ## panel 1: four versions of the p-value for the same fits
    pl <- list(dt$P, dt$p_chisq_cr2, dt$p_chisq_cr0, dt$p_chisq_san)
    EmptyQQ(pl, paste0(title, "\nsame fits, different reference / covariance"))
    cols <- c("red", "blue", "darkgreen", "orange")
    for (j in 1:4) AddQQ(pl[[j]], cols[j])
    legend("topleft", bty = "n", pch = 16, col = cols, cex = 0.8,
           legend = c("CR2 + HTZ (03a)", "CR2 + chi-square", "CR0 + chi-square", "san.se + chi-square"))

    ## panel 2: HTZ p-value by tertile of d2
    ok <- !is.na(dt$df_denom) & !is.na(dt$P)
    if (sum(ok) >= 30) {
        g <- ceiling(3 * rank(dt$df_denom[ok], ties.method = "first") / sum(ok))
        pl2 <- lapply(1:3, function(j) dt$P[ok][g == j])
        EmptyQQ(pl2, "CR2 + HTZ by tertile of d2")
        cols2 <- c("red", "purple", "blue")
        for (j in 1:3) AddQQ(pl2[[j]], cols2[j])
        legend("topleft", bty = "n", pch = 16, col = cols2, cex = 0.8,
               legend = sprintf("%s d2 (median %.1f)", c("low", "mid", "high"),
                                sapply(1:3, function(j) median(dt$df_denom[ok][g == j]))))
    } else plot.new()

    ## panel 3: proteins with / without a rare domain (< 5 cases)
    if ("n_rare" %in% names(dt) && sum(!is.na(dt$n_rare) & !is.na(dt$P)) >= 30) {
        rare <- dt$n_rare > 0
        pl3 <- list(dt$P[which(rare)], dt$P[which(!rare)])
        EmptyQQ(pl3, "CR2 + HTZ: rare domain (<5 cases) or not")
        AddQQ(pl3[[1]], "red"); AddQQ(pl3[[2]], "blue")
        legend("topleft", bty = "n", pch = 16, col = c("red", "blue"), cex = 0.8,
               legend = c(sprintf("any domain < 5 cases (%d proteins)", sum(rare, na.rm = TRUE)),
                          sprintf("all domains >= 5 cases (%d)", sum(!rare, na.rm = TRUE))))
    } else plot.new()
    dev.off()
}
