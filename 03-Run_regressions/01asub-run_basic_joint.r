









### Joint model with all NPS domains as predictors ###
## Estimates each domain's effect conditional on the others; the moderated F over
## all NPS coefficients is the joint "any domain" test.
message("=== Joint model")
### The test requires complete data for all NPS
complete_nps <- pdat[, Reduce(`&`, lapply(.SD, function(v) !is.na(v))),
                     .SDcols = npsvars.bin]
sub_j     <- pdat[complete_nps]
samples_j <- sub_j[[sample_id_col]]
message("Samples with complete data on all ", length(npsvars.bin),
        " domains: ", length(samples_j))

### Get some info
case_counts <- vapply(as.character(npsvars.bin), function(d) sum(sub_j[[d]]), numeric(1))
print(data.table(domain = npsvars.bin,
                 n_cases = case_counts,
                 n_controls = length(samples_j) - case_counts))

dom_j <- names(case_counts)[case_counts >= args$min_group &
                              (length(samples_j) - case_counts) >=
                              args$min_group]
dropped <- setdiff(npsvars.bin, dom_j)
if (length(dropped) > 0L) {
  message("Domains dropped from joint model: ",
          paste(dropped, collapse = ", "))
}

### Estimate SVs and define model
edat_j   <- expr[, samples_j, drop = FALSE]
sv_j     <- Estimate_svs(dom_j, sub_j, samples_j)
design_j <- Build_design(dom_j, covars, sub_j, samples_j, sv_j)

nps_coefs <- intersect(make.names(dom_j), colnames(design_j))
message("NPS coefficients in joint model: ",
        paste(nps_coefs, collapse = ", "))

## Correlated NPS domains can make the design rank deficient
qr_rank <- qr(design_j)$rank
if (qr_rank < ncol(design_j)) {
    warning("Joint design is rank deficient (rank ", qr_rank, " < ",
            ncol(design_j), "); check collinearity among NPS domains.")
}

fit_j <- Fit_limma(edat_j, design_j)


### Extract conditionally independent per-domain effects
joint_per_domain <- lapply(nps_coefs, function(cf) {
    topTable(fit_j, coef = cf, number = Inf, sort.by = "none") |>
    as.data.table(keep.rownames = "protein") |>
    Annotate_result(domain = cf,
                    cohort = args$cohort,
                    model = "joint",
                    n_samples = length(samples_j),
                    n_svs = if (is.null(sv_j)) 0L else ncol(sv_j))
}) |>
rbindlist(use.names = TRUE, fill = TRUE)
setorder(joint_per_domain, domain, P.Value)
joint_per_domain[fdr < 0.05, sig.qvalue05 := "sig"]


### Extract moderated F across all NPS coefficients
joint_F <- topTable(fit_j, coef = nps_coefs, number = Inf,
                    sort.by = "none") |>
as.data.table(keep.rownames = "protein") |>
Annotate_result(cohort = args$cohort,
                model = "joint_F",
                n_domains = length(nps_coefs),
                n_samples = length(samples_j))
setorder(joint_F, P.Value)
joint_F[fdr < 0.05, sig.qvalue := "sig"]

message("Joint F-test proteins at FDR < 0.05: ", joint_F[fdr < 0.05, .N])
