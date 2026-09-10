#!/usr/bin/env Rscript
## =====================================================================
## limma differential protein expression (DE) for neuropsychiatric
## symptom (NPS/BPSD) domains in TMT brain proteomics.
##
## Analysis A: one limma model per NPS domain (moderated t).
## Analysis B: one joint model with all NPS domains as simultaneous
##             predictors, giving (B1) conditionally independent
##             per-domain moderated t statistics and (B2) a moderated
##             F test across all NPS coefficients.
## =====================================================================

### SETUP ###
rm(list = ls())
library(argparse)
library(data.table)
library(limma)
library(sva)
library(qvalue)

source("~/comsv/Svattathil_Library/svattathil_functions.r")


### CONSTANTS
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")

covars_other_than_SVs <- c("Batch", "pmi", "age_death", "msex")
fac_covars <- c("Batch")
max_missing <- 0.5  ## max missing fraction allowed for proteins


### FUNCTIONS ###
source("1_Code/03-Run_regressions/01asub-limma_nps_DE_functions.r")


### Command line interface
parser <- ArgumentParser(
  description = "limma DE of log2 protein abundance on NPS domains"
)
parser$add_argument("--cohort", default = "OHSU",
                    help = "Cohort label used in output file names.")
parser$add_argument("--robust", action = "store_true", default = TRUE,
                    help = "Use eBayes(robust = TRUE).")
parser$add_argument("--trend", action = "store_true", default = FALSE,
                    help = "Use eBayes(trend = TRUE).")
parser$add_argument("--min-group", dest = "min_group", type = "integer",
                    default = 5L,
                    help = "Minimum cases and controls per domain [5].")
parser$add_argument("--seed", type = "integer", default = 20260101L)
parser$add_argument("--no-diagnostics", dest = "no_diagnostics",
                    action = "store_true", default = FALSE,
                    help = "Skip per-domain p-value / plotSA diagnostics.")
parser$add_argument("--plot-width", dest = "plot_width", type = "double",
                    default = 11,
                    help = "Diagnostic PNG width in inches [11].")
parser$add_argument("--plot-height", dest = "plot_height", type = "double",
                    default = 5,
                    help = "Diagnostic PNG height in inches [5].")

args <- parser$parse_args()

## Process some of the args
set.seed(args$seed)


### Files that exist
infiles <- list(
    phenos = paste0("2_Pipeline/02-Prepare_analysis_data/phenos_cleaned_svs_",
                    args$cohort, ".txt"),
    prot = "2_Pipeline/01-Explore_data/prot.pcafiltered.log2norm.txt"
)


### Files to be created
outdir <- "2_Pipeline/03-Run_regressions/Prot_perDomain/"
MyMkdir(outdir)

diag_dir <- file.path(outdir, "Diagnostics/")
if (!args$no_diagnostics) { MyMkdir(diag_dir) }


outfiles <- list(
    limma_basic     = paste0(args$cohort, "_limma_basic.tsv"),
    joint_perdomain = paste0(args$cohort, "_limma_joint_perDomain.tsv"),
    joint_ftest     = paste0(args$cohort, "_limma_joint_Ftest.tsv"),
    summary_dt      = paste0(args$cohort, "_limma_basic_summary.tsv"),
    diagnostics     = paste0(args$cohort, "_basic_diagnostics.tsv"),
    session_info    = "sessionInfo.txt"
    )

## This file is per domain so it is defined using a function
outfiles[["SA_plot"]] <- function(adomain) {
    paste0(diag_dir, paste0(args$cohort, "_", make.names(adomain), "_pvalue_plotSA.png"))
}


### MAIN ###
### Read data
pheno <- fread(infiles$pheno)
prot <- fread(infiles$prot)

sample_id_col <- "protsample"
protein_id_col <- names(prot)[1]


### Align samples and filter
expr_all <- as.matrix(prot[, !protein_id_col, with = FALSE])
rownames(expr_all) <- as.character(prot[[protein_id_col]])
storage.mode(expr_all) <- "double"

pheno[, (sample_id_col) := as.character(get(sample_id_col))]
setkeyv(pheno, sample_id_col)

shared <- intersect(colnames(expr_all), pheno[[sample_id_col]])
if (length(shared) == 0L) {
  stop("No sample IDs shared between --pheno and --prot.")
}
message("Samples shared between tables: ", length(shared))

expr_all <- expr_all[, shared, drop = FALSE]
pdat     <- pheno[J(shared)]                 # reordered to match columns
stopifnot(identical(as.character(pdat[[sample_id_col]]), colnames(expr_all)))


### Set factor covars
### and get SV columns
for (nm in fac_covars) pdat[, (nm) := factor(get(nm))]
svcols <- grep("SV", names(pdat), value = TRUE)
n_svs <- length(svcols)
covars <- c(covars_other_than_SVs, svcols)


### limma drops samples with missing covariates from every protein, so
### remove them once, up front, and report the count
covar_ok <- pdat[, Reduce(`&`, lapply(.SD, function(v) !is.na(v))),
                 .SDcols = covars]
nps_any  <- pdat[, Reduce(`|`, lapply(.SD, function(v) !is.na(v))),
                 .SDcols = npsvars.bin]
keep_samp <- covar_ok & nps_any
message("Samples dropped (missing covariate or no NPS observation): ",
        sum(!keep_samp))

expr_all <- expr_all[, keep_samp, drop = FALSE]
pdat     <- pdat[keep_samp]
for (nm in fac_covars) pdat[, (nm) := droplevels(get(nm))]


### Filter to proteins that pass missingness threshold
miss_frac <- rowMeans(is.na(expr_all))
keep_prot <- miss_frac <= max_missing
message("Proteins retained at missingness <= ", max_missing, ": ",
        sum(keep_prot), " of ", length(keep_prot))
expr <- expr_all[keep_prot, , drop = FALSE]


#### Analysis A: one model per NPS domain ####
## Set up listobjects to hold results per domain
per_domain <- vector("list", length(npsvars.bin))
names(per_domain) <- npsvars.bin

diag_list  <- vector("list", length(npsvars.bin))
names(diag_list) <- npsvars.bin


### Run for each domain in turn
for (dom in npsvars.bin) {
    ## Prepare data
    message("=== Per-domain model: ", dom)
    sub_pdat <- pdat[!is.na(get(dom))]   ## subset of donors with non-missing NPS for this domain
    samples <- sub_pdat[[sample_id_col]]
    n_case <- sum(sub_pdat[[dom]])
    n_ctrl <- length(samples) - n_case
    message("  n = ", length(samples), " (cases = ", n_case, ", controls = ", n_ctrl, ")")

    if (n_case < args$min_group || n_ctrl < args$min_group) {
        warning("Fewer than ", args$min_group, " cases or controls for ", dom, "; skipping.")
        next
    }

    edat <- expr[, samples, drop = FALSE]

    design <- Build_design(predictors = dom, covariates = covars,
                           pheno_dt = sub_pdat, sample_ids = samples)

    coef_name <- make.names(dom)
    if (!coef_name %in% colnames(design)) {
        warning("Coefficient ", coef_name, " not estimable; skipping ", dom)
        next
    }

    ## Fit model
    fit <- Fit_limma(edat, design)


    ## Extract results
    res <- fit |>
    topTable(coef = coef_name, number = Inf, confint = TRUE, sort.by = "none") |>
    as.data.table(keep.rownames = "protein")

    ## Add unmoderated and moderated SE columns
    res[, SE_unmod := (fit$sigma * fit$stdev.unscaled[, coef_name])]
    res[, SE_mod := (sqrt(fit$s2.post) * fit$stdev.unscaled[, coef_name])]
    res[, df_residual := fit$df.residual]
    res[, df_total := fit$df.residual + fit$df.prior]

    Annotate_result(dt = res,
                    domain = dom,
                    cohort = args$cohort,
                    model = "per_domain",
                    n_samples = length(samples),
                    n_cases = n_case,
                    n_svs = n_svs)
    setorder(res, P.Value)

    res[fdr < 0.05, sig.qvalue05 := "sig"]
    message("  proteins at FDR < 0.05: ", res[fdr < 0.05, .N])

    ## Save to cross-nps list
    per_domain[[dom]] <- res

    ## Run diagnostics and draw plot to file
    if (!args$no_diagnostics) {
        diag_list[[dom]] <- Pvalue_diagnostics(fit, res[["P.Value"]], dom,
                                               length(samples), n_case,
                                               n_svs)
        Write_diagnostic_plot(fit, res[["P.Value"]], dom, length(samples), n_case,
                              n_svs, outfiles$SA_plot(dom)
                              )
        message("  diagnostics: prop p < 0.05 = ",
                signif(diag_list[[dom]]$prop_p_lessthan_05, 3),
                ", pi0 = ", signif(diag_list[[dom]]$pi0, 3),
                ", df.prior = ", signif(diag_list[[dom]]$df_prior, 3))
    }
}


### Gather results across domains
per_domain_dt <- rbindlist(per_domain, use.names = TRUE, fill = TRUE)
setcolorder(per_domain_dt, c("cohort", "domain", "protein", "logFC", "CI.L", "CI.R",
                             "SE_unmod", "SE_mod",
                             "AveExpr", "t",
                             "P.Value", "fdr"))
per_domain_dt[, model := NULL] ## remove extra column


### Gather diagnostics across domains
if (!args$no_diagnostics) {
    diag_dt <- rbindlist(diag_list, use.names = TRUE, fill = TRUE)
    if (nrow(diag_dt) > 0L) {
        setorder(diag_dt, domain)
        print(diag_dt)

        ## Flag domains whose p-value distribution looks misspecified.
        flagged <- diag_dt[obs_exp_p_lessthan_05 > 2 | prop_p_lessthan_05 < 0.02 |
                           (!is.na(pi0) & pi0 > 1.05)]
        if (nrow(flagged) > 0L) {
            message("Domains with atypical p-value distributions (inspect the ",
                    "PNGs before interpreting FDR): ",
                    paste(flagged$domain, collapse = ", "))
        }
    }
}


### 4. Analysis B: joint model with all NPS domains as predictors
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
                              (length(samples_j) - case_counts) >= args$min_group]
dropped <- setdiff(npsvars.bin, dom_j)
if (length(dropped) > 0L) {
  message("Domains dropped from joint model: ", paste(dropped, collapse = ", "))
}


### Define data and model
edat_j   <- expr[, samples_j, drop = FALSE]
design_j <- Build_design(predictors = dom_j, covariates = covars,
                           pheno_dt = sub_j, sample_ids = samples_j)

nps_coefs <- intersect(make.names(dom_j), colnames(design_j))
message("NPS coefficients in joint model: ", paste(nps_coefs, collapse = ", "))

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
                    n_svs = n_svs)
}) |>
rbindlist(use.names = TRUE, fill = TRUE)
setorder(joint_per_domain, domain, P.Value)
joint_per_domain[fdr < 0.05, sig.qvalue05 := "sig"]

fwrite(joint_per_domain, file.path(outdir, outfiles$joint_perdomain), sep = "\t")


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
fwrite(joint_F, file.path(outdir, outfiles$joint_ftest), sep = "\t")


### 5. Run summary
if (nrow(per_domain_dt) > 0L) {
  summary_dt <- per_domain_dt[, .(n_proteins = .N,
                                  n_fdr05 = sum(fdr < 0.05, na.rm = TRUE),
                                  n_samples = n_samples[1],
                                  n_cases = n_cases[1],
                                  n_svs = n_svs[1]),
                              by = .(cohort, domain)]
  print(summary_dt)
  fwrite(summary_dt, file.path(outdir, outfiles$summary_dt), sep = "\t")
}

#writeLines(capture.output(sessionInfo()), file.path(outdir, outfiles$session_info))
message("Done.")


### FINISH ###
fwrite(per_domain_dt, file.path(outdir, outfiles$limma_basic), sep = "\t")

if (!args$no_diagnostics) {
        fwrite(diag_dt, file = outfiles$diagnostics, sep = "\t")
}






### Code used in exploration
if(0) {
## Check missingness v. batch
## This was in response to the warning about 'Partial NA coefficients for XX probes'
    any_missing_logi <- apply(expr, 1, function(x) { any(is.na(x)) } )

any_missing_prots <- names(any_missing_logi)[any_issing_logi]

formerge <- as.data.table(t(expr[any_missing_prots, ]), keep.rownames = "protsample")
merged <- merge(pdat, formerge, by = "protsample")

tabres <- sapply(any_missing_prots, function(aprot) { table(!is.na(merged[, get(aprot)]), merged$Batch) },
                 simplify = FALSE)

check_missing <- sapply(tabres, function(x) { any(x["TRUE", ] == 0) } )
}
