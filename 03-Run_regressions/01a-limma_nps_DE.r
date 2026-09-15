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
suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(limma))
suppressMessages(library(sva))
suppressMessages(library(qvalue))

source("~/comsv/Svattathil_Library/svattathil_functions.r")
source("~/comsv/Svattathil_Library/manhattan_and_qq_plot_functions.r")

### CONSTANTS
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")


### FUNCTIONS ###
source("1_Code/03-Run_regressions/01asub-limma_nps_DE_functions.r")


### Command line interface
parser <- ArgumentParser(description = "limma DE of log2 protein abundance on NPS domains")

## Control model and sample set
parser$add_argument("--run", default = "basic",
                    help = "Specify the run (basic, severity, males, females, e4_adj)")
parser$add_argument("--cohort", default = "OHSU",
                    help = "Cohort label used in output file names.")
## Control limma
parser$add_argument("--robust", action = "store_true", default = TRUE,
                    help = "Use eBayes(robust = TRUE).")
parser$add_argument("--trend", action = "store_true", default = FALSE,
                    help = "Use eBayes(trend = TRUE).")
## Other controls
parser$add_argument("--no-diagnostics", dest = "no_diagnostics", action = "store_true",
                    default = FALSE,
                    help = "Skip per-domain p-value / plotSA diagnostics.")
parser$add_argument("--plot-width", dest = "plot_width", type = "double", default = 11,
                    help = "Diagnostic PNG width in inches [11].")
parser$add_argument("--plot-height", dest = "plot_height", type = "double", default = 5,
                    help = "Diagnostic PNG height in inches [5].")
parser$add_argument("--seed", type = "integer", default = 20260101L)

args <- parser$parse_args()

## Process some of the args
set.seed(args$seed)


### Files that exist
infiles <- list(
    phenos = paste0("2_Pipeline/02-Prepare_analysis_data/phenos_cleaned_SVs_",
                    args$cohort, ".txt"),
    prot   = "2_Pipeline/01-Explore_data/prot.pcafiltered.log2norm.txt",
    e4_dat = "0_Data/From_v3/e4_phenos_allcohorts.txt"
)

## code to run joint model for conditional and joint F test (for basic run only)
script_run_basic_joint <- "1_Code/03-Run_regressions/01asub-run_basic_joint.r"


### Files to be created
outdir <- paste0("2_Pipeline/03-Run_regressions/Prot_perDomain/", Capwords(args$run), "/")
MyMkdir(outdir)

diag_dir <-  "Diagnostics/"
if (!args$no_diagnostics) {  MyMkdir(file.path(outdir, diag_dir)) }

fileid <- paste0(args$cohort, "_", args$run)
outfiles <- list(
    limma_per_domain = paste0(fileid, ".tsv"),
    summary_dt       = paste0(fileid, "_summary.tsv"),
    diagnostics      = paste0(diag_dir, fileid, "_diagnostics.tsv"),
    qq               = paste0(fileid, "_qqplots.pdf"),
    session_info     = "sessionInfo.txt"
)

if(args$run == "basic") {
    outfiles[["joint_per_domain"]] <- paste0(args$cohort, "_joint_perDomain.tsv")
    outfiles[["joint_ftest"]]     = paste0(args$cohort, "_joint_Ftest.tsv")
}

## These files are per domain, so define name using a function
outfiles[["SA_plot"]] <- function(adomain) {
    file.path(outdir, diag_dir, paste0(fileid, "_", make.names(adomain), "_pvalue_plotSA.png"))
}


### MAIN ###

### Set some values based on the run specifics
### ( More of this happens later also)
## npsvars for this run
if(args$run == "severity") {
    npsvars <- npsvars.sev
} else {
    npsvars <- npsvars.bin
}

## min non-missing donors to keep protein in analysis
mincount <- ifelse(args$run %in% c("males", "females"), 25, 50)


## covars for this run
## SVs will be added on the fly
if(args$run %in% c("males", "females")) {
    covars <- c("Batch", "pmi", "age_death")
    fac_covars <- c("Batch")
}

if(args$run %in% c("basic", "severity")) {
    covars <- c("Batch", "pmi", "age_death", "msex")
    fac_covars <- c("Batch")
}

if(args$run %in% c("e4_adj")) {
    covars <- c("Batch", "pmi", "age_death", "msex", "e4_expression")
    fac_covars <- c("Batch")
}


### Read data
pheno <- fread(infiles$pheno)
prot <- fread(infiles$prot)
e4_dat <- fread(infiles$e4_dat)

### Prepare data
sample_id_col <- "protsample"
protein_id_col <- names(prot)[1]

## Remove any existing SV columns
pheno.orig <- copy(pheno)
oldsvs <- grep("SV", names(pheno.orig), value = TRUE)
pheno[, (oldsvs) := NULL]


## Merge in e4 expression
pheno[e4_dat, e4_expression := i.e4_expression, on = "protsample"]


## Filter if necessary for this run
if(args$run == "males") {
    pheno <- pheno[msex == 1, ]
}

if(args$run == "females") {
    pheno <- pheno[msex == 1, ]
}

## Do some formatting
expr_all <- as.matrix(prot[, !protein_id_col, with = FALSE])
rownames(expr_all) <- as.character(prot[[protein_id_col]])
storage.mode(expr_all) <- "double"

pheno[, (sample_id_col) := as.character(get(sample_id_col))]
setkeyv(pheno, sample_id_col)


### Align samples between expression and phenos
shared <- intersect(colnames(expr_all), pheno[[sample_id_col]])
if (length(shared) == 0L) {
  stop("No sample IDs shared between --pheno and --prot.")
}
message("Samples shared between tables: ", length(shared))

expr_all <- expr_all[, shared, drop = FALSE]
pdat     <- pheno[J(shared)]                 # reordered to match columns
stopifnot(identical(as.character(pdat[[sample_id_col]]), colnames(expr_all)))

## Set covars as factors as necessary
for (nm in fac_covars) pdat[, (nm) := factor(get(nm))]


### limma drops samples with missing covariates from every protein, so
### remove them once, up front, and report the count
covar_ok <- pdat[, Reduce(`&`, lapply(.SD, function(v) !is.na(v))), .SDcols = covars]
nps_any  <- pdat[, Reduce(`|`, lapply(.SD, function(v) !is.na(v))), .SDcols = npsvars]
keep_samp <- covar_ok & nps_any
message("Samples dropped (missing covariate or no NPS observation): ", sum(!keep_samp))

expr_all <- expr_all[, keep_samp, drop = FALSE]
pdat     <- pdat[keep_samp]
for (nm in fac_covars) pdat[, (nm) := droplevels(get(nm))]


### Filter to proteins that pass missingness threshold
nonmiss_count <- rowSums(!is.na(expr_all))
keep_prot <- nonmiss_count >= mincount
message("Proteins retained at count >= ", mincount, " donors: ",
        sum(keep_prot), " of ", length(keep_prot))
expr <- expr_all[keep_prot, , drop = FALSE]

complete_rows <- rowSums(is.na(expr)) == 0L
message("Proteins with complete data (used for SVA): ", sum(complete_rows))


### Make list objects to hold results per domain
per_domain <- diag_list <- qq_list <- vector("list", length(npsvars))
names(per_domain) <- names(diag_list) <- names(qq_list) <- npsvars



### Run for each domain in turn
for (dom in npsvars) {
    ## Prepare data
    message("=== Per-domain model: ", dom)
    sub_pdat <- pdat[!is.na(get(dom))]   ## subset of donors with non-missing NPS for this domain
    samples  <- sub_pdat[[sample_id_col]]

    edat <- expr[, samples, drop = FALSE]

    ## Estimate SVs and define model
    sv     <- Estimate_svs(protected_vars = dom, pheno_dt = sub_pdat, sample_ids = samples)
    design <- Build_design(predictors = dom, covariates = covars,
                           pheno_dt = sub_pdat, sample_ids = samples, svs = sv)

    coef_name <- make.names(dom)
    if (!coef_name %in% colnames(design)) {
        warning("Coefficient ", coef_name, " not estimable; skipping ", dom)
        next
    }

    ## Fit model
    fit <- Fit_limma(edat, design)
    n_svs_used <- if (is.null(sv)) 0L else ncol(sv)

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
                    n_samples = length(samples),
                    n_svs = if (is.null(sv)) 0L else ncol(sv))
    setorder(res, P.Value)

    res[fdr < 0.05, sig.qvalue05 := "sig"]
    message("  proteins at FDR < 0.05: ", res[fdr < 0.05, .N])

    ## Save to cross-nps list
    per_domain[[dom]] <- res

    ## Run diagnostics and draw plot to file
    if (!args$no_diagnostics) {
        diag_list[[dom]] <- Pvalue_diagnostics(fit, res[["P.Value"]], dom,
                                               length(samples), n_svs_used)
        Write_diagnostic_plot(
            fit, res[["P.Value"]], dom, length(samples), n_svs_used,
            outfiles$SA_plot(dom)
        )
        message("  diagnostics: prop p < 0.05 = ",
                signif(diag_list[[dom]]$prop_p_lessthan_05, 3),
                ", pi0 = ", signif(diag_list[[dom]]$pi0, 3),
                ", df.prior = ", signif(diag_list[[dom]]$df_prior, 3))
    }


    ## calculate inflation factor
    chisq.all <- qchisq(res[, P.Value], 1, lower.tail=FALSE)
    lambda.median.all <- median(chisq.all) / qchisq(0.5,1)
    lambdatext.all <- bquote(lambda ~ '=' ~ .(signif(lambda.median.all, 4)))


    ##  Draw qq plot
    n.na <- sum(is.na(res$P.Value))
    titleqq <- paste0(dom, "\n", "N=", nrow(sub_pdat))
    if(n.na > 0) { titleqq <- paste0(titleqq, "; excluding ", n.na, " proteins with NA P-value") }
    qq_list[[dom]] <- qqunif.plot(res[!is.na(P.Value), P.Value],
                                  maintitle=titleqq, subxlab=lambdatext.all,
                                  par.settings=list(par.sub.text=list(cex=1)))
}


### Gather results across domains
per_domain_dt <- rbindlist(per_domain, use.names = TRUE, fill = TRUE)
setcolorder(per_domain_dt, c("cohort", "domain", "protein", "logFC", "CI.L", "CI.R",
                             "SE_unmod", "SE_mod",
                             "AveExpr", "t",
                             "P.Value", "fdr"))


### Gather diagnostics across domains
if (!args$no_diagnostics) {
    diag_dt <- rbindlist(diag_list, use.names = TRUE, fill = TRUE)
    if (nrow(diag_dt) > 0L) {
        setorder(diag_dt, domain)
        #print(diag_dt)

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

### Assemble summary
summary_dt <- per_domain_dt[, .(n_proteins = .N,
                                n_fdr05 = sum(fdr < 0.05, na.rm = TRUE),
                                n_samples = n_samples[1],
                                n_svs = n_svs[1]),
                            by = .(cohort, domain)]
print(summary_dt)



### Run joint model (basic run only)
if(args$run == "basic") {
    source(script_run_basic_joint)
}

message("Done.")


### FINISH ###
## Write results to files
fwrite(per_domain_dt, file.path(outdir, outfiles$limma_per_domain), sep = "\t")
fwrite(summary_dt, file.path(outdir, outfiles$summary_dt), sep = "\t")

if(args$run == "basic") {
    fwrite(joint_per_domain, file.path(outdir, outfiles$joint_per_domain), sep = "\t")
    fwrite(joint_F, file.path(outdir, outfiles$joint_ftest), sep = "\t")
}

if (!args$no_diagnostics) {
    fwrite(diag_dt, file.path(outdir, outfiles$diagnostics), sep = "\t")
}


pdf(file.path(outdir, outfiles$qq), height = 7, width = 7)
for (i in qq_list) { print(i) }
invisible(dev.off())
