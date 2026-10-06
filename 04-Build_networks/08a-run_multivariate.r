library(data.table)
library(argparse)
library(qvalue)
library(car)
library(geepack)
library(clubSandwich)

rm(list=ls())


### CONSTANTS ###
options(stringsAsFactors = FALSE)
options(nwarnings = 100000000)
source("1_Code/project_constants.r")

###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")
source("~/comsv/Svattathil_Library/manhattan_and_qq_plot_functions.r")

DomainCounts <- function(longtest) {
    ## Track per-domain case counts for any protein that errors, to help
    ## distinguish "aliased/constant domain" from other failure modes
    ## without having to re-run anything.
    tab <- longtest[, .(n = .N, cases = sum(outcome, na.rm = TRUE)), by = nps_id]
    paste(sprintf("%s:n=%d,cases=%d", tab$nps_id, tab$n, tab$cases), collapse = "; ")
}


###### SETUP ######
## Create parser
parser <- ArgumentParser()
parser$add_argument("--cohort", type="character", default = "OHSU")

## Read from parser
args <- parser$parse_args()


### Files that exist
indir <- "2_Pipeline/04-Build_networks/Wgcna_out/Consensus/Auto_softpower9/"
infiles <- list(phenos = paste0("2_Pipeline/02-Prepare_analysis_data/phenos_cleaned_SVs_",
                                args$cohort, ".txt"),
                eigen  = paste0(indir, "auto_softpower9_", args$cohort, "_eigengenes.txt")
                )


### Files to be created
outdir <- paste0(indir, "Multivariate/")
MyMkdir(outdir)

outdir.mvmeta <- paste0(outdir, "Mvmeta/In/")
MyMkdir(outdir.mvmeta)

fileid <- paste0("eigengene_multivar_", args$cohort)

outfiles <- list(
    summstats     = paste0(outdir, fileid, "_regression_stats.txt"),
    qq            = paste0(outdir, fileid, "_qq.png"),
    log           = paste0(outdir, fileid, ".log"),
    formeta_betas = paste0(outdir.mvmeta, fileid, "_betas.rds"),
    formeta_covs  = paste0(outdir.mvmeta, fileid, "_covariances.rds")
)


###### MAIN ######
#### Read in data ####
eigen.all <- fread(infiles$eigen)
phenos.all <- fread(infiles$phenos)


### Define model and model family
## For multivariate, we need to run for each NPS in turn
## so the modelform is a function
modelform <- function(currnps) { paste0(currnps, " ~ eigengene") }


#### Prepare data ####
### Filter samples
## Restrict phenos to samples in eigengene data (in case they did different filtering)
phenos.all <- phenos.all[protsample %in% eigen.all$protsample, ]


## ## Restrict phenos columns to phenotypes that will be used for analysis
representative_formula <- as.formula(modelform("agit"))
allvars <- setdiff(unique(c(npsvars.bin, all.vars(representative_formula))), "eigengene")
phenos <- phenos.all[, c("protsample", allvars), with = FALSE]


## Restrict eigen to samples in phenos and convert to matrix
setkey(eigen.all, protsample)
modcols <- setdiff(names(eigen.all), c("protsample", "MEgrey"))

eigen <- t(as.matrix(eigen.all[phenos.all$protsample, ..modcols]))
rownames(eigen) <- modcols
colnames(eigen) <- phenos.all$protsample


### Z-scale
print("Z-scaling data for each eigengene")
transformed <- t(apply(eigen, 1, scale))
rownames(transformed) <- rownames(eigen)
colnames(transformed) <- colnames(eigen)


### Iinitialize objects to hold results, and also coefficients and covariance matrices,
### which are necessary to run meta-analysis across cohorts
beta_vectors <- covariance_mats <- as.list(rep(NA, nrow(transformed)))   # initialize with NA
names(beta_vectors) <- names(covariance_mats) <- rownames(transformed)   # add eigengene names
results_list <- vector("list", nrow(transformed))


### Run regression for each eigengene in turn
print("Starting regressions for each eigengene")
## Set seed in case multivariate has some random function
set.seed(24983423)

for (i in seq_len(nrow(transformed))) {
    if (i %% 5 == 0) { print(paste0(pN(i), " of ", pN(nrow(transformed)))) }

    eigengene_name <- rownames(transformed)[i]
    domain_counts <- NA_character_   # populated inside tryCatch once longtest exists

    fit_result <- tryCatch({
        ### 1. Set up data
        ## Initiate test data for this eigengene (z-scaled residuals and sample ID)
        test <- data.frame(
            eigengene = unlist(transformed[i, ]),
            protsample = colnames(transformed))

        ## Add phenos to test data
        test <- merge(test, phenos, by="protsample")

        ## Restrict to complete data
        test <- test[complete.cases(test), ]
        setDT(test)

        ## Convert to long format
        idvars <- all.vars(as.formula(modelform("protsample"))) # cheaty way to get predictor
                                                                #  vars + sample id column
        longtest <- melt(test,
                         id.vars = idvars,
                         measure.vars = npsvars.bin,
                         variable.name = "nps_id",
                         value.name = "outcome")
        setDT(longtest)
        longtest[, protsample := factor(protsample)]
        longtest[, nps_id := factor(nps_id)]
        longtest[, wave := as.integer(nps_id)]
        setorder(longtest, protsample, wave)

        domain_counts <<- DomainCounts(longtest)

        ### 2. Run GEE model
        ## This formulation treats the NPS domains as repeated measures within each subject
        ## With interaction term, it estimates separate regression coefficients for each NPS
        ## Including nps_id term allows a different intercept per nps
        ## The working correlation captures the fact that the outcomes are correlated
        ## the default san.se method is downward-biased with small sample size and will
        ## cause inflation in the Wald p-value, but use it here for speed and we will
        ## replace that estimate with the one using CR2 bias-reduced sandwich covariates
        ##
        ## There will be 14 coefficients --
        ## intercept, 6 NPS main effects, 7 eigengene:NPS interaction effects
        gee_fit_int <- geeglm(update.formula(as.formula(modelform("outcome")),
                                             . ~ . -eigengene + nps_id + eigengene:nps_id),
                              id      = protsample,  # clustered by protsample
                              waves   = wave,        # used if "unstructured" corr structure
                              data    = longtest,
                              family  = binomial(link = "logit"),
                              corstr  = "unstructured",
                              std.err = "san.se")

        ## Detect aliasing explicitly rather than letting it surface later as an opaque
        ## singular-matrix error out of linearHypothesis
        if (anyNA(coef(gee_fit_int))) {
            stop("Aliased coefficient(s) in GEE fit (likely a constant/near-constant ",
                 "NPS domain for this eigengene): ",
                 paste(names(coef(gee_fit_int))[is.na(coef(gee_fit_int))], collapse = ", "))
        }

        ### 3. Extract coefficients for the interaction terms (eigengene:nps_id)
        coef_vec <- coef(gee_fit_int)
        prot_terms <- grep("eigengene", names(coef_vec), value = TRUE)
        beta_prot  <- coef_vec[prot_terms]   # length = n.nps
        names(beta_prot) <- sub(":eigengene", "", sub("nps_id", "", prot_terms))   # rename

        ### 4. Extract covariance matrix for *all* coefficients
        ## CR2 (bias-reduced) sandwich covariance — recomputed from residuals,
        ## independent of whatever std.err= was set in geeglm()
        V_full <- clubSandwich::vcovCR(
                                    gee_fit_int,     # square matrix (dim = #coeffs)
                                    cluster = longtest$protsample,
                                    type    = "CR2")

        ## Subset to the eigengene‑by‑NPS block
        ## V_prot is the sampling-error covariance matrix for the per-domain eigengene coefficients
        ## It accounts for the correlation among the outcomes through the GEE working correlation,
        ## and for any over-dispersion via the robust sandwich estimator
        V_prot <- V_full[prot_terms, prot_terms]   # n.nps x n.nps


        ### 5. Do joint hypothesis test with Wald test
        ## Small-sample-corrected joint test of the 7 eigengene:nps_id terms,
        ## using Hotelling-Thurston-Zhang (HTZ) F-test with adjusted df
        ## instead of a chi-square Wald test assuming K -> infinity
        wald_i <- clubSandwich::Wald_test(gee_fit_int,
                     constraints = constrain_zero(prot_terms),
                     vcov        = "CR2",
                     cluster     = longtest$protsample,
                     test        = "HTZ")


        ## Extract summary statistics
        n <- length(unique(gee_fit_int$id))
        stats_i <- data.frame(wald_i, nobs = n)

        list(beta = beta_prot, V = V_prot, stats = wald_i, error = NA_character_)

    }, error = function(e) {
        list(beta = NA, V = NA,
             stats = data.frame(test = NA_character_, Fstat = NA_real_, df_num = NA_integer_,
                       df_denom = NA_real_, pval = NA_real_, sig = NA_real_),
             error = conditionMessage(e))
    })

    beta_vectors[[i]]      <- fit_result$beta
    covariance_mats[[i]]   <- fit_result$V

    results_list[[i]] <- list(
        eigengene       = eigengene_name,
        Fstat         = unname(fit_result$stats$Fstat),
        P             = unname(fit_result$stats$p_val),
        nobs          = unname(fit_result$stats$nobs),
        error         = fit_result$error,
        domain_counts = domain_counts
    )
}

n.errors <- sum(!sapply(results_list, function(x) is.na(x$error)))
if (n.errors > 0) {
    print(paste0(pN(n.errors), " of ", pN(nrow(transformed)),
                 " eigengenes failed and were skipped (see 'error' column in output)."))
}


summstats <- rbindlist(lapply(results_list, as.data.table), fill = TRUE)
setcolorder(summstats, "eigengene")


### Add columns for multiple testing adjustment
summstats[, P.BH := p.adjust(P, method="BH")]
summstats[, P.bonf := p.adjust(P, method="bonferroni")]
summstats[, sig.nominal05 := ifelse(P < 0.05, "sig", NA_character_)] ## for easy grepping


### Sort by p-value
setorder(summstats, P)


### Add rank
summstats[, p.rank := rank(P, ties.method="max")]
setcolorder(summstats, "p.rank")


### Draw qq plot
n.na <- sum(is.na(summstats$P))
titleqq <- paste0(fileid, "\n", "N=", nrow(phenos))
if(n.na > 0) { titleqq <- paste0(titleqq, "; excluding ", n.na, " eigengenes with NA P-value") }
forqq <- summstats[!is.na(P), P] # get non-NA p-values
forqq[which(forqq == 0)] <- 1e-300  # make non-zero so qqplot works

qqplot <- qqunif.plot(forqq, maintitle=titleqq)


### Make log table 1 - basic count stats
logtab <- data.table(cohort            = args$cohort,
                     NPS               = args$nps,
                     model             = paste0(deparse(modelform), collapse = ""),
                     n.eigengenes.tested = nrow(transformed),
                     n.sig.nominal05    = sum(summstats$sig.nominal05 == "sig", na.rm=TRUE))



###### FINISH ######
### Write summary statistics
write.table(summstats, file = outfiles$summstats, row.names = FALSE, quote = FALSE, sep = "\t")


### Write log file
write("\nModel and significance summary:", file = outfiles$log)
fwrite(logtab, file = outfiles$log, sep = "\t", quote = FALSE, col.names = TRUE, append = TRUE)


### Write objects to do meta-analysis
saveRDS(beta_vectors, file = outfiles$formeta_betas)
saveRDS(covariance_mats, file = outfiles$formeta_covs)


### Draw qqplot
png(outfiles$qq, height=7, width=7, units="in", res=300)
print(qqplot)
dev.off()
