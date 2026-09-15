library(data.table)
library(argparse)
library(qvalue)
library(car)
library(geepack)

rm(list=ls())


### CONSTANTS ###
options(stringsAsFactors = FALSE)
options(nwarnings = 100000000)
source("1_Code/project_constants.r")

###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")
source("~/comsv/Svattathil_Library/manhattan_and_qq_plot_functions.r")


GetStats <- function(regression.out, modelfamily = modelfamily, var="protein") {
    if(length(regression.out) == 1) {
        ## this can happen when robust regression fails and tryCatch returns NA
        return(c(Estimate = NA, SE = NA, P = NA, nobs = regression.out))
    }

    stats <- c(Chisq = regression.out[2, "Chisq"],
               P     = regression.out[2, "Pr(>Chisq)"])

    ## Add number of samples that were included in analysis
    nobs = regression.out[1, "Res.Df"] / length(npsvars.bin)
    toreturn <- c(stats, nobs = nobs)

    return(toreturn)
}


###### SETUP ######
## Create parser
parser <- ArgumentParser()
parser$add_argument("--cohort", type="character", default = "OHSU")

## Read from parser
args <- parser$parse_args()


### Files that exist
indir <- "2_Pipeline/02-Prepare_analysis_data/"
infiles <- list(phenos = paste0(indir, "phenos_cleaned_SVs_", args$cohort, ".txt"),
                resid  = paste0(indir, "resid_regress_covars_", args$cohort, ".txt")
                )

### Files to be created
outdir <- paste0("2_Pipeline/03-Run_regressions/Prot_multivar/")
MyMkdir(outdir)

outdir.mvmeta <- "Mvmeta/In/"
MyMkdir(outdir.mvmeta)

fileid <- paste0("prot_multivar_", args$cohort)

outfiles <- list(
    summstats     = paste0(outdir, fileid, "_regression_stats.txt"),
    qq            = paste0(outdir, fileid, "_qq.png"),
    log           = paste0(outdir, fileid, ".log"),
    formeta_betas = paste0(outdir, outdir.mvmeta, fileid, "_betas.rds"),
    formeta_covs  = paste0(outdir, outdir.mvmeta, fileid, "_covariances.rds")
)


###### MAIN ######
#### Read in data ####
resid.all <- fread(infiles$resid)
phenos.all <- fread(infiles$phenos)


### Define model and model family
## For multivariate, we need to run for each NPS in turn
## so the modelform is a function
modelfamily <- "binomial"
predictorvar <- "protein"
modelform <- function(currnps) { paste0(currnps, " ~ protein") }

#### Prepare data ####
### Filter samples
## Restrict phenos to samples in resid data (because they did different filtering)
phenos.all <- phenos.all[protsample %in% colnames(resid.all), ]


## ## Restrict phenos columns to phenotypes that will be used for analysis
representative_formula <- as.formula(modelform("agit"))
allvars <- setdiff(unique(c(npsvars.bin, all.vars(representative_formula))), "protein")
phenos <- phenos.all[, c("protsample", allvars), with = FALSE]


## Restrict resid to samples in phenos and convert to matrix
resid <- as.matrix(resid.all[, phenos.all$protsample, with = FALSE])
rownames(resid) <- resid.all$protein


### Select proteins with data for at least 50 participants
tokeep <- apply(resid, 1, function(x) {sum(!is.na(x)) > 50 })
for_transformation <- resid[tokeep, ]
rm(tokeep)


### Z-scale
print("Z-scaling data for each protein")
transformed <- t(apply(for_transformation, 1, scale))
rownames(transformed) <- rownames(for_transformation)
colnames(transformed) <- colnames(for_transformation)


### Iinitialize objects to hold coefficients and covariance matrix,
### which are necessary to run meta-analysis across cohorts
beta_vectors <- covariance_mats <- as.list(rep(NA, nrow(transformed)))   # initialize with NA
names(beta_vectors) <- names(covariance_mats) <- rownames(transformed)   # add protein names


### Run regression for each protein in turn
print("Starting regressions for each protein")
## Set seed since multivariate seems to have some random function
set.seed(24983423)
summstats <- data.table(data.frame(t(sapply(1:nrow(transformed), function(i) {
    if(i %% 1000 == 0) { print(paste0(pN(i), " of ", pN(nrow(transformed)))) }

    ## 1. Set up data
    ## Initiate test data for this protein (z-scaled residuals and sample ID)
    test <- data.frame(
        protein = unlist(transformed[i, ]),
        protsample = colnames(transformed))

    ## Add phenos to test data
    test <- merge(test, phenos, by="protsample")


    ## Restrict to complete data
    test <- test[complete.cases(test), ]
    setDT(test)

    ## Convert to long format
    ## This is a cheaty way to get predictor variables + sample id column
    idvars <- all.vars(as.formula(modelform("protsample")))
    longtest <- melt(test,
                     id.vars = idvars,
                     measure.vars = npsvars.bin,
                     variable.name = "nps_id",
                     value.name = "outcome")
    setDT(longtest)
    longtest[, protsample := factor(protsample)]
    longtest[, nps_id := factor(nps_id)]

    ## 2. Run GEE model
    ## This formulation treats the nine NPS as repeated measures within each subject
    ## With interaction term, it estimates separate regression coefficients for each NPS
    ## The working correlation captures the fact that the outcomes are correlated
    gee_fit_int <- geeglm(update.formula(as.formula(modelform("outcome")),
                                         . ~ . -protein + protein:nps_id),
                          id     = protsample,  # clustered by protsample
                          data   = longtest,
                          family = binomial(link = "logit"),
                          corstr = "exchangeable")

    ## 3. Extract coefficients for the interaction terms (protein:nps_id)
    coef_vec <- coef(gee_fit_int)
    prot_terms <- grep("^protein:", names(coef_vec), value = TRUE)
    beta_prot  <- coef_vec[prot_terms]   # length = 9
    names(beta_prot) <- sub("^protein:nps_id", "", prot_terms)   # rename

    ## 4. Extract robust (sandwich) covariance matrix for *all* coefficients
    V_full <- vcov(gee_fit_int)   # square matrix (dim = #coeffs)

    ## Subset to the protein‑by‑NPS block
    ## V_prot is the sampling-error covariance matrix for the nine protein coefficients
    ## It accounts for the correlation among the outcomes through the GEE working correlation,
    ## and for any over-dispersion via the robust sandwich estimator
    V_prot <- V_full[prot_terms, prot_terms]   # 9 × 9

    ## 5. Do joint hypothesis test with Wald test
    ## Wald test (Chi‑square approximation)
    ## Save as lm_i for compatibility
    lm_i <- car::linearHypothesis(gee_fit_int,
                                  hypothesis.matrix = prot_terms)

    beta_vectors[[i]] <<- beta_prot
    covariance_mats[[i]] <<- V_prot

    ## Extract and return summary statistics
    stats_i <- GetStats(lm_i, modelfamily = modelfamily, var = predictorvar)
    return(stats_i)
}))))


## Add protein column
summstats[, protein := rownames(transformed)]
setcolorder(summstats, "protein")


### Add columns for qvalue
summstats[, qvalue := qvalue(P)$qvalues]
summstats[, lfdr := qvalue(P)$lfdr]
summstats[, P.BH := p.adjust(P, method="BH")]
summstats[, P.bonf := p.adjust(P, method="bonferroni")]
summstats[, sig.qvalue05 := ifelse(qvalue < 0.05, TRUE, FALSE)]
summstats[, sig := ifelse(qvalue < 0.05, "sig", NA_character_)] ## for easy grepping


### Sort by p-value
setorder(summstats, P)


### Add rank
summstats[, p.rank := rank(P, ties.method="max")]
setcolorder(summstats, "p.rank")


### Draw qq plot
n.na <- sum(is.na(summstats$P))
titleqq <- paste0(fileid, "\n", "N=", nrow(phenos))
if(n.na > 0) { titleqq <- paste0(titleqq, "; excluding ", n.na, " proteins with NA P-value") }
qqplot <- qqunif.plot(summstats[!is.na(P), P], maintitle=titleqq)


### Make log table 1 - basic count stats
logtab <- data.table(cohort            = args$cohort,
                     NPS               = args$nps,
                     model             = paste0(deparse(modelform), collapse = ""),
                     n.proteins.tested = nrow(transformed),
                     n.sig.qvalue05    = sum(summstats$sig.qvalue05, na.rm=TRUE))



###### FINISH ######
### Write summary statistics
write.table(summstats, file = outfiles$summstats, row.names = FALSE, quote = FALSE, sep = "\t")


### Write log file
write("\nModel and significance summary:", file = outfiles$log)
fwrite(logtab, file = outfiles$log, sep = "\t", quote = FALSE, col.names = TRUE, append = TRUE)


### Write objects to do  meta-analysis
saveRDS(beta_vectors, file = outfiles$formeta_betas)
saveRDS(covariance_mats, file = outfiles$formeta_covs)


### Draw qqplot
png(outfiles$qq, height=7, width=7, units="in", res=300)
print(qqplot)
dev.off()
