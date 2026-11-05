library(data.table)
library(argparse)
library(ggplot2)
library(qvalue)
rm(list = ls())



###### GOAL ######



##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")

sampsizes <- c(OHSU = 194, rush = 104, emory = 73)

###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")
source("~/comsv/Svattathil_Library/manhattan_and_qq_plot_functions.r")


StoufferCombine <- function(pvals, N, weighting = "sqrtN") {
    ## Combine one-sided p-values (F-test p-values are upper-tail by construction)
    ## Z_i = qnorm(1 - p_i);  Z = sum(w_i * Z_i) / sqrt(sum(w_i^2))
    keep <- !is.na(pvals)
    pvals <- pvals[keep]; N <- N[keep]
    if(length(pvals) == 0) return(list(Z = NA_real_, Pvalue = NA_real_,
                                       n_cohorts = 0L, N_total = 0L))

    ## Keep p strictly inside (0, 1) so Z is finite
    pvals <- pmin(pmax(pvals, .Machine$double.xmin), 1 - .Machine$double.eps)
    z <- qnorm(pvals, lower.tail = FALSE)
    w <- if(weighting == "sqrtN") sqrt(N) else rep(1, length(z))
    Z <- sum(w * z) / sqrt(sum(w^2))

    list(Z = Z, Pvalue = pnorm(Z, lower.tail = FALSE),
         n_cohorts = length(pvals), N_total = sum(N))
}


###### SETUP ######
### Command line arguments
## Create parser
parser <- ArgumentParser()
#parser$add_argument("--", type="character", help="")
#parser$add_argument("--", action="store_true", default=TRUE, help="")
#parser$add_argument("--", action="store_true", default=FALSE, help="")
#parser$add_argument("--", action="store_false", dest="", help="")
#parser$add_argument("--", type="integer", default=10, help="")
#parser$add_argument("--", type="double", default=10, help="")

## Read from parser
args <- parser$parse_args()


### Files that exist
indir1 <- "2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/PermP_fromHang/"
indir2 <- "2_Pipeline/03-Run_regressions/Prot_perDomain/Basic_withSV/" ## to get protein-specific sample counts
infiles_pperm <- sapply(cohorts, function(x) {
    paste0(indir1, "NPS_7df_pperm_", Capwords(x), "_forQ.txt")
}, simplify = FALSE)

infiles_n <- sapply(cohorts, function(x) {
    paste0(indir2, "basic_withSV_", x, ".tsv")
}, simplify = FALSE)


### Files to be created
outdir <- "2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Postcompute/"
MyMkdir(outdir)
outfiles <- list(
    txt = paste0(outdir, "stouffer_combination.txt"),
    hist = paste0(outdir, "stouffer_combination_p_histogram.png"),
    qq = paste0(outdir, "stouffer_combination_qqplot.png")
)


###### MAIN ######
### Read in data
dat_list <- sapply(infiles_pperm, fread, simplify = FALSE)
n_list <- sapply(infiles_n, fread, simplify = FALSE)


### Prepare data
## Get protein-specific sample sizes
## take unique because there is one row per BPSD domain, but since we had no
## missing phenotype data the counts are identical for each domain
protein_n <- unique(rbindlist(n_list, idcol = "cohort")[, .(cohort, protein, n_obs)])


## Extract relevant columns to single cross-cohort tables, for simplicity
pperm <- rbindlist(dat_list, idcol = "cohort")[, .(cohort, protein, pperm)]
pperm[, cohort := factor(cohort, levels = cohorts)]


## Add sample sizes
pperm[protein_n, N := i.n_obs, on = c("cohort", "protein")]


stouffer_summ <- pperm[, {
    sc <- StoufferCombine(pperm, N)
    list(n_cohorts = sc$n_cohorts,
         N_total = sc$N_total,
         weighting = "sqrtN",
         Stouffer_Z = sc$Z,
         Pvalue     = sc$Pvalue)
}, by = "protein"]

setorder(stouffer_summ, Pvalue)


p1 <- qqunif.plot(stouffer_summ$Pvalue)
p2 <- ggplot(stouffer_summ, aes(x = Pvalue)) +
    geom_histogram(bins = 100, fill = "orange3", linewidth = 1) +
    theme_minimal()


###### FINISH ######
write.table(pperm, file = outfiles$txt, row.names = FALSE, quote = FALSE, sep = "\t")
ggsave(file = outfiles$hist, width = 5, height = 5, unit = "in")
png(outfiles$qq)
print(p1)
dev.off()
