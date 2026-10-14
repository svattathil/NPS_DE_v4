library(data.table)
library(argparse)
library(qvalue)

rm(list = ls())


###### GOAL ######
## Run cross-cohort meta-analysis for specified outcome


##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")
source("~/comsv/Svattathil_Library/manhattan_and_qq_plot_functions.r")


###### SETUP ######
### Command line arguments
## Create parser
parser <- ArgumentParser()
parser$add_argument("--SE_type", type="character", default="SE_mod", help = "SE_mod or SE_unmod")


## Read from parser
args <- parser$parse_args()


### Files that exist
resdir <- paste0("2_Pipeline/03-Run_regressions/Prot_perDomain/Basic/")
Statsfile <- function(acohort) { paste0(resdir, acohort, "_joint_Ftest.tsv") }


### Files to be created
outdir <- paste0(resdir, "Cross-cohort_joint/")
MyMkdir(outdir)
outfiles<- list(stats = paste0(outdir, "ftest_cross-cohort_fisher_stouffer.txt"),
                qq    = paste0(outdir, "ftest_cross-cohort_fisher_stouffer_qqplots.pdf")
                )


###### MAIN ######
### Read in data
res.list.cohort <- sapply(cohorts, function(x) { fread(Statsfile(x)) }, simplify = FALSE)


### Prepare data
pvals.list <- sapply(cohorts, function(acohort) {
    curr_res <- res.list.cohort[[acohort]][, .(protein, P.Value, n_obs)]
    setnames(curr_res, c("P.Value", "n_obs"), c(acohort, paste0(acohort, "_n")))

    return(curr_res)
}, simplify = FALSE)

## This merge keeps only proteins present in all three cohorts
merged <- Reduce(function(x,y) { merge(x,y, by = "protein") }, pvals.list)


Combine_omnibus_p <- function(dt, cohorts = c("OHSU", "rush", "emory"),
                              n_suffix = "_n", min_cohorts = 3) {
    ## dt needs: protein, one p-value column per cohort, and one donor-count column per cohort
    ## named paste0(cohort, n_suffix)  (OHSU_n, RADC_n, emory_n)
    dt <- as.data.table(copy(dt))
    ncols <- paste0(cohorts, n_suffix)
    stopifnot(all(cohorts %in% names(dt)), all(ncols %in% names(dt)))

    pmat <- as.matrix(dt[, ..cohorts])
    nmat <- as.matrix(dt[, ..ncols])

    ## a cohort is usable only if it has a p-value AND a positive donor count
    bad <- !is.na(pmat) & (is.na(nmat) | nmat <= 0)
    if (any(bad)) {
        warning(sum(bad), paste0(" protein-cohort entries have a p-value",
                                 " but a missing/non-positive donor count; ",
                                 "those entries are excluded."))
        pmat[bad] <- NA
    }
    nmat[is.na(pmat)] <- NA
    pmat <- pmin(pmax(pmat, 1e-300), 1 - 1e-16)    # avoid exact 0 or 1; NA stays NA
    dim(pmat) <- dim(nmat)

    dt[, n_cohorts := rowSums(!is.na(pmat))]

    ## Fisher: unweighted; df = 2 x number of cohorts with usable data
    dt[, p_fisher := pchisq(-2 * rowSums(log(pmat), na.rm = TRUE),
                            df = 2 * n_cohorts, lower.tail = FALSE)]

    ## Stouffer: per-protein weights sqrt(n donors used for that protein in that cohort)
    z <- qnorm(pmat, lower.tail = FALSE)
    w <- sqrt(nmat); w[is.na(w)] <- 0
    dt[, p_stouffer := pnorm(-rowSums(z * w, na.rm = TRUE) / sqrt(rowSums(w^2)))]

    ## require a minimum number of contributing cohorts
    dt[n_cohorts < min_cohorts, c("p_fisher", "p_stouffer") := NA_real_]

    ## BH across proteins (NA p-values are excluded from the count)
    dt[, `:=`(BH_p_fisher = p.adjust(p_fisher, "BH"),
              BH_p_stouffer = p.adjust(p_stouffer, "BH"))]

    dt[, `:=`(qvalue_fisher = qvalue(p_fisher)$qvalue,
              qvalue_stouffer = qvalue(p_stouffer)$qvalue)]
    dt[]
}

res <- Combine_omnibus_p(merged)
setorder(res, p_stouffer)
res[qvalue_fisher < 0.05, sig_fisher := "sig"]
res[qvalue_stouffer < 0.05, sig_stouffer := "sig"]


###  Make qq plots
qq_list <- list()

for(amethod in c("fisher", "stouffer")) {
    pcol <- paste0("p_", amethod)
    pvals <- unlist(res[, ..pcol])

    ## calculate inflation factor
    chisq.all <- qchisq(pvals, 1, lower.tail=FALSE)
    lambda.median.all <- median(chisq.all) / qchisq(0.5,1)
    lambdatext.all <- bquote(lambda ~ '=' ~ .(signif(lambda.median.all, 4)))

    ## draw plot
    n.na <- sum(is.na(pvals))
    titleqq <- paste0(Capwords(amethod), " combined p-values")
    if(n.na > 0) { titleqq <- paste0(titleqq, "; excluding ", n.na, " proteins with NA P-value") }
    qq_list[[amethod]] <- qqunif.plot(pvals[!is.na(pvals)],
                                  maintitle=titleqq, subxlab=lambdatext.all,
                                  par.settings=list(par.sub.text=list(cex=1)))
}


###### FINISH ######
## Draw QQ plots to file
pdf(outfiles$qq)
par(mfrow = c(1, 2))
invisible(sapply(qq_list, print))
dev.off()

## Write results table to file
fwrite(res, file = outfiles$stats, quote = FALSE, sep = "\t", na = NA)

