library(data.table)
library(argparse)
library(readxl)
library(stringr)

rm(list = ls())



###### GOAL ######
# Make table summarizing module eigengene associations to use as main table in manuscript


##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")


nps_map <- data.table(nps = npsvars.bin, nps_text = names(npsvars.bin))
setkey(nps_map, nps)

###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")

MergeRes <- function(alist) { Reduce(function(x, y) {
    merge(x, y, by = "module", all = TRUE) }, alist) }

###### SETUP ######


### Files that exist
indir <- "2_Pipeline/04-Build_networks/Wgcna_out/Consensus/Auto_softpower9/"
inpref <- "auto_softpower9"

infiles <- list()
infiles[["module_counts"]] <- paste0(indir, inpref, "_module_protein_counts.txt")
infiles[["kme"]] <- paste0(indir, inpref, "_mean_kME.txt")
infiles[["sig_nps"]] <- paste0(indir, "nps_regression_sig_", inpref, ".txt")
infiles[["enrichment"]] <- paste0(indir, "Enrichment/consensus_modules_",
                                  inpref, "_pthresh0.1_gene_set_enrichment.xlsx")
infiles[["sig_latfac"]] <- sapply(cohorts, function(acohort) {
    paste0(indir, "DEA_latent_factors/", inpref, "_", acohort,
           "_associated_modules_latent_factors.txt") })

infiles[["manual_order"]] <- "0_Data/module_manual_order.txt"
infiles[["multi"]] <- sapply(cohorts, function(acohort) {
    paste0(indir, "Multivariate/",
           "eigengene_multivar", "_", acohort, "_regression_stats.txt") })

infiles[["mvmeta"]] <- paste0(indir, "Multivariate/Mvmeta/Out/wald_results.txt")


### Files to be created
outdir <- "2_Pipeline/Manuscript/"
MyMkdir(outdir)
outfiles <- list(txt = paste0(outdir, "table_consensus_module_association.txt"))


###### MAIN ######
### Read in data
module_counts <- fread(infiles[["module_counts"]])
kme <- fread(infiles[["kme"]])
sig_nps <- fread(infiles[["sig_nps"]])

Readsheets <- function() {
    sheets  <- excel_sheets(infiles[["enrichment"]])
    #print(sheets)   # just to verify the names

    res_list <- lapply(sheets[-length(sheets)], function(sht) {
        ## last sheet is summary of gene counts, don't read that into this object
        ## read_excel returns a tibble → convert to data.table instantly
        setDT(read_excel(infiles[["enrichment"]], sheet = sht))
    })

    ## Name the list elements for easy reference
    names(res_list) <- sheets[-length(sheets)]

    ## Flatten to single table with column for module name
    ## Use fill=TRUE in case of empty tables
    restab <- rbindlist(res_list, idcol = "module_color", fill = TRUE)
    return(restab)
}
res_enrichment <- Readsheets()

sig_latfac_list <- sapply(infiles[["sig_latfac"]], fread, simplify = FALSE)
sig_multi_list <- sapply(infiles[["multi"]], fread, simplify = FALSE)
mvmeta_allres <- fread(infiles[["mvmeta"]])
manual_order <- fread(infiles[["manual_order"]])


### Prepare data
### Initiate list to hold objects to bind for overview table
### Should have one row per module
forbind <- vector(mode = "list")


### Prepare hub gene counts
forbind[["module_counts"]] <- module_counts[, .(module, hub_genes)]


### Prepare top enriched BP
setnames(res_enrichment, "module_color", "module")
res_enrichment[, module := factor(module, levels = module_counts$module)]
res_bp <- res_enrichment[Source == "GO_BP", -c("geneID")]
setorder(res_bp, pvalue)
forbind[["enrichment"]] <- res_bp[, head(.SD, 1), by = "module",
                                  .SDcols = c("ID", "Description")]


### Prepare top hub genes
kme[, gene := sapply(protein, function(x) {
    unlist(strsplit(x, split = "\\|"))[1] })]

setorder(kme, -kME_avg)
forbind[["top_hub_genes"]] <- kme[kME_avg > 0.7,
                                  list(top_hub_genes = paste0(head(gene, 10),
                                                              collapse = ", ")),
                                  by = "module"]


### Prepare associations
### First get per-cohort NPS associations
setnames(sig_nps, "NPS", "nps")
sig_nps <- merge(sig_nps, nps_map, by = "nps")

sigstrings_by_cohort <- list()
sigstrings_by_cohort[["nps"]] <-
    sapply(cohorts, function(acohort) {
        tmp <- sig_nps[cohort == acohort,
                       paste0(nps_text, collapse = ", "),
                       by = c("module")]
        setnames(tmp, "V1", paste0("sig_nps_", acohort))
        return(tmp)
    }, simplify = FALSE)


### Then get per-cohort latent factor associations
## rename f1 to latent factor 1
sigstrings_by_cohort[["factor"]] <- sapply(cohorts, function(acohort) {
    sig_latfac_list[[acohort]][, latentfactor_sig05 := gsub("f",
                                                            "latent factor ",
                                                            latentfactor_sig05)]
    sig_latfac_list[[acohort]][, cohort := acohort]

    sigres <- sig_latfac_list[[acohort]][grep("actor", latentfactor_sig05),
                                         paste0(latentfactor_sig05,
                                                collapse = ", "),
                                         by = c("module")]
    setnames(sigres, "V1", paste0("sig_factor_", acohort))
    return(sigres)
}, simplify = FALSE)


### Then get per-cohort multivariate associations
sigstrings_by_cohort[["multi"]]<- sapply(cohorts, function(acohort) {
    sig_multi_list[[acohort]][, module := gsub("ME", "", eigengene)]
    sigres <- sig_multi_list[[acohort]][P < .05, .(module,
                                                   sig_multi = "multivariate")]
    setnames(sigres, "sig_multi", paste0("sig_multi_", acohort))
    return(sigres)
}, simplify = FALSE)


### Flatten per-cohort results
## to one table per outcome
sigstrings_by_outcome <- sapply(sigstrings_by_cohort, MergeRes, simplify = FALSE)

## to one table with a column per cohort*outcome
sigstrings_assoc <- MergeRes(sigstrings_by_outcome)

## Make columns with all sig outcomes per cohort
col_prefixes <- c("sig_nps_", "sig_factor_", "sig_multi_")

sigstrings_assoc[, assoc_OHSU := apply(.SD, 1, function(x) {
    y <- x[!is.na(x)]; paste(y, collapse = ", ") } ),
    .SDcol = paste0(col_prefixes, "OHSU")]
sigstrings_assoc[, assoc_rush := apply(.SD, 1, function(x) {
    y <- x[!is.na(x)]; paste(y, collapse = ", ") } ),
    .SDcol = paste0(col_prefixes, "rush")]
sigstrings_assoc[, assoc_emory := apply(.SD, 1, function(x) {
    y <- x[!is.na(x)]; paste(y, collapse = ", ") } ),
    .SDcol = paste0(col_prefixes, "emory")]

forbind[["assoc"]] <- sigstrings_assoc[, .(module, assoc_OHSU,
                                           assoc_rush, assoc_emory)]

## Merge into single table
forprint <- MergeRes(forbind)


### Assemble meta-analysis associations
## from per-domain tests
sig_meta_perdomain <- sig_nps[fe_p < 0.05, paste0(nps_text, collapse = ", "),
                             by = "module"]
setnames(sig_meta_perdomain, "V1", "domains")


## from mvmeta on multivariate results
mvmeta_allres[, module := gsub("ME", "", eigengene)]
setkey(mvmeta_allres, module)
sig_mvmeta <- mvmeta_allres[sig.nominal05 == "sig", ]
sig_mvmeta[, multivariate := "multivariate"]
forbind[["meta"]] <- merge(sig_meta_perdomain, sig_mvmeta[, .(module, multivariate)],
                           by = "module", all = TRUE)

## make combined string of outcomes
CombineStrings <- function(x, y) {
    if(is.na(x)) { return(y) }
    if(is.na(y)) { return(x) }
    paste0(x, ", ", y)
}

forbind[["meta"]][, assoc_meta := mapply(CombineStrings, domains, multivariate)]
forbind[["meta"]][, c("domains", "multivariate") := NULL]

### Combine meta-analysis and per-cohort results
forprint <- merge(forprint, forbind[["meta"]], by = "module", all = TRUE)


### Make it pretty
## order using my manual order
forprint <- merge(manual_order, forprint, by = "module", all = TRUE)
setorder(forprint, manual_order, na.last = TRUE)

## Edit case
## Use sentence case for description
## Then fix because sentence case messes up 'mRNA'
forprint[, Description := str_to_sentence(Description)]
forprint[, Description := gsub("Mrna", "mRNA", Description)]

setcolorder(forprint, c("manual_order", "ID", "module",
                        "Description", "hub_genes", "top_hub_genes"))


###### FINISH ######
fwrite(forprint, file = outfiles[["txt"]],
       row.names = FALSE, quote = FALSE, sep = "\t")
