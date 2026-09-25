library(data.table)
library(argparse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(openxlsx)

rm(list = ls())

###### GOAL ######
# For each module, run GO BP, GO MF, and Reactome enrichment using gene symbols.

###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######
### Command line arguments
parser <- ArgumentParser()
parser$add_argument("--network_set", type = "character",
                    help = "Consensus or Cohort-specific",
                    default = "Consensus")
parser$add_argument("--network_group", type = "character",
                    help = "e.g., Auto_softpower7",
                    default = "Auto_softpower9")
parser$add_argument("--runid", type = "character",
                    help = "",
                    default = "auto_softpower9_OHSU")
parser$add_argument("--pthresh", type = "double",
                    help = "threshold for padjust for enrichment",
                    default = 0.10)
parser$add_argument("--truncate", action="store_true", default=FALSE, help="")


args <- parser$parse_args()


### Constants
options(stringsAsFactors = FALSE)


### Files that exist
if(args$network_set == "Cohort-specific") {
    indir <- "2_Pipeline/04-Build_networks/Cohort-specific/Selected_results/"

}

if(args$network_set == "Consensus") {
    indir <- paste0("2_Pipeline/04-Build_networks/Wgcna_out/Consensus/", args$network_group, "/")
}

files <- list(gene_stats = paste0(indir, args$runid, "_gene_conn_stats.txt")
              )


### Files to be created
out_dir <- paste0(indir, "Enrichment/")
MyMkdir(out_dir)

if(args$network_set == "Cohort-specific") {
    excel_outpref <- file.path(out_dir, paste0(args$runid, "_GO_Reactome_enrichment.xlsx"))
}
if(args$network_set == "Consensus") {
    excel_outpref <- file.path(out_dir, paste0("consensus_modules_",
                                               tolower(args$network_group), "_pthresh",
                                               args$pthresh,
                                               "_gene_set_enrichment.xlsx"))
}

if(args$truncate) {
    excel_outpref <- paste0(excel_outpref, "_truncated")
}

excel_outfile <- paste0(excel_outpref, ".xlsx")

wb <- createWorkbook()


###### MAIN ######
### Read data
gene_stats <- fread(files[["gene_stats"]])


### Prepare data
### Add gene column
gene_stats[, symbol := sapply(protein, function(x) {
    unlist(strsplit(x, split = "\\|"))[1] })]

### Extract background genes
bg_symbols <- gene_stats$symbol

## Remove NAs (proteins that were missing gene symbols)
## And convert to Entrez ID
bg_symbols <- bg_symbols[!is.na(bg_symbols)]
bg_ids_df <- clusterProfiler::bitr(bg_symbols, fromType = "SYMBOL",
                                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
bg_entrez <- unique(bg_ids_df$ENTREZID)


### Define modules to test for enrichment
testmodules <- sort(setdiff(unique(gene_stats$module), "grey"))


### Initiate table for module gene counts
gene_counts <- data.frame(set = "background",
                          n_symbol = length(bg_symbols),
                          n_entrez = length(bg_entrez))


### Loop over each module gene set
for (amodule in testmodules) {
    cat("--- Processing module", amodule, "---\n")

    ## Extract gene list
    ## And convert to Entrez ID
    symbols <- unlist(gene_stats[module == amodule, symbol])
    ids_df <- bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    entrez <- unique(ids_df$ENTREZID)

    ## Add gene count to summary table
    gene_counts <- rbind(gene_counts, data.table(
                                          set = amodule,
                                          n_symbol = length(symbols),
                                          n_entrez = length(entrez)))


    ### Store results in a list of data.tables
    results_list <- list()


    ## ---- GO Biological Process ----
    cat("...Running GO BP enrichment\n")
    ego_bp <- enrichGO(
        gene = entrez,
        universe = bg_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = args$pthresh,
        qvalueCutoff = 1.0
    )

    results_list[["GO_BP"]] <- as.data.table(ego_bp)

    ## ---- GO Molecular Function ----
    cat("...Running GO MF enrichment\n")
    ego_mf <- enrichGO(
        gene = entrez,
        universe = bg_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "MF",
        pAdjustMethod = "BH",
        pvalueCutoff = args$pthresh,
        qvalueCutoff = 1.0
    )
    results_list[["GO_MF"]] <- as.data.table(ego_mf)

    ## ---- Reactome Pathways ----
    cat("...Running Reactome enrichment\n")
    react <- enrichPathway(
        gene = entrez,
        universe = bg_entrez,
        organism = "human",
        pAdjustMethod = "BH",
        pvalueCutoff = args$pthresh,
        qvalueCutoff = 1.0,
        readable = TRUE
    )
    results_list[["Reactome"]] <- as.data.table(react)


    ### Create one table per module by combining results across sources
    all_results <- rbindlist(lapply(names(results_list), function(name) {
        dt <- results_list[[name]]
        if (nrow(dt) > 0) { dt[, Source := name] }


        ## if more than 100 pathways, then just output the top 100
        ## to keep the output file from being too huge
        setorder(dt, pvalue)
        if(args$truncate) {
            if( nrow(dt) > 100) { dt <- head(dt, 100) }
        }
        return(dt)
    }), use.names = TRUE, fill = TRUE)

    ## order by source and pvalue
    setorder(all_results, Source, pvalue)

    ### Write table to Excel sheet
    addWorksheet(wb, amodule)
    writeDataTable(wb, sheet = amodule, x = all_results)


    ### Clean up before next module
    rm(list = c("all_results", "results_list",
                "ego_bp", "ego_mf", "react", "entrez"))
}
rm(amodule)


addWorksheet(wb, "gene_counts")
writeDataTable(wb, sheet = "gene_counts", x = gene_counts)


### FINISH ###
saveWorkbook(wb, excel_outfile, overwrite = TRUE)


