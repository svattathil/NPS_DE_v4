library(data.table)
library(argparse)
library(readxl)

rm(list = ls())



###### GOAL ######
## Make table of enrichment statistics for single top enriched gene set per module.
##
## Currently focusing on GO BP sets, but code is
## flexible for also including GO MF or Reactome sets.

##### CONSTANTS #####
options(stringsAsFactors = FALSE)
cohorts <- c("OHSU", "rush", "emory")
npsvars.bin <- c(agitation = "agit", anxiety = "anx", apathy = "apa", delusion = "del",
                 depression = "depd", disinhibition = "disn", hallucination = "hall",
                 irritability = "irr", sleep = "nite")


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### COMMAND LINE ARGUMENTS ######
## Create parser
parser <- ArgumentParser()
parser$add_argument("--ntop", type="integer", default=1,
                    help="number of gene sets to keep per source")

parser$add_argument("--thresh", type="double", default=-1, help="adjusted p-value threshold")

## Read from parser
args <- parser$parse_args()


###### SETUP ######
### Files that exist
dir_enrichment <- paste0("2_Pipeline/04-Build_networks/Consensus/",
                         "Auto_softpower8/Enrichment/")
infiles <- list(enrichment = paste0(dir_enrichment,
                                    "consensus_modules_auto_softpower8_",
                                    "pthresh1_gene_set_enrichment.xlsx"))


### Files to be created
## if args$thresh was provided at command line, use that to select gene sets
## Otherwise, use number of top gene sets from args$ntop
if(args$thresh > 0) {
outfiles <- list(summtab = paste0(dir_enrichment,
                                      "consensus_modules_auto_softpower8_",
                                      "adjp", args$thresh, "_enriched_processes.txt"))
} else {
    outfiles <- list(summtab = paste0(dir_enrichment,
                                      "consensus_modules_auto_softpower8_",
                                      "top", args$ntop, "_enriched_processes.txt"))
}


###### MAIN ######
### Read in data
wb_path <- infiles[["enrichment"]]
sheets  <- excel_sheets(wb_path)

res_list <- lapply(sheets[-length(sheets)], function(sht) {
    ## last sheet is summary of gene counts, don't read that into this object
    ## read_excel returns a tibble → convert to data.table instantly
    setDT(read_excel(wb_path, sheet = sht))
})

## Name the list elements for easy reference
names(res_list) <- sheets[-length(sheets)]

## Flatten to single table with column for module name
## Use fill=TRUE in case of empty tables
res_enrich <- rbindlist(res_list, idcol = "module_color", fill = TRUE)

## Read in gene counts
res_genecounts <- setDT(read_excel(wb_path, sheet = tail(sheets, 1)))


### Prepare data
### Tables to merge should have
### column module_color
### And any other columns

### Prepare gene counts
setnames(res_genecounts, "set", "module_color")


### Prepare enrichment resuts
## Define which databases to include sets from
dbs_to_keep <- c("GO_BP")  ## c("GO_BP", "GO_MF", "Reactome")

## Grab the desired top sets per database
top_sets_allcol <- rbindlist(lapply(dbs_to_keep, function(asource) {
    ## Extract results for this source
    this_source <- res_enrich[Source == asource, ]

    ## Extract top enriched sets based on p-values
    setorder(this_source, pvalue)
    this_source[, p.rank := 1:nrow(.SD), by = c("module_color")]

    if(args$thresh > 0) {
        top <- this_source[p.adjust < args$thresh, .SD, by = c("module_color")]
    } else {
        top <- this_source[, head(.SD, args$ntop), by = c("module_color")]
    }
    return(top)
}))

## Make version with selected few columns
keepcols <- c("ID", "Description", "pvalue", "p_BH",
              "GeneRatio", "BgRatio", "RichFactor", "FoldEnrichment",
              "geneID", "Count")
setnames(top_sets_allcol, c("p.adjust"), c("p_BH"))
top_sets <- copy(top_sets_allcol[, c("module_color", "p.rank", "Source", keepcols),
                                 with = FALSE])

## Convert to wide format
enrich_for_merge <- dcast(
    top_sets, module_color + p.rank ~ Source,
    value.var = keepcols,
    sep = "_"
    )


### Merge in gene counts
### Will be repeated across rows per module; that's fine
forprint <- merge(res_genecounts, enrich_for_merge,
                  all = TRUE, by = "module_color")


### Reorder columns
enrichcols_sorted <- paste(
    rep(keepcols, times = length(dbs_to_keep)),
    rep(dbs_to_keep, each = length(keepcols)),
    sep = "_")

setcolorder(forprint, c("module_color", "n_symbol", "n_entrez", "p.rank",
                        enrichcols_sorted))


### Update column names
setnames(forprint, names(forprint),
         gsub("pvalue", "p-value", names(forprint)))
setnames(forprint, names(forprint),
         gsub("p_BH", "BH adjusted p-value", names(forprint)))


### Update row order
## order by set size, like WGCNA does
setorder(forprint, -n_symbol, na.last = TRUE)

### Remove grey 'module'
forprint <- forprint[module_color != "grey", ]


###### FINISH ######
write.table(forprint, file = outfiles[["summtab"]],
            row.names = FALSE, quote = FALSE, sep = "\t")
