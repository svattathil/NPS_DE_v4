library(data.table)
library(argparse)
library(WGCNA)
library(dplyr)
library(stringr)

rm(list = ls())



###### GOAL ######
## Build cross-cohort consensus co-expression network


##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")
softpowers_using_shared_genes <- c(OHSU = 9, rush = 9, emory = 9)

nps <- npsvars.bin


##### COMMAND LINE ARGUMENTS #####
## Create parser
parser <- ArgumentParser()
parser$add_argument("--auto_softpower", type="integer", default=9, help="")
parser$add_argument("--use_powervec1", action="store_true", default=FALSE, help="")
parser$add_argument("--residset", type="character", default="ignoreSVs", help="")

## Read from parser
args <- parser$parse_args()

## Define some values based on command line arguments
if(args$use_powervec1) {
    network_group  <- "powervec1"
    softpower <- softpowers_using_shared_genes
} else{
    network_group <- paste0("auto_softpower", args$auto_softpwer)
    softpower <- args$auto_softpower
}


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")

GetCoefficientAndSE <- function(regression.out, type, var="miRNA") {
    if(type=="ordinal") {
        return(summary(regression.out)$coefficients[var, c("Value", "Std. Error")]) }
    else if(type %in% c("gaussian", "binomial", "lmm")) {
        return(summary(regression.out)$coefficients[var, c("Estimate", "Std. Error")])
    }
}

GetPvalue <- function(regression.out, type, var="miRNA") {
    if(type=="ordinal") {
        return(pnorm(abs(summary(regression.out)$coefficient[var, "t value"]),
                     lower.tail=FALSE) * 2) }
    else if(type %in% c("binomial", "neg.binomial")) {
            summary(regression.out)$coefficients[var, "Pr(>|z|)"] }
    else if(type %in% c("gaussian", "lmm")) {
                summary(regression.out)$coefficients[var, "Pr(>|t|)"] }
}

PrepareDataME <- function(modulename, aphenodata, amedata) {
    ## Extract eigengene for current module and add protsample as column
    mes.to.merge <- data.frame(protsample=rownames(amedata),
                               mod.eigen=amedata[,modulename])

    ## Add eigengene column to phenos
    merged <- merge(mes.to.merge, aphenodata, by = "protsample")

    ## Convert to data frame and add protsample as rownames
    samps <- merged$protsample
    setDF(merged)
    rownames(merged) <- samps
    return(merged)
}


###### SETUP ######
### Files that exist
infiles <- list()
indir <- "2_Pipeline/04-Build_networks/Wgcna_out/Using_shared_genes/"
infiles[["shared_genes"]] <- "2_Pipeline/04-Build_networks/shared_genes_for_consensus_module.txt"

infiles[["toms"]] <- sapply(cohorts, function(acohort) {
    paste0(indir,
           "for_consensus_", acohort, "_", args$residset, "_softpower",
           softpowers_using_shared_genes[acohort],
           "_cutheight0.30_", "tom-block.1.RData") })


infiles[["expr"]] <- sapply(cohorts, function(acohort) {
    paste0(indir,
           "for_consensus_", acohort, "_", args$residset, "_softpower",
           softpowers_using_shared_genes[acohort],
           "_cutheight0.30_", "expr_data.txt") })

infiles[["phenos"]] <- sapply(cohorts, function(acohort) {
    paste0("2_Pipeline/02-Prepare_analysis_data/phenos_cleaned_SVs_", acohort, ".txt") })


### Files to be created
outdir   <- paste0("2_Pipeline/04-Build_networks/Wgcna_out/Consensus/", stringr::str_to_title(network_group), "/")
MyMkdir(outdir)

outfiles <- list()
Outfilepref <- function() { paste0(outdir, network_group, "_", acohort) }
outfiles[["Gene_conn_stats"]] <- function(acohort) { paste0(Outfilepref(), "_gene_conn_stats.txt") }
outfiles[["Eigengenes"]]      <- function(acohort) { paste0(Outfilepref(), "_eigengenes.txt") }
outfiles[["kME_mean"]]        <- paste0(outdir, network_group, "_mean_kME.txt")
outfiles[["gene_counts"]]     <- paste0(outdir, network_group, "_module_protein_counts.txt")
outfiles[["regression"]]      <- paste0(outdir, "nps_regression_results_per-cohort.txt")


###### MAIN ######
### Read in data ###
## Shared genes/proteins
## sort alphabetically; the proteins have been sorted that way in the TOMs
shared_genes <- sort(scan(infiles[["shared_genes"]], what = ""))

## TOMs
## These are stored in .rdata files with one object named TOM
## Read the .rdata for one cohort, save to tom to the list, then move to next cohort
toms_list <- list()   ## initiate list

for(acohort in cohorts) {
    load(infiles[["toms"]][[acohort]])
    toms_list[[acohort]] <- TOM
    rm(TOM)
}
rm(acohort)

## Expression data
expr_vec <- vector()
for(acohort in cohorts) {
    ## Read
    tmp <- fread(infiles[["expr"]][[acohort]])

    ## Extract sample names and remove sample column
    tmp_sampnames <- tmp$protsample
    tmp[, protsample := NULL]

    ## Order proteins alphabetically
    setcolorder(tmp, sort(names(tmp)))

    ## Convert to data frame and add sample names as row names
    setDF(tmp)
    rownames(tmp) <- tmp_sampnames

    ## save as named 'data' element per WGCNA specifications
    expr_vec[[acohort]] <- list(data = tmp)
}

## Phenotype data
phenos_list <- sapply(infiles[["phenos"]], fread, simplify = FALSE)


### Build consensus network -- Method 1 ###
### Start from cohort-specific expression data
### Use blockwiseConsensusModules function to build consensus modules in one step

### Build consensus modules

net <- blockwiseConsensusModules(
    multiExpr         = expr_vec,
    power             = softpower,
    corType           = "bicor",
    networkType       = "signed",
    TOMType           = "signed",
    maxBlockSize      = 10000,
    maxPOutliers      = 0.10,
    minModuleSize     = 20,
    reassignThreshold = 0,
    deepSplit         = 4,
    mergeCutHeight    = 0.3,
    numericLabels     = FALSE,
    pamRespectsDendro = FALSE
)


### Get module assignments per protein
moduleColors <- net$colors
module_assignments  <- cbind(moduleColors)
module_assignments <- data.table(protein = rownames(module_assignments),
                                 module = as.character(module_assignments))
nmodules <- length(unique(moduleColors))


### Plot gene dendrogram colored by module
if(0) {
plotDendroAndColors(dendro = net$dendrograms[[1]],
                    colors = moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
}

### Get eigengenes for each donor -- List with element per cohort
## The moduleEigengenes function returns a list with 12 elements
MEList_per_cohort <- sapply(cohorts, function(acohort) {
    moduleEigengenes(expr_vec[[acohort]]$data, colors = moduleColors)
}, simplify = FALSE)

## Extract the module eigengenes for each cohort
MEs_per_cohort <- sapply(MEList_per_cohort, function(x) {
    MEs0 <- x$eigengenes
    return(MEs0)
}, simplify = FALSE)


### Get kME -- List with element per cohort
## Calculate kME for each protein to every eigengene
kMEs_all_modules_list <- sapply(cohorts, function(acohort) {
    tmp <- signedKME(datExpr = expr_vec[[acohort]]$data, datME = MEs_per_cohort[[acohort]])
    data.table(protein = rownames(tmp), tmp)
    }, simplify = FALSE)

## Extract kME for module to which the protein was assigned
kMEs_list <- sapply(cohorts, function(acohort) {
    currdat <- kMEs_all_modules_list[[acohort]]  ## extract for this cohort
    currdat <- merge(module_assignments, currdat, by = "protein") ## add module column
    setkey(currdat, protein)

    for(amodule in unique(currdat$module)) {
        this_kme <- paste0("kME", amodule)  ## define which kME to get for these proteins
        currdat[module == amodule, kME := get(this_kme)]  ## get that kME
    }

    return(currdat[, .(protein, module, kME)])  ## keep only necessary columns
}, simplify = FALSE)


### Get weighted average kME to identify hub genes
### and make table of module gene counts
## get sample count per cohort to use as weights
sampcounts <- sapply(expr_vec, function(x) { nrow(x$data) })
weights <- sampcounts / sum(sampcounts)

## get proteins and modules from the first cohort -- they're the same across cohorts
## drop grey 'module'
kMEs_merged <- kMEs_list[[1]][module != "grey", .(protein, module)]

## add kME column for each cohort
for(acohort in cohorts) {
    kMEs_merged[, paste0("kME_", acohort) := kMEs_list[[acohort]][kMEs_merged$protein, kME]] }


## Calculate weighted average
kMEs_merged[, kME_avg := (kME_OHSU*weights["OHSU"]) +
                     (kME_rush*weights["rush"]) +
                     (kME_emory*weights["emory"])]

## Make module gene count table
module_gene_counts <- kMEs_merged[,
                                  list(total_genes = .N,
                                       hub_genes = sum(kME_avg > 0.7))
                                      ,
                                       by = "module"]
setorder(module_gene_counts, -total_genes)
module_gene_counts[, modnum := 1:nrow(module_gene_counts)]
setcolorder(module_gene_counts, "modnum")


### Test association between NPS and module eigengenes
## Set up object to hold results
## Nested list one list per cohort, each is a list of results per nps
allsumms_list <- vector(mode = "list")
for(acohort in cohorts) { allsumms_list[[acohort]] <- vector(mode = "list")}

## Iterate over cohorts and NPS
for(acohort in cohorts) {
    for(anps in nps) {
        aform <- paste0(anps, "~ mod.eigen + msex")
        atype <- "binomial"

        MEs_thiscohort <- MEs_per_cohort[[acohort]]
        ## run regression for each module (exclude grey)
        testmodules <- setdiff(colnames(MEs_thiscohort), "MEgrey")
        regression.results <- lapply(testmodules, function(amodule) {
            regression_data <- PrepareDataME(amodule, phenos_list[[acohort]], MEs_thiscohort)
            res.onemodule <- glm(form = aform, data = regression_data, family = "binomial")
            return(res.onemodule)
        })

        ## Summarize results for this outcome
        regression.summ <- data.frame(coeff.stderr = t(sapply(regression.results,
                                                              GetCoefficientAndSE, type=atype, var="mod.eigen")),
                                      pvalues = sapply(regression.results, GetPvalue, type=atype, var="mod.eigen"))
        colnames(regression.summ) <- c("Estimate", "Std.Err", "Pvalue")
        setDT(regression.summ)
        regression.summ[, eigengene :=  testmodules] ## add column with eigengene name
        setcolorder(regression.summ, "eigengene")

        ## Clean up
        print(paste("model", anps, "--", "Cleaning up"))
        allsumms_list[[acohort]][[anps]] <- regression.summ
        rm(list=c("regression.results", "regression.summ", "aform", "atype"))
    }
    names(allsumms_list[[acohort]]) <- nps
}


### Prepare flat table of association results
## First collapse across NPS within each cohort
allsumms_list2 <- sapply(allsumms_list, function(x) { rbindlist(x, idcol = "NPS") }, simplify = FALSE)

## Then collapse across cohorts to get a single data table of results (long format)
allsumms_long <- rbindlist(allsumms_list2, idcol = "cohort")
setorder(allsumms_long, eigengene, NPS)
setcolorder(allsumms_long, c("eigengene", "NPS"))


### FINISH ###
for(acohort in cohorts) {
    ## Write gene connectivity stats (kMEs)
    write.table(kMEs_list[[acohort]], file = outfiles[["Gene_conn_stats"]](acohort),
                row.names = FALSE, quote = FALSE, sep = "\t")

    ## Write module eigengenes
    ## include protsample column
    eig_for_print <- data.frame(protsample = rownames(MEs_per_cohort[[acohort]]), MEs_per_cohort[[acohort]])
    write.table(eig_for_print, file = outfiles[["Eigengenes"]](acohort),
                row.names = FALSE, quote = FALSE, sep = "\t")
    rm(eig_for_print)
}

## Write average kME
write.table(kMEs_merged, file = outfiles[["kME_mean"]], row.names = FALSE, quote = FALSE, sep = "\t")

## Write module gene counts
write.table(module_gene_counts, file = outfiles[["gene_counts"]], row.names = FALSE, quote = FALSE, sep = "\t")

## Write regression per-cohort results
write.table(allsumms_long, file = outfiles[["regression"]], row.names = FALSE, quote = FALSE, sep = "\t")
