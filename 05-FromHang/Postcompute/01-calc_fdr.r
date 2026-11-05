library(data.table)
library(argparse)
library(ggplot2)
library(qvalue)
rm(list = ls())



###### GOAL ######



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
#parser$add_argument("--", type="character", help="")
#parser$add_argument("--", action="store_true", default=TRUE, help="")
#parser$add_argument("--", action="store_true", default=FALSE, help="")
#parser$add_argument("--", action="store_false", dest="", help="")
#parser$add_argument("--", type="integer", default=10, help="")
#parser$add_argument("--", type="double", default=10, help="")

## Read from parser
args <- parser$parse_args()


### Files that exist
indir <- "2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/PermP_fromHang/"
infiles <- sapply(cohorts, function(x) {
    paste0(indir, "NPS_7df_pperm_", Capwords(x), "_forQ.txt")
}, simplify = FALSE)


### Files to be created
outdir <- "2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Postcompute/"
MyMkdir(outdir)
outfiles <- list(
    txt = paste0(outdir, "pperm_fdr_allcohorts.txt"),
    png = paste0(outdir, "pperm_histograms_allcohorts.png")
)


###### MAIN ######
### Read in data
dat_list <- sapply(infiles, fread, simplify = FALSE)


### Prepare data
## Extract relevant columns to single table, for simplicity
pperm <- rbindlist(dat_list, idcol = "cohort")[, .(cohort, protein, pperm)]
pperm[, cohort := factor(cohort, levels = cohorts)]

### Make quick p-value histograms for each cohort
p_hists <- ggplot(pperm, aes(x = pperm)) +
    facet_wrap(~ cohort, nrow = 2) +
    geom_histogram(bins = 100, fill = "SteelBlue", color = "SteelBlue4") +
    theme_bw()


### Add FDR columns
for(acohort in cohorts) {
    pperm[cohort == acohort, p.BH := p.adjust(pperm, method = "BH")]
    pperm[cohort == acohort, qvalue := qvalue(pperm)$qvalues]
}
rm(acohort)

pperm[p.BH < 0.05, sig.BH05 := "sig_BH05"]
pperm[qvalue < 0.05, sig.qvalue05 := "sig_qvalue05"]
setorder(pperm, cohort, pperm)


###### FINISH ######
ggsave(file = outfiles$png, p_hists, width = 7, height = 7, unit = "in", dpi = 300)
write.table(pperm, file = outfiles$txt, row.names = FALSE, quote = FALSE, sep = "\t")
