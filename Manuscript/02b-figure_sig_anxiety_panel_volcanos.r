library(data.table)
library(argparse)
library(ggplot2)
library(scales)
library(patchwork)
library(stringr)
rm(list = ls())


###### GOAL ######
## Draw a forest plot of the differential expression results with FDR < 10%
## for anxiety from the last visit NPS logistic regressions.


### Constants
options(stringsAsFactors = FALSE)
cohorts_code <- c("OHSU", "rush", "emory")
cohorts <- c("OHSU", "RADC", "Emory")
sets <- factor(c(cohorts, "Meta"), ordered = TRUE)


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######

### Files that exist
indir <- "2_Pipeline/03-Run_regressions/Prot_perDomain/Basic/"
infiles <- list()
for(cohort in cohorts_code) {
infiles[[cohort]] <- paste0(indir, cohort, "_basic.tsv")
}
rm(cohort)
names(infiles) <- cohorts ## update with cohort names to be used in plot

infiles[["Meta"]] <- paste0(indir, "Metafor/", "metafor_basic_using_SE_mod.txt")


### Files to be created
outdir <- "2_Pipeline/Manuscript/"
MyMkdir(outdir)

outfiles <- list(plotpref = paste0(outdir, "figure_anxiety_volcanos"))


###### MAIN ######
### Read in data
dats_all <- sapply(infiles, fread, simplify = FALSE)


### Prepare data
## Harmonize column names
setnames(dats_all[["Meta"]],
    c("fe_beta", "fe_se", "fe_p"),
    c("LogFC", "SE", "P")
    )

## Add gene column
for(aset in sets) {
    dats_all[[aset]][, gene := sapply(protein, function(x) {
        unlist(strsplit(x, split = "\\|"))[1] })]
}
rm(aset)

for(cohort in cohorts) {
    setnames(dats_all[[cohort]],
             c("domain", "logFC", "SE_mod", "P.Value", "fdr"),
             c("nps", "LogFC", "SE", "P", "qvalue"))
}
rm(cohort)


## Restrict to anxiety results
dats <- sapply(dats_all, function(dat) { dat[nps == "anx", ] })


### Identify proteins to plot
## Set as factor to preserve order by p-value
genes_to_plot <- factor(dats[["Meta"]][qvalue < 0.05, gene])


### Make long table of results for all proteins
meta_label <- "Cross-cohort Meta-analysis" ## label to use in the plot
longdat <- rbindlist(sapply(dats, function(currdat) {
    currdat[, .(gene, LogFC, SE, P, qvalue)]
}, simplify = FALSE), idcol = "set")
longdat[set == "Meta", set := meta_label]
longdat[, set := factor(set, levels = c(cohorts, meta_label))]
longdat[, sig := ifelse(qvalue < 0.05, TRUE, FALSE)]


### Build the volcano plots
### Set up
dt <- copy(longdat)
dt[, negLog10P := -log10(P)]  ## add -log10(P) column

## Make one plot per set
xlims <- c(min(dt$LogFC), max(dt$LogFC))
ylims <- c(0, max(dt$negLog10P))

plotlist <- sapply(unique(longdat$set), function(aset) {
    ggplot(dt[set == aset, ],
           aes(x = LogFC,
               y = negLog10P,
               colour = sig)) +               ## colour by significance
        geom_point(alpha = 0.7, size = 1.2) +
        geom_vline(xintercept = 0,
                   linetype = "dashed",
                   colour = "grey70") +
        scale_colour_manual(values = c("FALSE"  = "grey60",
                                       "TRUE" = "#D55E00"),
                            name = "Significant") +
        scale_x_continuous(breaks = function(lim) {
            seq(floor(lim[1]), ceiling(lim[2]), by = 0.5)
        }) +
        coord_cartesian(xlim = xlims, ylim = ylims) +
        labs(title = stringr::str_wrap(aset, 13),
             x     = "LogFC",
             y     = expression(-log[10]~"p-value")) +
        theme_minimal(base_size = 10) +
        theme(legend.position = "none",
              plot.title      = element_text(hjust = 0, vjust = 0),
              panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.5),
              panel.spacing   = unit(0.8, "lines")
              )
}, simplify = FALSE)


### Modify to prepare for plotting together
## Turn off x-axis title for plots 1 and 3
plotlist[[1]] <- plotlist[[1]] + theme(axis.title.x = element_blank())
plotlist[[3]] <- plotlist[[3]] + theme(axis.title.x = element_blank())

## Turn off y-axis title and labels for plots 2 and 3
plotlist[[2]] <- plotlist[[2]] + theme(axis.title.y = element_blank(),
                                       axis.text.y = element_blank())
plotlist[[3]] <- plotlist[[3]] + theme(axis.title.y = element_blank(),
                                       axis.text.y = element_blank())


## Add extra margin to right of plot 4
plotlist[[4]] <- plotlist[[4]] + theme(plot.margin = margin(0,0,0,20))


### Draw with patchwork
paneled <- (plotlist[[1]] | plotlist[[2]] | plotlist[[3]]) + plotlist[[4]]


###### FINISH ######
for(atype in c(".png", ".pdf")) {
    ggsave(paneled, file = paste0(outfiles[["plotpref"]], atype),
           dpi = 300,
           width = 7, height = 3)
}
rm(atype)


