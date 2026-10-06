library(data.table)
library(argparse)
library(ggplot2)

rm(list = ls())


###### GOAL ######
## Draw a forest plot of the differential expression results with FDR < 10%
## for anxiety from the last visit NPS logistic regressions.


### Constants
options(stringsAsFactors = FALSE)
cohorts <- rev(c("OHSU", "rush", "emory"))
sets <- factor(c("meta", cohorts), ordered = TRUE)


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######

### Files that exist
indir <- "2_Pipeline/03-Run_regressions_main_resid2/"
infiles <- list(
    meta = paste0(indir, "Metal/metal_qval_anxiety.txt"),
    OHSU = paste0(indir, "Stats/regression_stats_main_resid2_OHSU_anx.txt"),
    rush = paste0(indir, "Stats/regression_stats_main_resid2_rush_anx.txt"),
    emory = paste0(indir, "Stats/regression_stats_main_resid2_emory_anx.txt")
    )


### Files to be created
outdir <- "2_Pipeline/Manuscript/"
MyMkdir(outdir)

outfiles <- list(plotpref = paste0(outdir, "figure_anxiety_forest"))


###### MAIN ######
### Read in data
dats <- sapply(infiles, fread, simplify = FALSE)


### Prepare data
## Harmonize column names
setnames(
    dats[["meta"]],
    c("MarkerName", "Effect", "StdErr"),
    c("protein", "Estimate", "SE")
)

## Add gene column
for(aset in sets) {
    dats[[aset]][, gene := sapply(protein, function(x) { unlist(strsplit(x, split = "\\|"))[1] })]
}
rm(aset)

### Identify proteins to plot
## Set as factor to preserve order by p-value
genes_to_plot <- factor(dats[["meta"]][qvalue < 0.1, gene])


### Make long table of results for proteins of interest
longdat <- rbindlist(sapply(dats, function(currdat) {
    currdat[gene %in% genes_to_plot, .(gene, Estimate, SE, P, qvalue)]
}, simplify = FALSE), idcol = "set")
longdat[, set := factor(set, levels = sets)]
longdat[, gene := factor(gene, levels = genes_to_plot)]
setorder(longdat, gene, set)


### Compute the 95% confidence interval limits
longdat[, `:=`(
  lower = Estimate - (1.96*SE),
  upper = Estimate + (1.96*SE)
)]


### Add columns for OR, SE(OR), and 95% CI of OR
### Approximate using delta method formula
longdat[, OR := exp(Estimate)]
longdat[, SE_OR := exp(Estimate) * SE]

longdat[, `:=`(
  lower_OR = OR - (1.96*SE_OR),
  upper_OR = OR + (1.96*SE_OR)
)]


### Build the forest plot
### Set up
## define alternating background colors for panels
bg <- data.frame(
  gene = genes_to_plot,
  fill = rep(c("white", "grey95"),
             length.out = length(genes_to_plot))
)

## define point colors, shapes, and sizes
colorvec <- c(meta  = "#D55E00",   # strong orange/red for meta
              emory = "#9ecae1",  # blue for cohorts
              rush  = "#6baed6",
              OHSU  = "#3182bd")
shapevec <- c(meta  = 18,
              emory  = 15,
              rush  = 15,
              OHSU = 15)
sizevec <- c(meta  = 04,
             emory  = 02,
             rush  = 02,
             OHSU = 02)

## Get sample sizes
sampsizes <- sapply(cohorts, function(acohort) { max(dats[[acohort]]$nobs) })


## define set labels
labelvec <-  c(meta   = "Meta-analysis",
               OHSU   = paste0("OHSU (N=", sampsizes["OHSU"], ")"),
               rush   = paste0("RADC (N=", sampsizes["rush"], ")"),
               emory  = paste0("Emory (N=", sampsizes["emory"], ")"))

### Part 1
p <- ggplot(longdat,
            aes(x = Estimate,
                y      = set,
                colour = set,
                shape  = set
                )) +

    ## background fill
    geom_rect(data = bg,
              aes(xmin = -Inf, xmax = Inf,
                  ymin = -Inf, ymax = Inf,
                  fill = fill),          # use the colour from `bg$fill`
              inherit.aes = FALSE) +    # don’t inherit the main aes mapping
    scale_fill_identity() +

    ## vertical lines at null and other reference points
    geom_vline(xintercept = 0, linetype = "solid", colour = "grey50") +
    geom_vline(xintercept = setdiff(seq(-1, 1, 0.5), 0), linetype = "dashed", colour = "grey80") +

    ## horizontal error bar = 95% CI
    geom_errorbarh(aes(xmin = lower,
                       xmax = upper),
                   height = 0.2) +

    ## point for the estimate
    geom_point(aes(size = set))


### Part 2
p <- p +
    scale_colour_manual(breaks = sets, values = colorvec, labels = labelvec) +
    scale_shape_manual(values = shapevec, breaks = rev(sets)) +
    scale_size_manual(values = sizevec, breaks = rev(sets)) +

    scale_y_discrete(labels = function(x) { x }) +

    ## facet by gene
    facet_wrap(~ gene,
               ncol = 1,
               strip.position = "left",
               scales = "free_y") +

    ## tidy up the axis and labels
    labs(x = "Beta estimate (95% confidence interval)",
         y = NULL,
         colour = NULL,
         title = NULL) +

    ## combine into one legend (drop shape legend)
    guides(colour = guide_legend(override.aes = list(
                                     color = rev(colorvec),
                                     shape = rev(shapevec),
                                     size = rev(sizevec)),
                                 reverse = TRUE
                                 ),
           shape = "none",
           size = "none") +
    ## theme tweaks for a clean look
    theme_minimal(base_size = 10) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0, "pt")) +

    theme(axis.text.y       = element_blank(),
          strip.text.y.left = element_text(angle = 0)) +
    theme(
        legend.position = "right"

    )
#p


###### FINISH ######
for(atype in c(".png", ".pdf")) {
    ggsave(p, file = paste0(outfiles[["plotpref"]], atype),
           dpi = 300,
           width = 5, height = 2.5)
}
rm(atype)
