library(data.table)
library(argparse)
library(ggplot2)
library(reshape2)
library(scales)

rm(list = ls())



###### GOAL ######




###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######
### Command line arguments
## Create parser
parser <- ArgumentParser()
parser$add_argument("--cohort", type="character", default = "OHSU", help="OHSU or emory")

## Read from parser
args <- parser$parse_args()


### Constants
options(stringsAsFactors = FALSE)
npsvars.bin <- c(agitation = "agit", anxiety = "anx", apathy = "apa", delusion = "del",
                 depression = "depd", disinhibition = "disn", hallucination = "hall",
                 irritability = "irr", sleep = "nite")

### Files that exist
infiles <- list(long = paste0("2_Pipeline/01-Explore_data/longitudinal_cleaned_", args$cohort, ".txt"),
              basic = paste0("2_Pipeline/01-Explore_data/phenos_cleaned_", args$cohort, ".txt"))


### Files to be created
outfiles <- list(plot_lastvisit = paste0("2_Pipeline/01-Explore_data/nps_corr_lastvisit_", args$cohort, ".png"))


###### MAIN ######
### Read in data
dats <- sapply(infiles, fread, simplify = FALSE)


### Prepare data
trait_mat <- as.matrix(dats[["basic"]][, ..npsvars.bin, with = FALSE])



### Calculate correlations
corr_mat <- cor(trait_mat, method = "pearson", use = "pairwise.complete.obs")

# Optional: also compute a p‑value matrix (useful for masking non‑significant links)
p_mat <- {
  n    <- nrow(trait_mat)
  pval <- matrix(NA_real_, ncol = ncol(trait_mat), nrow = ncol(trait_mat))
  colnames(pval) <- rownames(pval) <- colnames(corr_mat)
  for (i in 1:ncol(trait_mat)) {
    for (j in i:ncol(trait_mat)) {
      test <- cor.test(trait_mat[, i], trait_mat[, j],
                       method = "pearson")
      pval[i, j] <- pval[j, i] <- test$p.value
    }
  }
  pval
}

### Reshape to long format for ggplot
corr_long <- melt(corr_mat, varnames = c("Trait1", "Trait2"),
                  value.name = "Correlation")

## Merge in p-values
if (exists("p_mat")) {
  p_long <- melt(p_mat, varnames = c("Trait1", "Trait2"),
                 value.name = "p.value")
  corr_long <- merge(corr_long, p_long,
                     by = c("Trait1", "Trait2"))
}

###  Build the heat‑map
p1 <- ggplot(corr_long, aes(x = Trait1, y = Trait2, fill = Correlation)) +
    geom_tile(colour = "white") +                         # tiles with white borders

    ## Colour scale – blue (negative) → white (zero) → red (positive)
    scale_fill_gradient2(low = muted("steelblue"),
                         mid = "white",
                         high = muted("tomato"),
                         midpoint = 0,
                         limits = c(-1, 1),
                         name = "Φ‑coeff.") +

    ## Optional: mask non‑significant cells (e.g., p > .05) by turning them grey
    ## Uncomment the next two lines if you added the p‑value matrix
    geom_tile(data = subset(corr_long, p.value >= 0.05),
              fill = "grey", colour = NA) +

    ## Print the numeric value inside each tile
    geom_text(aes(label = sprintf("%.2f", Correlation)),
              size = 3, colour = "black") +

    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          legend.position = "right") +
    labs(title = paste0("Pairwise binary‑trait correlations - last visit ", args$cohort),
         subtitle = "Pearson φ‑coefficient (color scale shown if p-value < 0.05)",
         x = NULL, y = NULL)


###### FINISH ######
ggsave(p1, file = outfiles[["plot_lastvisit"]])
