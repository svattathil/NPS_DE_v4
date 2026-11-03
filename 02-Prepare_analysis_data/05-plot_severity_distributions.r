library(data.table)
library(ggplot2)
library(patchwork)



rm(list = ls())
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### GOAL ######

##### CONSTANTS #####
options(stringsAsFactors = FALSE)
source("1_Code/project_constants.r")

cohorts <- c("OHSU", "emory")


###### FUNCTIONS ######
source("~/comsv/Svattathil_Library/svattathil_functions.r")


###### SETUP ######
### Command line arguments

### Files that exist
indir <- "2_Pipeline/01-Explore_data/"
infiles <- sapply(cohorts, function(cohort) {
                paste0(indir, "phenos_cleaned_", cohort, ".txt")
                }, simplify = FALSE)


### Files to be created
outdir <- "2_Pipeline/02-Prepare_analysis_data/"
MyMkdir(outdir)

outfiles <- list()
outfiles$pdf <- paste0(outdir, "/severity_distributions.pdf")


###### MAIN ######
### Read in data
phenos_list <- sapply(infiles, fread, simplify = FALSE)


### Prepare data
phenos <- rbindlist(phenos_list, idcol = "cohort", fill = TRUE)




## Draw plot for each BPSD domain
## 3 x 3 grid per cohort
## Written by Claude 2026-10-28
## long format: one row per subject x domain
long <- melt(phenos[, c("cohort", npsvars.sev), with = FALSE],
             id.vars       = "cohort",
             measure.vars  = npsvars.sev,
             variable.name = "domain",
             value.name    = "severity")
long[cohort == "emory", cohort := "Emory"]

long <- long[!is.na(severity)]
long[, severity := factor(severity, levels = 0:3)]      # keep all 4 levels
long[, domain   := factor(domain, levels = npsvars.sev)] # panel order = vector order

plot_cohort <- function(coh, ttl = coh) {
  ggplot(long[cohort == coh], aes(x = severity)) +
    geom_bar(fill = "steelblue4", width = 0.8) +
    scale_x_discrete(drop = FALSE) +   # show empty severity levels
    facet_wrap(~ domain, ncol = 4) +
    labs(title = ttl, x = "Severity", y = "N") +
    theme_bw() +
    theme(strip.background = element_rect(fill = "grey90"),
          plot.title = element_text(face = "bold"))
}

p <- plot_cohort("OHSU") / plot_cohort("Emory")
#p



###### FINISH ######
ggsave(p, file = outfiles$pdf, width = 8, height = 7)
