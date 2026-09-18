npsvars.bin <- c(agitation = "agit", anxiety = "anx", apathy = "apa",
                      depression = "depd", irritability = "irr", psychosis = "psych",
                       sleep = "nite")
npsvars.sev <- paste0(npsvars.bin, "sev")

cohorts <- c("OHSU", "rush", "emory")




Residfile <- function(option, cohort) {
  valid_options <- c("regressSVs", "ignoreSVs")
  if (!option %in% valid_options) {
    stop(sprintf(
      "Invalid 'option': '%s'. Must be one of: %s",
      option, paste(valid_options, collapse = ", ")
    ))
  }

  setlabel <- switch(option,
    regressSVs = "BatchPmiAgeMsexSVs",
    ignoreSVs  = "BatchPmiAgeMsex"
  )

  dir <- "2_Pipeline/02-Prepare_analysis_data/"
  return(paste0(dir, "resid_regress_", setlabel, "_", cohort, ".txt"))
}
