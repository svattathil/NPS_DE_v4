npsvars.bin <- c(agitation = "agit", anxiety = "anx", apathy = "apa",
                      depression = "depd", irritability = "irr", psychosis = "psych",
                       sleep = "nite")
npsvars.sev <- paste0(npsvars.bin, "sev")

cohorts <- c("OHSU", "rush", "emory")




Residfile <- function(option, cohort) {
    valid_options <- c("regressSVs", "ignoreSVs", "noProtect", "resid2_multi")
    if (!option %in% valid_options) {
        stop(sprintf(
            "Invalid 'option': '%s'. Must be one of: %s",
            option, paste(valid_options, collapse = ", ")
        ))
    }

    if(option %in% c("regressSVs", "ignoreSVs", "noProtect")) {
        setlabel <- switch(option,
                           regressSVs = "BatchPmiAgeMsexSVs",
                           ignoreSVs  = "BatchPmiAgeMsex",
                           noProtect  = "BatchPmiAgeMsex_noProtect"
                           )

        dir <- "2_Pipeline/02-Prepare_analysis_data/"
        return(paste0(dir, "resid_regress_", setlabel, "_", cohort, ".txt"))
    }

    if(option %in% c("resid2_multi")) {
        dir <- paste0("Older_versions/2_Pipeline_v3/02-Prepare_analysis_data/",
                      "Resid2_regressSVs_starting_from_resid1/")
        return(paste0(dir, "resid_regressBatchPmiAgeSVs_", cohort, "_multi.txt"))
    }
}
