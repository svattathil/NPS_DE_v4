DIAGNOSTIC CODES 1-5a   (plan: NPS_Analysis_rev4.tex, Section 3)                    11 Oct 2026  (v2)
=============================================================================================
This folder supersedes diagnostics_code/ (v1). Put all .r files in 1_Code/ next to 03a-run_multivariate.r.

  diag_functions.r             helpers; contains the per-protein analysis of 03a
  03a_code1_diagnostics.r      Code 1; Code 2 (--residset regressSVs); Code 3 runs (--residfile);
                               Code 4 (--corstr, --scalefix)
  03a_code5a_permutation.r     Code 5a (permutation null, test only)
  01a_code3_cleaning_variants.r  Code 3: new residual files (missing-aware cleaning, noProtect, replica) + bookkeeping
  compare_runs.r               compares two runs protein by protein (Code 1 vs Code 2, 3, 4)
  test_make_synthetic_data.R   generator of a small synthetic cohort (used by us to test the code)

Nothing of the earlier results is overwritten. Output of the 03a-type scripts goes to
  2_Pipeline/03-Run_regressions/Prot_multivar_diag/<tag>_<Residset>/
01a_code3_cleaning_variants.r writes NEW files in 2_Pipeline/02-Prepare_analysis_data/ and stops if a
file already exists.

---------------------------------------------------------------------------------------------
STEP 0. Smoke test (a few minutes), from the project directory as for 03a. Please send us the output.
  Rscript 1_Code/03a_code1_diagnostics.r  --cohort emory --maxprot 20
  Rscript 1_Code/03a_code5a_permutation.r  --cohort emory --B 2 --nprot 30
(cohort names as used by 03a/Residfile: OHSU, rush, emory). With the earlier 03a summary file present,
Code 1 prints "Reproduction check: max |diff log10 P| = ..." (should be ~1e-6 or smaller).

STEP 1. Code 1 (about 25% more time than the original run):   [per cohort]
  Rscript 1_Code/03a_code1_diagnostics.r --cohort emory --residset ignoreSVs

STEP 2. Code 5a (about (B+1) x nprot fits; the time per fit is printed and stored in elapsed_sec):
  Rscript 1_Code/03a_code5a_permutation.r --cohort emory --B 20 --nprot 1000
  (if the job stops, run the same command again: finished replicates are skipped)

STEP 3. Code 2:
  Rscript 1_Code/03a_code1_diagnostics.r --cohort emory --residset regressSVs

STEP 4. Code 3 (cleaning variants).
  3b bookkeeping only (fast, no new residual file):
    Rscript 1_Code/01a_code3_cleaning_variants.r --cohort emory --protect yes --na zero --bookkeeping_only
  The check of the cleaning function (our replica of jaffelab::cleaningY against your existing file):
    Rscript 1_Code/01a_code3_cleaning_variants.r --cohort emory --protect yes --na zero
    -> prints "Replica of cleaningY vs existing file: max |diff| = ...": ~1e-8 or smaller means that the
       jaffelab version used sets missing values to 0 (as the current GitHub version does).
  The new residual files (2 x 2 design; 'protect yes, na zero' and 'protect no, na zero' already exist):
    Rscript 1_Code/01a_code3_cleaning_variants.r --cohort emory --protect yes --na aware
    Rscript 1_Code/01a_code3_cleaning_variants.r --cohort emory --protect no  --na aware
  Then Code 1 on each cell (the 'noProtect' file already exists):
    Rscript 1_Code/03a_code1_diagnostics.r --cohort emory --residset noProtect
    Rscript 1_Code/03a_code1_diagnostics.r --cohort emory --residset naAware \
        --residfile 2_Pipeline/02-Prepare_analysis_data/resid_regress_BatchPmiAgeMsex_naAware_emory.txt
    Rscript 1_Code/03a_code1_diagnostics.r --cohort emory --residset noProtect_naAware \
        --residfile 2_Pipeline/02-Prepare_analysis_data/resid_regress_BatchPmiAgeMsex_noProtect_naAware_emory.txt

STEP 5. Code 4 (working model; each a separate run):
    Rscript 1_Code/03a_code1_diagnostics.r --cohort emory --corstr independence
    Rscript 1_Code/03a_code1_diagnostics.r --cohort emory --corstr exchangeable
    Rscript 1_Code/03a_code1_diagnostics.r --cohort emory --scalefix
    Rscript 1_Code/03a_code1_diagnostics.r --cohort emory --corstr independence --scalefix
    Rscript 1_Code/03a_code1_diagnostics.r --cohort emory --corstr exchangeable --scalefix

STEP 6. Comparison of any run with Code 1:
    Rscript 1_Code/compare_runs.r --a <Code1 ...stats_diag.txt> --b <other run ...stats_diag.txt> \
        --labela Code1 --labelb Code4_independence --out <prefix>

---------------------------------------------------------------------------------------------
WHAT TO SEND BACK (protein-level and aggregate files only; no individual-level data)
 Each 03a-type run:  *_regression_stats_diag.txt, *_diag_cohort_summary.txt, *_diag_once.txt,
                     *_diag_qq.png, *_qq.png, (*_diag_errors.txt if some fits failed)
 Code 5a:            *_null_lambda_by_replicate.txt, *_null_pooled_summary.txt, *_null_qq.png,
                     *_real_with_empirical_p.txt, *_rep###.txt (rep000 = real data)
 Code 3:             *_code3_bookkeeping.txt, *_na_per_protein.txt, *_batch_by_nps.txt (cells <5 are masked),
                     *_code3_check.txt   (not the residual files themselves, unless asked)
 compare_runs:       *_compare_summary.txt, *_compare.png

---------------------------------------------------------------------------------------------
COLUMNS ADDED TO THE PER-PROTEIN TABLE
 nobs, n_g              subjects used for the protein (in 03a the column nobs seems not to be written)
 df_num, df_denom       q and d2 of the HTZ test; eta_from_df = d2 + q - 1
 Q_cr2, Q_from_F        Wald quadratic form with the CR2 block, and Q recomputed from the HTZ F statistic
                        (rel_diff_Q ~ 0 checks our reading of the HTZ formula)
 p_chisq_cr2/cr0/san    p-values with the chi-square(7) reference and the CR2, CR0, san.se covariances
                        (CR0 and san.se coincide: a check of the two software routes)
 phi_hat, phi_pearson   dispersion from geeglm / recomputed from Pearson residuals
 min_eig_R, cond_R      smallest eigenvalue / condition number of the working correlation used by geeglm
 min_eig_Remp, cond_Remp, max_abs_rho_emp   the same for the empirical correlation of Pearson residuals
 alpha_min, alpha_max   range of the estimated working correlations
 min_cases, min_noncases, n_rare   case counts per domain (n_rare = domains with < 5 cases)
 kurt_x, max_abs_x, n_distinct_x   shape of the protein's expression (normal: kurtosis 3)
 n_na_resid             missing values of the protein in the cleaned file
 min_mu, max_mu, max_abs_beta, max_se_cr2, max_se_san   near-separation indicators
 n_warn, warn_first, geese_error   warnings and geepack's convergence value
 delta_wald             column 'delta' of Wald_test (NOT eta; kept for reference)
 elapsed_sec            seconds for this protein

TESTED by us on synthetic cohorts (73 subjects; R 4.3.3, geepack 1.3.10, clubSandwich 0.5.10):
 - Code 1 reproduces 03a to 1e-15; Code 5a runs and resumes; Codes 2/3/4 run through the same script;
 - the missing-aware cleaning agrees with an independent lm() per protein (1e-15) and equals the zero
   version for proteins without missing values; the replica of cleaningY reproduces a file made by it.
NOTE: with geepack 1.3.10 + clubSandwich 0.5.10 geeglm returns linear.predictors as an N x 1 matrix and
vcovCR stops ("non-conformable arrays"); diag_functions.r flattens it (changes no number; a no-op for
versions that already return a vector). geeglm(scale.value = 1) is not accepted; --scalefix uses scale.fix = TRUE only.
