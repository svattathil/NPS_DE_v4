README_Round2.txt  --  diagnostics_code_Round2   (Emory only)
=====================================================================================================
What this round asks
  Round 1 showed (Emory, n = 73, 7,959 proteins):  the HTZ correction over-corrects (null Q is about 1.23-1.27 x
  chi-square(7) at all quantiles, independent of d2, yet the HTZ null rejection rate rises with d2); the real Q is
  another 1.3-1.55 x the null Q at every quantile (a shift of the whole bulk).  Round 2 tests the most plausible
  explanations and generates the "san.se / jackknife" curves of the earlier analysis.

  R2-A  Reference curves: the SAME model, scored with the other tests (CR2 chi-square, Naive-F, HTA/HTB/HTZ, EDF/EDT)
        and with the san.se / jackknife covariances (chi-square Wald test); plus the 7 single-domain tests.   [Step 1]
  R2-B  Does the cleaning with protected NPS columns create the real-vs-null shift?  (permutation null that REPEATS
        the cleaning, Code 5b; the existing residual files noProtect / regressSVs)                          [Steps 2-4]
  R2-C  Does the batch-NPS association matter?  (permutation within Batch, Code 5c)                          [Step 5]
  R2-D  Is the problem the working correlation, the dispersion or the x distribution?  (exchangeable /
        independence, scale.fix, rank-inverse-normal x)                                                      [Step 6]

Folder contents (copy 1_Code/* into the project's 1_Code/; tests/ is optional)
  1_Code/cleaning_functions_r2.r     replica of the script-V cleaning (reproduces your residual file to 5e-14 in Emory)
  1_Code/diag_functions_r2.r         Round 1 diag_functions.r + the Round 2 additions (all tests, jackknife, per-domain)
  1_Code/03a_run_r2.r                one full run (Round 1 Code 1 + options)
  1_Code/03a_perm_null_r2.r          permutation null: 5a (fixedx) / 5b (all) / 5c (batch)
  1_Code/summarize_calibration_r2.r  tables and plots from any number of the runs above
  1_Code/compare_runs_r2.r           protein-by-protein comparison of two runs
  tests/test_make_synthetic_data_r2.R synthetic Emory-like data for a dry run (never use for results)
Nothing in Round 1 or in your existing files is overwritten; all output goes to
  2_Pipeline/03-Run_regressions/Prot_multivar_diag_R2/<tag>/
Prerequisites: the packages used before (data.table, argparse, car, geepack, clubSandwich); qvalue is optional.
Run everything from the project directory, as for 03a.  Use "emory" (lower case) for --cohort.

-----------------------------------------------------------------------------------------------------
Step 0 (2 min)  Smoke test
  Rscript 1_Code/03a_run_r2.r --cohort emory --maxprot 20 --outtag R2smoke
  Check <outdir>/..._diag_cohort_summary.txt:
     repro_max_abs_diff_log10P_03a / _Round1   ~ 0 (printed if the earlier summary files exist)
     check_max_diff_log10P_vs_pHTZ              ~ 0 (HTZ of Wald_test = P)
     check_max_diff_log10P_chisq_vs_cr2         ~ 0
  If any is not ~ 0, STOP and send the three lines back.

Step 1 (about 1-2.5 h; ~1 s per protein with all extras, slower with --stderr_variants jack,j1s,fij)
  Rscript 1_Code/03a_run_r2.r --cohort emory --residset ignoreSVs --stderr_variants jack,j1s,fij
  -> prot_multivar_ignoreSVs_emory_R2_regression_stats_diag.txt  and  ..._refs_qq.png
  The *_refs_qq.png has the Q-Q curves of the same 7,959 proteins under: CR2+HTZ (the analysis), CR2+chi-square,
  CR2 + Naive-F / HTA / HTB / EDF / EDT, and san.se + jackknife (jack, j1s, fij) chi-square.  These are the "earlier run"
  curves, so Selina does not need to find the old plots.  The 7 single-domain tests are the pdom_* columns.
  (For speed, --stderr_variants jack alone is enough; the default is "jack".)

Step 2 (about 8-10 min per replicate of 500 proteins;  B = 20 -> about 3 h per design; start with --B 5)
  The key test:  does the permutation null that repeats the cleaning reproduce the real-data shift?
  Rscript 1_Code/03a_perm_null_r2.r --cohort emory --design fixedx --B 20 --nprot 500                 # 5a (same 500 proteins as Round 1)
  Rscript 1_Code/03a_perm_null_r2.r --cohort emory --design all --protect yes --B 20 --nprot 500       # 5b, analysis cleaning
  Rscript 1_Code/03a_perm_null_r2.r --cohort emory --design all --protect no  --B 20 --nprot 500       # 5b, noProtect cleaning
  Replicate 0 of each is the real data processed the same way; with --protect yes it must agree with 5a's replicate 0.
  If a run is stopped, repeat the same command: finished replicates are kept.  If the settings differ from the
  existing replicate files, the old files are moved to stale_<time>/ (never reused).

Step 3 (about 1-2.5 h each)  Real data with the existing residual files (needs the files to exist)
  Rscript 1_Code/03a_run_r2.r --cohort emory --residset noProtect --stderr_variants none
  Rscript 1_Code/03a_run_r2.r --cohort emory --residset regressSVs --stderr_variants none

Step 4 (a minute)  Summaries
  D=2_Pipeline/03-Run_regressions/Prot_multivar_diag_R2
  Rscript 1_Code/summarize_calibration_r2.r --out $D/summary/step4 \
     --perm 5a=$D/Perm5a_fixedx --perm 5b_protect=$D/Perm5b_all_protectyes --perm 5b_noprot=$D/Perm5b_all_protectno \
     --full base=$D/R2_IgnoreSVs/prot_multivar_ignoreSVs_emory_R2_regression_stats_diag.txt

Step 5 (optional, about 3 h)  within-Batch permutation
  Rscript 1_Code/03a_perm_null_r2.r --cohort emory --design batch --protect yes --B 20 --nprot 500

Step 6 (optional, about 1-2.5 h each)  working correlation / dispersion / x transform
  Rscript 1_Code/03a_run_r2.r --cohort emory --corstr exchangeable --stderr_variants none
  Rscript 1_Code/03a_run_r2.r --cohort emory --corstr independence --stderr_variants none
  Rscript 1_Code/03a_run_r2.r --cohort emory --scalefix --stderr_variants none
  Rscript 1_Code/03a_run_r2.r --cohort emory --xtransform rankint --stderr_variants none
  Protein-by-protein comparison with the baseline:
  Rscript 1_Code/compare_runs_r2.r --a <baseline stats_diag> --b <other stats_diag> --labela base --labelb other --out $D/cmp_x

What to send back (small files only; no individual-level data are in them)
  from <outdir> of each run:  *_diag_cohort_summary.txt, *_diag_once.txt, *_refs_qq.png, *_diag_qq.png
  from each permutation folder:  *_manifest.txt, *_null_lambda_by_replicate.txt, *_null_pooled_summary.txt, *_null_qq.png
  from summary/:  step4_lambda_table.txt, step4_Q_ratio.txt, step4_by_d2_tertile.txt, step4_perdomain.txt,
                  step4_scaled_chisq.txt, step4_qq.png
  (the large *_regression_stats_diag.txt and *_repNNN.txt files are only needed on request)

How to read the results (what each outcome would mean; none of these is established yet)
  * real_over_null in step4_Q_ratio.txt is the bulk shift.  Compare it across 5a / 5b_protect / 5b_noprot:
      - 5b_protect null Q is larger than 5a's and close to the real Q  -> the cleaning with protected NPS columns
        explains the shift; 5b_noprot null Q ~ 5a and the real noProtect run (Step 3) has a smaller real Q
        -> use a cleaning that does not protect the NPS outcomes, or calibrate with 5b.
      - 5b ~ 5a  -> the shift is not from the cleaning; look at 5c (batch) and at Step 6.
      - 5c null Q larger than 5b  -> the batch-NPS association contributes; adjust for Batch in the permutation scheme.
  * null_over_chisq7 ~ 1.23-1.27 at all quantiles in every design -> the over-correction of HTZ is a separate issue;
    the p_chisq_* / NaiveF / HTA / EDF columns show which reference is calibrated (Step 1 and by_d2_tertile).
  * The san.se / jackknife curves (Step 1): compare lambda_GC of p_chisq_jack, p_chisq_j1s, p_chisq_fij with the CR2 ones.
  * *_scaled_chisq.txt is a post-hoc recalibration for orientation only (the c comes from the same null); it is not a
    valid analysis, and its BH counts are not trustworthy if the shift is not explained.

Dry run without real data (optional, any scratch folder with a copy of 1_Code and ~/comsv/Svattathil_Library):
  Rscript tests/test_make_synthetic_data_r2.R      # then the Step 0 command, and Step 2 with --B 2 --nprot 30
  This only tests that the code runs; the numbers mean nothing.
