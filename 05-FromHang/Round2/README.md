---
output:
  pdf_document: default
  html_document: default
---
# Round 2 diagnostics: instructions for running the code (Emory only)

## How to run

1. **Copy** the six `.r` files from `diagnostics_code_Round2/1_Code/` **into** the project's `1_Code/` folder.
2. **Open the Terminal tab** in RStudio (next to the Console tab).
3. **Go to the project folder** (the same place you run `03a`), for example:
   ```
   cd /path/to/31_NPS_proteomics_project_analysis
   ```
4. **Run the commands below in order**, Step 0 first, then Step 1, and so on. Wait for each command to finish before starting the next one (Step 2 commands can run in parallel if you have the cores).

Notes:

- Use `emory` (lower case) for `--cohort`.
- Packages needed: `data.table`, `argparse`, `car`, `geepack`, `clubSandwich` (`qvalue` is optional).

All output goes to a new folder, so nothing existing is overwritten:

```
2_Pipeline/03-Run_regressions/Prot_multivar_diag_R2/
```

If a run stops, run the same command again. Finished replicates are kept.

## Step 0: quick check (about 2 minutes)

```
Rscript 1_Code/03a_run_r2.r --cohort emory --maxprot 20 --outtag R2smoke
```

Open the file ending in `_diag_cohort_summary.txt` in the output folder. These three lines should be about 0 (smaller than 1e-6):

- `repro_max_abs_diff_log10P_03a` (and `..._Round1`, if shown)
- `check_max_diff_log10P_vs_pHTZ`
- `check_max_diff_log10P_chisq_vs_cr2`

If any of them is not about 0, **stop** and send us those three lines.

## Step 1: main run (about 1 to 2.5 hours)

```
Rscript 1_Code/03a_run_r2.r --cohort emory --residset ignoreSVs --stderr_variants jack,j1s,fij
```

## Step 2: permutation runs (about 3 hours each; 20 replicates of 500 proteins)

Run the three commands one after another (or in parallel if you have the cores).

```
Rscript 1_Code/03a_perm_null_r2.r --cohort emory --design fixedx --B 20 --nprot 500
Rscript 1_Code/03a_perm_null_r2.r --cohort emory --design all --protect yes --B 20 --nprot 500
Rscript 1_Code/03a_perm_null_r2.r --cohort emory --design all --protect no  --B 20 --nprot 500
```

## Step 3: runs with the existing residual files (about 1 to 2.5 hours each)

```
Rscript 1_Code/03a_run_r2.r --cohort emory --residset noProtect --stderr_variants none
Rscript 1_Code/03a_run_r2.r --cohort emory --residset regressSVs --stderr_variants none
```

## Step 4: summary (about 1 minute)

```
D=2_Pipeline/03-Run_regressions/Prot_multivar_diag_R2
Rscript 1_Code/summarize_calibration_r2.r --out $D/summary/step4 \
   --perm 5a=$D/Perm5a_fixedx \
   --perm 5b_protect=$D/Perm5b_all_protectyes \
   --perm 5b_noprot=$D/Perm5b_all_protectno \
   --full base=$D/R2_IgnoreSVs/prot_multivar_ignoreSVs_emory_R2_regression_stats_diag.txt
```

## Optional steps (only if time allows)

Within-batch permutation (about 3 hours):

```
Rscript 1_Code/03a_perm_null_r2.r --cohort emory --design batch --protect yes --B 20 --nprot 500
```

Other model settings (about 1 to 2.5 hours each):

```
Rscript 1_Code/03a_run_r2.r --cohort emory --corstr exchangeable --stderr_variants none
Rscript 1_Code/03a_run_r2.r --cohort emory --corstr independence --stderr_variants none
Rscript 1_Code/03a_run_r2.r --cohort emory --scalefix --stderr_variants none
Rscript 1_Code/03a_run_r2.r --cohort emory --xtransform rankint --stderr_variants none
```

## What to send back

Please send these small files (they contain no individual-level data). Everything is under `2_Pipeline/03-Run_regressions/Prot_multivar_diag_R2/`.

| From | Files |
|---|---|
| Each run folder (Steps 0, 1, 3 and optional) | `*_diag_cohort_summary.txt`, `*_diag_once.txt`, `*_refs_qq.png`, `*_diag_qq.png` |
| Each permutation folder (Step 2 and optional) | `*_manifest.txt`, `*_null_lambda_by_replicate.txt`, `*_null_pooled_summary.txt`, `*_null_qq.png` |
| `summary/` (Step 4) | `step4_lambda_table.txt`, `step4_Q_ratio.txt`, `step4_by_d2_tertile.txt`, `step4_perdomain.txt`, `step4_scaled_chisq.txt`, `step4_qq.png` |

Please do not send the large `*_regression_stats_diag.txt` or `*_repNNN.txt` files unless we ask for them.

If any command gives an error, send us the last 20 lines of the console output.
