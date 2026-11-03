# Round 4 (v2): what to run (Emory only)

**What is new compared with the earlier Round 3 / Round 4 packages.** The permutation run ("5d") is the same as before (Freedman-Lane shuffle of the protein side inside each Batch, the whole pipeline repeated, B = 100, 500 proteins), but now **the seven single-domain tests are switched on** (`--perdomain`). Each test is a 1-df test of one NPS domain (agitation, anxiety, apathy, depression, irritability, psychosis, night-time) from the same GEE fit, so the permutations that calibrate the global 7-df test also calibrate the single-domain tests. Three additions:

1. `--perdomain` in every run command (the earlier packages did not have it; without it the replicate files contain no domain-specific p-values).
2. **A new summary script `summarize_perdomain_r4.r`**: lambda_GC of the domain-specific p-values (raw and permutation-adjusted), Q-Q plots per domain, direction of the effects, correlation between domains, adjusted p-values per protein and domain.
3. **Two small extra runs**: the real data with *all* proteins (Step 3) and a phenotype table (Step 4).

The output goes to **new folders** (`..._pdom_emory` and `..._pdom_emory_realall`), so nothing from earlier runs is touched, moved or overwritten. If you already ran an earlier Round 3/4 version for Emory, you do not need to delete it; the 7-df results of this run should reproduce it (same seed).

Run time: the seven extra tests add time to every protein fit, so expect **longer than the earlier run** (about 3 hours for Emory on 4 cores without them). Step 1 gives a first idea. **Please tell me how long Steps 1 and 2 took.** Choose Option A or Option B, not both.

## Step 0. Copy the files (1 minute)

Copy these files from `1_Code/` of this package **into the project's `1_Code/` folder**. The first five are the same as in Round 3; **please copy all of them, and replace `run_5d_4cores.sh`** (it has changed).

| File | What it is |
|---|---|
| `03a_perm_null_r3.r` | the permutation run (unchanged) |
| `fl_functions_r3.r` | the 5d shuffle (unchanged) |
| `diag_functions_r3.r` | the per-protein analysis functions (unchanged) |
| `cleaning_functions_r2.r` | the cleaning step of script V (unchanged) |
| `summarize_perm_r3.r` | the summary of the global 7-df test (unchanged) |
| `summarize_perdomain_r4.r` | **new**: the summary of the single-domain tests |
| `phenotype_check_r4.r` | **new**: phenotype table (no proteins) |
| `run_5d_4cores.sh` | optional, for Option A on Mac or Linux; now with `--perdomain` and the new output tag |

Packages: the same as before (`data.table`, `argparse`, `car`, `geepack`, `clubSandwich`).

Output folders:
```
2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes_pdom_emory/          (main run)
2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes_pdom_emory_realall/    (real data, all proteins)
```
Open the **Terminal** tab in RStudio and go to the project folder (the same place where you run `03a`): `cd /path/to/31_NPS_proteomics_project_analysis`. Use `emory` in lower case. **Do not remove `--outtag` or `--perdomain`.** Every command is on one line; the full list is also in `README_commands_Round4_Emory_v2.txt`.

## Step 1. Quick check (about 2 minutes or more)

```
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 2 --nprot 30 --perdomain --outtag R4smoke_pdom_emory
Rscript 1_Code/summarize_perm_r3.r --dir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/R4smoke_pdom_emory
Rscript 1_Code/summarize_perdomain_r4.r --dir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/R4smoke_pdom_emory
```

- The first command prints a small table **"Subjects per Batch"**. **Please send me this table.** It must end with `Done: replicates 0 to 2 ...`.
- The third command must print a table with one row per domain plus a row "global 7-df (HTZ)". (With B = 2 the numbers mean nothing; we only check that it runs.)
- If you see **"No pdom_ columns"**, a message **"perdomain_error"**, an error, or a warning that a Batch has fewer than 8 subjects, **stop and send me the last 20 lines.**
- Please note how long the first command took.

## Step 2. The main run: choose Option A or Option B

Settings in both options: B = 100 shuffles, 500 proteins, seven single-domain tests on.

### Option A: four cores

We split the 101 replicates (0 to 100) into four blocks, each a separate R process. **Check first that your computer has at least 4 free cores** (Mac: `sysctl -n hw.ncpu`; Linux: `nproc`). If not, use Option B or tell me.

**A1. Mac or Linux: one command**
```
bash 1_Code/run_5d_4cores.sh
```
It starts the four blocks by itself (the first one a minute ahead of the others) and waits. Messages go to `logs_5d/block1.log` to `block4.log`.

**A2. Any computer: four Terminal tabs**

Open four Terminal tabs in RStudio, `cd` to the project folder in each, and paste one command per tab. Start Tab 1 first and wait about one minute until you see `Proteins used: 500`.

Tab 1:
```
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --perdomain --outtag Perm5d_fl_protectyes_pdom_emory --rep_from 0 --rep_to 25
```
Tab 2:
```
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --perdomain --outtag Perm5d_fl_protectyes_pdom_emory --rep_from 26 --rep_to 51
```
Tab 3:
```
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --perdomain --outtag Perm5d_fl_protectyes_pdom_emory --rep_from 52 --rep_to 77
```
Tab 4:
```
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --perdomain --outtag Perm5d_fl_protectyes_pdom_emory --rep_from 78 --rep_to 100
```

### Option B: one core
```
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --perdomain --outtag Perm5d_fl_protectyes_pdom_emory
```
On Mac or Linux, in the background (closing RStudio does not stop it):
```
nohup Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --perdomain --outtag Perm5d_fl_protectyes_pdom_emory > run_5d.log 2>&1 &
```

### Both options: if something stops
Run **the same command again**. Finished replicates are kept and skipped. 101 files `*_repNNN.txt` in the output folder means finished. **Do not mix B values or change options inside the same folder.**

## Step 3. Real data with ALL proteins (one core, one replicate)

The permutation run uses 500 random proteins. This run gives the single-domain test results of the real data for **all** proteins (needed for the adjusted p-values of all proteins and for comparing cohorts protein by protein). It can be started at any time, also while Step 2 is running if a free core is available; it takes about as long as 15-20 replicates of Step 2 on one core, so **please note how long it takes**.
```
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 1 --nprot 0 --rep_from 0 --rep_to 0 --perdomain --outtag Perm5d_fl_protectyes_pdom_emory_realall
```
It must end with `Done: replicates 0 to 0 ...`. (The message "1 of 2 replicate files exist so far" at the end is expected.) If it stops with an error, send me the last 20 lines.

## Step 4. Phenotype table (about 1 minute)
```
Rscript 1_Code/phenotype_check_r4.r --cohort emory --outdir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes_pdom_emory
```
It reads only the phenotype columns; no protein results.

## Step 5. Summaries (about 1 minute each)

After all 101 replicate files of Step 2 exist and Step 3 has finished:
```
Rscript 1_Code/summarize_perm_r3.r --dir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes_pdom_emory
Rscript 1_Code/summarize_perdomain_r4.r --dir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes_pdom_emory --realall 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes_pdom_emory_realall
```
If Step 3 is not finished yet, run the second command without `--realall ...` and run it again later with it.

## Step 6. What to send back

All files are small and contain no individual-level data. In `2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes_pdom_emory/`:

| Files | What it is |
|---|---|
| `*_manifest.txt` | the settings |
| `*_r3_table.txt`, `*_r3_by_d2.txt`, `*_r3_by_replicate.txt`, `*_r3_qq_curves.txt`, `*_r3_pvalues.txt`, `*_r3_qq_2x3.png` | the global 7-df test, as in Round 3 |
| `*_r4_domain_table.txt` | **main table of the single-domain tests**: lambda_GC real and null, raw and permutation-adjusted, share p < 0.05 |
| `*_r4_domain_qq_2x4.png` | Q-Q plots of the seven domains and the distribution of the largest lambda |
| `*_r4_domain_by_df.txt`, `*_r4_domain_by_replicate.txt` | calibration of each domain by df tertile; values per replicate |
| `*_r4_domain_signed.txt`, `*_r4_domain_corr.txt`, `*_r4_domain_cases.txt` | direction of the effects; correlation between domains; cases per domain |
| `*_r4_domain_pcal.txt`, `*_r4_domain_top.txt` | per protein and domain: effect, p-value, adjusted p-value, q-value; top 10 per domain |
| `emory_r4_pheno_domains.txt`, `emory_r4_pheno_cases_by_batch.txt`, `emory_r4_pheno_phi.txt` | phenotype table (Step 4) |

Also send: the **"Subjects per Batch"** table from Step 1, the **times** of Steps 1, 2 and 3, and the last 20 lines of the console output if there was an error. Please **do not** send the large `*_repNNN.txt` files.

## For me (not needed for running)

- Settings are those of the Emory Round 3 run: `--protect yes`, `--na zero`, z-scaling, unstructured working correlation, seed 20261011, no jackknife; the 7-df results should reproduce the earlier run of this cohort. `--perdomain` is stored in the manifest, which is why the output folder has a new tag (an old folder would be moved to `stale_...`).
- Single-domain test: CR2 sandwich + Satterthwaite t (`clubSandwich::coef_test`), one coefficient of the same GEE fit as the 7-df test.
- Permutation adjustment of a domain: calibration within tertiles of the Satterthwaite df of that domain against the pooled null (leave-one-replicate-out for the null replicates), the same rule as for the global test and d2.
- Subjects: n = 73. The 500 proteins are a random subset of this cohort's proteins (same seed, not the same proteins in the other cohorts).
