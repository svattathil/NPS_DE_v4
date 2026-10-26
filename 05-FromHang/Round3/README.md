---
output:
  pdf_document: default
  html_document: default
---
# Round 3: what to run (Emory only)

**One thing to run: a permutation test called "5d". Everything else is a short check.**
It takes about **4 hours on 4 cores** (Option A) or about **16 hours on 1 core** (Option B). Choose one option, not both.

<!--
## What is 5d? (read once, 1 minute)

In Round 2 we shuffled the NPS vectors between subjects. This time we shuffle the **protein side** instead:

1. For each protein we remove what Batch, PMI, age and sex explain (a normal regression).
2. What is left over is shuffled between subjects **inside the same Batch**.
3. We add the explained part back, and run **exactly the same pipeline as before** (cleaning of script V, z-scaling, GEE, all tests).

The NPS data are **not** touched, so their relation to Batch, PMI, age and sex stays as in the real data.
The same shuffle is used for all proteins of one replicate. Replicate 0 is the real data.

We use B = 100 shuffles and the same 500 proteins as in Round 2. The slow extras of Round 2 (jackknife, single-domain tests) are switched off.
-->

## Step 0. Copy the files (1 minute)

Copy these files from `diagnostics_code_Round3/1_Code/` **into the project's `1_Code/` folder**:

| File | What it is |
|---|---|
| `03a_perm_null_r3.r` | the permutation run (new) |
| `fl_functions_r3.r` | the 5d shuffle (new) |
| `diag_functions_r3.r` | same as the Round 2 file, only faster (new) |
| `cleaning_functions_r2.r` | **unchanged** Round 2 file. If you still have it in `1_Code/`, you do not need to copy it again |
| `summarize_perm_r3.r` | makes the tables and the figure (new) |
| `run_5d_4cores.sh` | optional, for Option A on Mac or Linux |

Packages: the same as before (`data.table`, `argparse`, `car`, `geepack`, `clubSandwich`). Nothing existing is overwritten. All output goes to a new folder:

```
2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/
```

How to run (every time): open the **Terminal** tab in RStudio, then go to the project folder (the same place where you run `03a`):

```
cd /path/to/31_NPS_proteomics_project_analysis
```

Use `emory` in lower case.

## Step 1. Quick check (about 2 minutes)

```
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 2 --nprot 30 --outtag R3smoke
```

Look at the console output:

- It prints a small table **"Subjects per Batch"**. **Please send me this table.**
- It must end with `Done: replicates 0 to 2 ...`.
- If there is an error, or a warning that a Batch has fewer than 8 subjects, **stop and send me the last 20 lines**.

## Step 2. The main run: choose Option A or Option B

Settings are the same in both options: B = 100 shuffles, 500 proteins.

### Option A: four cores (about 4 hours)

We split the 101 replicates (0 to 100) into four blocks. Each block is a separate R process. **Check first that your computer has at least 4 free cores** (Mac: `sysctl -n hw.ncpu`; Linux: `nproc`). If not, use Option B or tell me.

**A1. Mac or Linux: one command**

```
bash 1_Code/run_5d_4cores.sh
```

It starts the four blocks by itself (the first one a minute ahead of the others) and waits. Messages of each block go to `logs_5d/block1.log` to `block4.log`.

**A2. Any computer (Mac, Linux or Windows): four Terminal tabs**

Open **four** Terminal tabs in RStudio (the **+** button next to the Terminal tab). In each tab go to the project folder (`cd ...`). Then paste **one** command per tab:

Tab 1 (start this one first, wait about one minute until you see `Proteins used: 500`):
```
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --rep_from 0 --rep_to 25
```
Tab 2:
```
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --rep_from 26 --rep_to 51
```
Tab 3:
```
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --rep_from 52 --rep_to 77
```
Tab 4:
```
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --rep_from 78 --rep_to 100
```

(Windows: the commands are the same. Type them on one line, without any `\`.)

When all four tabs print `Done:`, go to Step 3.

### Option B: one core (about 16 hours)

One command, one tab:

```
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500
```

Tips: start it in the evening. Keep the computer awake. On Mac or Linux you can run it in the background so that closing RStudio does not stop it:

```
nohup Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 > run_5d.log 2>&1 &
```

### Both options: if something stops

Run **the same command again**. Finished replicates are kept and skipped. To check the progress, count the files `*_repNNN.txt` in the output folder: 101 files means finished.

If you would like a shorter run first, use `--B 50` instead of `--B 100` in every command (about half the time). **Do not mix B values** in the same folder.

## Step 3. Summary (about 1 minute)

Run this **after all 101 replicate files exist**:

```
Rscript 1_Code/summarize_perm_r3.r \
   --dir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes \
   --full 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R2/R2_IgnoreSVs/prot_multivar_ignoreSVs_emory_R2_regression_stats_diag.txt
```

(On Windows put it on one line and remove the `\`.) The `--full` file is the Round 2 main run; it gives the red dots of all proteins. If you no longer have it, leave the `--full` part out.

**Check:** in `Perm5d_fl_protectyes_r3_table.txt`, the row `HTZ` should have `lambda_real` about **1.87** (the same as replicate 0 of Round 2). If it is not close to 1.87, **stop and tell me**.

## Step 4. What to send back

All files are small and contain no individual-level data. 

They are in
`2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes/`:

| File | What it is |
|---|---|
| `*_manifest.txt` | the settings |
| `*_r3_table.txt` | **main table**: one row per reference method |
| `*_r3_by_d2.txt` | the same, by tertile of the denominator df |
| `*_r3_by_replicate.txt` | one row per replicate and method |
| `*_r3_qq_curves.txt` | numbers of the Q-Q curves |
| `*_r3_pvalues.txt` | p-values of all proteins and replicates (about 5 MB) |
| `*_r3_qq_2x3.png` | the Q-Q figure (2 x 3 panels) |

Also send:
- the **"Subjects per Batch"** table from Step 1,
- the **time** the run took (the last lines of the log show it),
- if there was an error: the last 20 lines of the console output.

Please **do not** send the large `*_repNNN.txt` files. <!-- unless we ask. -->

## For me (not needed for running)

- Replicate 0 reproduces the real data: the script stops with an error if the rebuilt expression differs from the real one by more than 1e-8.
- The settings are the analysis settings: `--protect yes`, `--na zero` (missing set to 0 in the cleaning, as in script V), z-scaling, unstructured working correlation. The effect of missing values is a later question and is **not** tested in Round 3.
- Seed 20261011 and the same 500 proteins as Round 2.
- Test on synthetic data (never use for results): replicate 0 of 5d equals replicate 0 of 5b to 1e-15, and the p-values of `diag_functions_r3.r` equal those of `diag_functions_r2.r`.
