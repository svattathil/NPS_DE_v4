ROUND 3: TERMINAL COMMANDS ONLY (from README.md; the explanations are in README.md)
=====================================================================================
- Use the Terminal tab in RStudio. Run every command from the project folder.
- Every command is on ONE line. Copy the whole line.
- Lines that start with # are comments, not commands.
- Use Step 2 Option A OR Option B, not both.


# ---- Before every step: go to the project folder (the same place where you run 03a) ----
cd /path/to/31_NPS_proteomics_project_analysis


# ---- Step 1. Quick check (about 2 minutes) ----
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 2 --nprot 30 --outtag R3smoke


# ---- Step 2, Option A: four cores (about 4 hours) ----

# A1. Mac or Linux: one command
bash 1_Code/run_5d_4cores.sh

# A2. Any computer: four Terminal tabs, one command per tab.
#     Start Tab 1 first and wait about one minute (until it prints "Proteins used: 500").

# Tab 1
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --rep_from 0 --rep_to 25

# Tab 2
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --rep_from 26 --rep_to 51

# Tab 3
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --rep_from 52 --rep_to 77

# Tab 4
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --rep_from 78 --rep_to 100


# ---- Step 2, Option B: one core (about 16 hours) ----

# B1. Normal
Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500

# B2. Mac or Linux, in the background (closing RStudio does not stop it)
nohup Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 > run_5d.log 2>&1 &

# If a run stops, run the same command again: finished replicates are kept.


# ---- Step 3. Summary (about 1 minute; run after all 101 replicate files exist) ----
Rscript 1_Code/summarize_perm_r3.r --dir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes --full 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R2/R2_IgnoreSVs/prot_multivar_ignoreSVs_emory_R2_regression_stats_diag.txt

# (If you no longer have the Round 2 main-run file, use the same command without the "--full ..." part.)
