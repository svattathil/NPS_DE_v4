ROUND 4 (v2) for Rush: TERMINAL COMMANDS ONLY (explanations in README_Rush.md)
=====================================================================================
- Use the Terminal tab in RStudio. Run every command from the project folder.
- Every command is on ONE line. Copy the whole line.
- Lines that start with # are comments, not commands.
- Use Step 2 Option A OR Option B, not both.
- Do not remove --outtag (keeps the Rush results apart from the other cohorts and from earlier rounds)
  and do not remove --perdomain (it switches on the seven single-domain tests).
- NEW compared with the earlier Round 3/4 command sheets: --perdomain, Steps 3-5, and a new output folder.


# ---- Before every step: go to the project folder (the same place where you run 03a) ----
cd /path/to/31_NPS_proteomics_project_analysis


# ---- Step 0. Copy the files of 1_Code/ of this package into the project's 1_Code/ (README Step 0) ----


# ---- Step 1. Quick check (about 2 minutes or more; please note how long it takes) ----
Rscript 1_Code/03a_perm_null_r3.r --cohort rush --design fl --B 2 --nprot 30 --perdomain --outtag R4smoke_pdom_rush
Rscript 1_Code/summarize_perm_r3.r --dir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/R4smoke_pdom_rush
Rscript 1_Code/summarize_perdomain_r4.r --dir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/R4smoke_pdom_rush


# ---- Step 2, Option A: four cores ----

# A1. Mac or Linux: one command
bash 1_Code/run_5d_4cores.sh

# A2. Any computer: four Terminal tabs, one command per tab.
#     Start Tab 1 first and wait about one minute (until it prints "Proteins used: 500").

# Tab 1
Rscript 1_Code/03a_perm_null_r3.r --cohort rush --design fl --B 100 --nprot 500 --perdomain --outtag Perm5d_fl_protectyes_pdom_rush --rep_from 0 --rep_to 25

# Tab 2
Rscript 1_Code/03a_perm_null_r3.r --cohort rush --design fl --B 100 --nprot 500 --perdomain --outtag Perm5d_fl_protectyes_pdom_rush --rep_from 26 --rep_to 51

# Tab 3
Rscript 1_Code/03a_perm_null_r3.r --cohort rush --design fl --B 100 --nprot 500 --perdomain --outtag Perm5d_fl_protectyes_pdom_rush --rep_from 52 --rep_to 77

# Tab 4
Rscript 1_Code/03a_perm_null_r3.r --cohort rush --design fl --B 100 --nprot 500 --perdomain --outtag Perm5d_fl_protectyes_pdom_rush --rep_from 78 --rep_to 100


# ---- Step 2, Option B: one core ----

# B1. Normal
Rscript 1_Code/03a_perm_null_r3.r --cohort rush --design fl --B 100 --nprot 500 --perdomain --outtag Perm5d_fl_protectyes_pdom_rush

# B2. Mac or Linux, in the background (closing RStudio does not stop it)
nohup Rscript 1_Code/03a_perm_null_r3.r --cohort rush --design fl --B 100 --nprot 500 --perdomain --outtag Perm5d_fl_protectyes_pdom_rush > run_5d.log 2>&1 &

# If a run stops, run the same command again: finished replicates are kept.


# ---- Step 3. Real data, ALL proteins, with the single-domain tests (one core; may be started at the same time as Step 2 if a 5th core is free) ----
Rscript 1_Code/03a_perm_null_r3.r --cohort rush --design fl --B 1 --nprot 0 --rep_from 0 --rep_to 0 --perdomain --outtag Perm5d_fl_protectyes_pdom_rush_realall


# ---- Step 4. Phenotype table (about 1 minute) ----
Rscript 1_Code/phenotype_check_r4.r --cohort rush --outdir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes_pdom_rush


# ---- Step 5. Summaries (about 1 minute each; after all 101 replicate files of Step 2 exist, and Step 3 is finished) ----
Rscript 1_Code/summarize_perm_r3.r --dir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes_pdom_rush
Rscript 1_Code/summarize_perdomain_r4.r --dir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes_pdom_rush --realall 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/Perm5d_fl_protectyes_pdom_rush_realall

# (If Step 3 is not finished yet, run the last line without "--realall ..." and run it again later with it.)
