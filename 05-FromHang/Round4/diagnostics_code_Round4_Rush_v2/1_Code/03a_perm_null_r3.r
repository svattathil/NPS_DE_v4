###############################################################################
## 03a_perm_null_r3.r        (Round 3: permutation null, design 5d = Freedman-Lane on the protein side)
##
## Same script as Round 2 (03a_perm_null_r2.r) with these changes (it uses diag_functions_r3.r):
##   * NEW design  --design fl  (5d, the default): shuffle the part of the raw expression that the
##     nuisance variables (Batch, PMI, age, sex) do not explain, inside each Batch, keep the real NPS
##     vectors, then repeat the cleaning (NPS protected), z-scaling, GEE and tests.  See fl_functions_r3.r.
##   * the designs of Round 2 are still available: fixedx (5a), all (5b), batch (5c)
##   * the slow extras are OFF by default: no jackknife refits, no single-domain tests
##     (the HTZ / HTA / HTB / EDF / naive F / chi-square p-values come from the same fit, so they are kept)
##   * --rep_from / --rep_to : run only some of the replicates, so that several R processes
##     can work at the same time (see README.md, Option A)
##   * output goes to 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/<tag>/
## Replicate 0 is the REAL data processed the same way (for design fl: the residuals are not shuffled).
##
## Usage (from the project directory):
##   Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500                 # one core
##   Rscript 1_Code/03a_perm_null_r3.r --cohort emory --design fl --B 100 --nprot 500 --rep_from 0  --rep_to 25
## Smoke test first:  --B 2 --nprot 30
##
## Each replicate is written to its own file; if the job is stopped, running the same command again
## continues. A manifest records the settings; if the settings or the protein subset differ from the
## existing replicate files, those files are MOVED to a stale_<time>/ sub-folder, never reused.
## (Only the process that starts at replicate 0 does this check; the others stop with an error instead.)
###############################################################################
library(data.table)
setDTthreads(1)          # one thread per R process (matters when several processes run together)
library(argparse)
library(car)
library(geepack)
library(clubSandwich)

rm(list = ls())
options(stringsAsFactors = FALSE)
options(nwarnings = 100000000)
source("1_Code/project_constants.r")
source("~/comsv/Svattathil_Library/svattathil_functions.r")
source("~/comsv/Svattathil_Library/manhattan_and_qq_plot_functions.r")
source("1_Code/diag_functions_r3.r")
source("1_Code/fl_functions_r3.r")

parser <- ArgumentParser()
parser$add_argument("--cohort", type = "character", default = "emory")
parser$add_argument("--design", type = "character", default = "fl", help = "fl (5d, Freedman-Lane on the protein side), fixedx (5a), all (5b), batch (5c)")
parser$add_argument("--protect", type = "character", default = "yes", help = "yes (ignoreSVs) or no (noProtect); designs all/batch")
parser$add_argument("--na", type = "character", default = "zero", help = "zero (jaffelab) or aware; designs all/batch")
parser$add_argument("--B", type = "integer", default = 100, help = "number of permutations")
parser$add_argument("--rep_from", type = "integer", default = 0, help = "first replicate to run (0 = the real data)")
parser$add_argument("--rep_to", type = "integer", default = -1, help = "last replicate to run (-1 = B)")
parser$add_argument("--nprot", type = "integer", default = 500, help = "random subset of proteins (0 = all)")
parser$add_argument("--seed", type = "integer", default = 20261011)
parser$add_argument("--residset", type = "character", default = "ignoreSVs", help = "design fixedx only: which residual file")
parser$add_argument("--residfile", type = "character", default = NULL)
parser$add_argument("--corstr", type = "character", default = "unstructured")
parser$add_argument("--xtransform", type = "character", default = "z")
parser$add_argument("--stderr_variants", type = "character", default = "none", help = "none or comma list of jack,j1s,fij (slower; not used in Round 3)")
parser$add_argument("--perdomain", action = "store_true", help = "also run the seven single-domain tests (slower; not used in Round 3)")
parser$add_argument("--outtag", type = "character", default = NULL)
args <- parser$parse_args()
stopifnot(args$design %in% c("fl", "fixedx", "all", "batch"), args$protect %in% c("yes", "no"), args$na %in% c("zero", "aware"))
sev <- if (args$stderr_variants %in% c("none", "")) NULL else strsplit(args$stderr_variants, ",")[[1]]
extra <- list(alltests = TRUE, stderr_variants = sev, perdomain = args$perdomain)
stopifnot(args$B >= 1)
if (args$rep_to < 0) args$rep_to <- args$B
stopifnot(args$rep_from >= 0, args$rep_from <= args$rep_to, args$rep_to <= args$B)
reclean <- args$design != "fixedx"

if (is.null(args$outtag))
    args$outtag <- if (!reclean) "Perm5a_fixedx" else
        paste0(c(all = "Perm5b_all", batch = "Perm5c_batch", fl = "Perm5d_fl")[[args$design]], "_protect", args$protect,
               ifelse(args$na == "aware", "_naAware", ""))
if (args$corstr != "unstructured") args$outtag <- paste0(args$outtag, "_", args$corstr)
if (args$xtransform != "z") args$outtag <- paste0(args$outtag, "_", args$xtransform)

outdir <- paste0("2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/", args$outtag, "/")
MyMkdir(outdir)
fileid <- paste0("prot_multivar_", args$cohort, "_", args$outtag)
repfile <- function(b) sprintf("%s%s_rep%03d.txt", outdir, fileid, b)
manifest_file <- paste0(outdir, fileid, "_manifest.txt")

dat <- PrepareData(args$cohort, args$residset, args$residfile, xtransform = args$xtransform)
transformed <- dat$transformed
phenos      <- dat$phenos

### the random subset of proteins (same set and seed as Round 1 Code 5a; the same in every replicate)
set.seed(args$seed)
nsel <- if (args$nprot > 0) min(args$nprot, nrow(transformed)) else nrow(transformed)
sel  <- sort(sample.int(nrow(transformed), nsel))
selnames <- rownames(transformed)[sel]
print(paste0("Proteins used: ", nsel, " of ", nrow(transformed)))

### data for the cleaning designs
if (reclean) {
    cin <- LoadCleaningInputs(args$cohort)
    ## align the GEE phenos (PrepareData) to the order of the cleaning phenos; the sample sets must agree
    stopifnot(setequal(cin$phenos$protsample, phenos$protsample),
              all(selnames %in% rownames(cin$edata)))
    phenos_g <- phenos[match(cin$phenos$protsample, phenos$protsample), ]
    stopifnot(all(phenos_g$protsample == cin$phenos$protsample))
    ## real-data NPS vs the NPS in the two tables must agree
    for (v in npsvars.bin) stopifnot(isTRUE(all.equal(as.numeric(phenos_g[[v]]), as.numeric(cin$phenos[[v]]))))
    edata_sel <- cin$edata[selnames, , drop = FALSE]
    bt <- table(cin$phenos$Batch)
    print("Subjects per Batch:"); print(bt)
    if (args$design %in% c("fl", "batch") && any(bt < 8))
        warning("A Batch has fewer than 8 subjects; permuting inside it gives few distinct shuffles. Tell us.")
    if (args$design == "fl") {
        fl <- FLSetup(edata_sel, cin$phenos)
        print(sprintf("Design fl: nuisance model has %d columns (rank %d); %d proteins", fl$n_par_Z, fl$rank_Z, nrow(edata_sel)))
        stopifnot(fl$rank_Z == fl$n_par_Z)
    }
}

### manifest and stale-replicate guard
cur <- c(design = args$design, protect = ifelse(reclean, args$protect, NA), na = ifelse(reclean, args$na, NA),
         nprot = nsel, seed = args$seed, corstr = args$corstr, xtransform = args$xtransform,
         n_samples = nrow(phenos), residset = ifelse(reclean, NA, args$residset),
         proteins_md5 = as.character(tools::md5sum(local({ f <- tempfile(); writeLines(selnames, f); f }))),
         stderr_variants = paste(sev, collapse = ","), perdomain = extra$perdomain)
cur <- setNames(ifelse(is.na(cur) | cur == "", "-", as.character(cur)), names(cur))   # "-" = not applicable (survives fwrite/fread)
if (file.exists(manifest_file)) {
    old <- fread(manifest_file); oldv <- setNames(old$value, old$item)
    same <- all(names(cur) %in% names(oldv)) && all(mapply(function(a, b) identical(as.character(a), as.character(b)),
                                                           cur, oldv[names(cur)]))
    if (!same && args$rep_from > 0)
        stop("The manifest of this folder differs from the settings of this command. Start the block with --rep_from 0 first, or remove the folder.")
    if (!same) {
        stale <- paste0(outdir, "stale_", format(Sys.time(), "%Y%m%d_%H%M%S"), "/")
        MyMkdir(stale)
        ff <- list.files(outdir, pattern = "_rep[0-9]+\\.txt$", full.names = TRUE)
        file.rename(ff, paste0(stale, basename(ff)))
        warning("Settings differ from the existing replicate files: ", length(ff), " files moved to ", stale)
    }
}
if (args$rep_from == 0) fwrite(data.table(item = names(cur), value = cur), manifest_file, sep = "\t", quote = FALSE)

RunOneRep <- function(b) {
    if (b > 0) set.seed(args$seed + b)
    if (!reclean) {
        ph <- copy(phenos)
        if (b > 0) ph <- PermuteNPS(ph, MakeDonor(ph, "all"))
        xm <- transformed[sel, , drop = FALSE]
    } else if (args$design == "fl") {
        ## 5d: shuffle the residual of the raw expression inside Batch; the NPS vectors stay as they are
        donor <- if (b > 0) MakeDonor(cin$phenos, "batch") else seq_len(nrow(cin$phenos))
        xstar <- FLBuild(fl, donor)
        if (b == 0) stopifnot(max(abs(xstar - edata_sel), na.rm = TRUE) < 1e-8)   # replicate 0 = real expression
        ph <- phenos_g
        cl <- CleanOnce(xstar, cin$phenos, protect = args$protect, na = args$na)
        xm <- ScaleRows(cl, args$xtransform)
        keep <- intersect(selnames, rownames(xm)); xm <- xm[keep, , drop = FALSE]
    } else {
        donor <- if (b > 0) MakeDonor(cin$phenos, args$design) else seq_len(nrow(cin$phenos))
        ph <- PermuteNPS(phenos_g, donor)
        cp <- PermuteNPS(cin$phenos, donor)
        cl <- CleanOnce(edata_sel, cp, protect = args$protect, na = args$na)
        xm <- ScaleRows(cl, args$xtransform)
        keep <- intersect(selnames, rownames(xm)); xm <- xm[keep, , drop = FALSE]
    }
    rows <- vector("list", nrow(xm))
    for (j in seq_len(nrow(xm))) {
        res <- RunProteinAnalysis(x_vec = xm[j, ], sample_ids = colnames(xm), phenos = ph,
                                  npsvars = npsvars.bin, corstr = args$corstr, extra = extra)
        rows[[j]] <- ResultRow(rownames(xm)[j], res, n_na_resid = sum(is.na(xm[j, ])))
    }
    dtb <- rbindlist(lapply(rows, as.data.table), fill = TRUE)
    dtb[, rep := b]
    dtb
}

t_start <- proc.time()[["elapsed"]]
for (b in args$rep_from:args$rep_to) {
    if (file.exists(repfile(b))) {
        chk <- tryCatch(nrow(fread(repfile(b), select = "protein")), error = function(e) -1)
        if (chk > 0 && chk <= nsel) { print(paste0("replicate ", b, " exists, skipped")); next }
        file.rename(repfile(b), paste0(repfile(b), ".bad"))
    }
    dtb <- RunOneRep(b)
    fwrite(dtb, file = repfile(b), sep = "\t", quote = FALSE)
    el <- proc.time()[["elapsed"]] - t_start
    print(sprintf("replicate %d done; %.1f sec per protein; elapsed %.1f min", b,
                  sum(dtb$elapsed_sec, na.rm = TRUE) / nrow(dtb), el / 60))
}

print(paste0("Done: replicates ", args$rep_from, " to ", args$rep_to, " are in ", outdir))
nfiles <- sum(file.exists(sapply(0:args$B, repfile)))
if (nfiles == args$B + 1) {
    print("All replicates 0..B exist. Now run summarize_perm_r3.r (see README.md, Step 3).")
} else {
    print(sprintf("%d of %d replicate files exist so far; run summarize_perm_r3.r when all %d are there.", nfiles, args$B + 1, args$B + 1))
}
