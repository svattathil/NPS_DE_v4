###############################################################################
## summarize_perdomain_r4.r                                           (Round 4)
##
## Reads the replicate files of ONE permutation folder made by 03a_perm_null_r3.r WITH --perdomain
## (rep000 = real data, rep001... = Freedman-Lane permutations within Batch) and summarises the
## seven SINGLE-DOMAIN tests (CR2 + Satterthwaite t) in the same way as summarize_perm_r3.r summarises
## the global 7-df test. The same permutations serve both: the permuted data are drawn under "no protein-NPS
## association given Batch, PMI, age, sex", which is the null of the global test and of every domain.
##
## Files written (prefix = <dir>/<folder name>):
##   <prefix>_r4_domain_table.txt      one row per domain (+ one row "global7df", HTZ, for comparison):
##        lambda_GC of the raw p-values (real; null pooled; range over the permutations; empirical p of the real value),
##        the same for the PERMUTATION-ADJUSTED p-values (calibrated within tertiles of the Satterthwaite df, as the
##        global test is calibrated within tertiles of d2), share of p<0.05, and the Westfall-Young type p-value
##        that accounts for looking at 7 domains (max over domains of lambda_GC in each permutation)
##   <prefix>_r4_domain_by_df.txt      null calibration of each domain by tertile of its Satterthwaite df
##   <prefix>_r4_domain_by_replicate.txt  lambda_GC (raw and adjusted) and share p<0.05, per replicate and domain
##   <prefix>_r4_domain_signed.txt     direction: mean / median / sd of the t statistics, share positive (real, null)
##   <prefix>_r4_domain_corr.txt       correlation of the domain-specific t statistics across proteins (real, null)
##   <prefix>_r4_domain_cases.txt      number of subjects and cases per domain
##   <prefix>_r4_domain_pcal.txt       per protein and domain: beta, se, df, t, p, permutation-adjusted p, BH q
##   <prefix>_r4_domain_top.txt        the 10 smallest adjusted p-values of each domain
##   <prefix>_r4_domain_qq_2x4.png     Q-Q plots of the seven domains (red real, black pooled null, grey single
##                                     permutations) and the null distribution of the maximum lambda_GC
##
## Usage (from the project directory):
##   Rscript 1_Code/summarize_perdomain_r4.r --dir 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/<outtag>
##   optional:  --realall 2_Pipeline/03-Run_regressions/Prot_multivar_diag_R3/<outtag>_realall
##              (folder or file of the real-data run with ALL proteins; then the red dots, the adjusted p-values and
##               the table pcal use all proteins, the null always uses the 500 proteins of the permutation run)
###############################################################################
library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--dir", type = "character", required = TRUE, help = "folder with the *_repNNN.txt files (run with --perdomain)")
parser$add_argument("--realall", type = "character", default = NULL, help = "optional: folder or file of the real-data run with all proteins")
parser$add_argument("--prefix", type = "character", default = NULL, help = "output prefix (default: <dir>/<folder name>)")
args <- parser$parse_args()

dir <- sub("/*$", "/", args$dir)
files <- list.files(dir, pattern = "_rep[0-9]+\\.txt$", full.names = TRUE)
stopifnot(length(files) >= 2)
repno <- as.integer(sub(".*_rep([0-9]+)\\.txt$", "\\1", files))
files <- files[order(repno)]; repno <- sort(repno)
stopifnot(repno[1] == 0)
if (any(diff(repno) != 1)) warning("Replicate numbers are not consecutive: ", paste(setdiff(0:max(repno), repno), collapse = ","), " missing")
prefix <- if (is.null(args$prefix)) paste0(dir, basename(sub("/$", "", dir))) else args$prefix
all <- rbindlist(lapply(files, fread), fill = TRUE)

pcols <- grep("^pdom_", names(all), value = TRUE)
if (!length(pcols)) stop("No pdom_ columns in the replicate files: this run was made WITHOUT --perdomain. Please rerun with --perdomain.")
doms <- sub("^pdom_", "", pcols)
domlab <- c(agit = "agitation", anx = "anxiety", apa = "apathy", depd = "depression", irr = "irritability", psych = "psychosis", nite = "night-time")
lab <- function(d) ifelse(d %in% names(domlab), domlab[d], d)
for (d in doms) for (pre in c("beta_", "se_", "df_")) if (!paste0(pre, d) %in% names(all)) stop("column ", pre, d, " is missing")
if ("perdomain_error" %in% names(all)) cat("NOTE: perdomain_error was recorded for", sum(!is.na(all$perdomain_error) & all$perdomain_error != ""), "fits; first message:",
                                           head(na.omit(all$perdomain_error[all$perdomain_error != ""]), 1), "\n")
real <- all[rep == 0]; null <- all[rep > 0]
B <- length(unique(null$rep)); reps <- sort(unique(null$rep))
cat(sprintf("%d replicate files; real data: %d proteins; permutations: B = %d; domains: %s\n", length(files), nrow(real), B, paste(doms, collapse = ", ")))

realall <- NULL
if (!is.null(args$realall)) {
    ra <- args$realall
    if (dir.exists(ra)) { ff <- list.files(ra, pattern = "_rep0*\\.txt$", full.names = TRUE); ff <- ff[grepl("_rep000\\.txt$", ff)]
                          if (length(ff) != 1) stop("expected exactly one *_rep000.txt in ", ra); ra <- ff }
    realall <- fread(ra); realall <- realall[rep == 0]
    if (!all(pcols %in% names(realall))) stop("the --realall file has no pdom_ columns (it must be run with --perdomain)")
    cat(sprintf("real data, all proteins: %d\n", nrow(realall)))
}

lam <- function(p) { p <- p[!is.na(p)]; if (!length(p)) return(NA_real_); p[p == 0] <- 1e-300
                     median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, df = 1) }
shr <- function(p, a) { p <- p[!is.na(p)]; if (!length(p)) NA_real_ else mean(p < a) }
rng <- function(x) if (all(is.na(x))) c(NA_real_, NA_real_) else range(x, na.rm = TRUE)

## ---- the tests: seven domains (1 df, Satterthwaite t) and, for comparison, the global 7-df HTZ ----
tests <- lapply(doms, function(d) list(name = d, p = paste0("pdom_", d), df = paste0("df_", d)))
if (all(c("P", "df_denom") %in% names(all))) tests <- c(tests, list(list(name = "global7df", p = "P", df = "df_denom")))

## ---- permutation adjustment: calibrate within tertiles of the df, using the pooled null ----
## adjusted p = (1 + #{null p <= p in the same tertile}) / (1 + #null in the tertile); for a null replicate b the
## replicate itself is left out (leave-one-out), so that null and real data are treated alike.
CalPrep <- function(nl, tt) {
    edges <- unique(quantile(nl[[tt$df]], c(1/3, 2/3), na.rm = TRUE))
    tn <- findInterval(nl[[tt$df]], edges) + 1L; tn[is.na(nl[[tt$df]]) | is.na(nl[[tt$p]])] <- NA_integer_
    S <- lapply(seq_len(length(edges) + 1L), function(t) sort(nl[[tt$p]][which(tn == t)]))
    list(edges = edges, tn = tn, S = S)
}
CalReal <- function(dt, tt, cp) {
    tr <- findInterval(dt[[tt$df]], cp$edges) + 1L; tr[is.na(dt[[tt$df]])] <- NA_integer_
    out <- rep(NA_real_, nrow(dt)); p <- dt[[tt$p]]
    for (t in seq_along(cp$S)) { i <- which(tr == t & !is.na(p))
        if (length(i)) out[i] <- (1 + findInterval(p[i], cp$S[[t]])) / (1 + length(cp$S[[t]])) }
    out
}
CalNull <- function(nl, tt, cp) {
    out <- rep(NA_real_, nrow(nl)); p <- nl[[tt$p]]
    for (b in reps) for (t in seq_along(cp$S)) {
        i <- which(nl$rep == b & cp$tn == t & !is.na(p)); if (!length(i)) next
        Sb <- sort(p[i])
        out[i] <- (1 + findInterval(p[i], cp$S[[t]]) - findInterval(p[i], Sb)) / (1 + length(cp$S[[t]]) - length(Sb))
    }
    out
}
cps <- list()
for (tt in tests) {
    cp <- CalPrep(null, tt); cps[[tt$name]] <- cp
    null[[paste0("adj_", tt$name)]] <- CalNull(null, tt, cp)
    real[[paste0("adj_", tt$name)]] <- CalReal(real, tt, cp)
    if (!is.null(realall) && all(c(tt$p, tt$df) %in% names(realall))) realall[[paste0("adj_", tt$name)]] <- CalReal(realall, tt, cp)
}

## ---- per-replicate lambda_GC ----
byrep <- rbindlist(lapply(tests, function(tt) {
    null[, .(domain = tt$name, n = sum(!is.na(get(tt$p))), lambda = lam(get(tt$p)), lambda_adj = lam(get(paste0("adj_", tt$name))),
             share05 = shr(get(tt$p), 0.05)), by = rep]
}))
fwrite(byrep, paste0(prefix, "_r4_domain_by_replicate.txt"), sep = "\t", quote = FALSE)

## ---- main table ----
rl <- if (!is.null(realall)) realall else NULL
tab <- rbindlist(lapply(tests, function(tt) {
    nm <- tt$name; br <- byrep[domain == nm]; adj <- paste0("adj_", nm)
    lr <- lam(real[[tt$p]]); la <- lam(real[[adj]])
    r <- data.table(domain = nm, label = if (nm == "global7df") "global 7-df (HTZ)" else lab(nm),
        n_real = sum(!is.na(real[[tt$p]])), n_null = sum(!is.na(null[[tt$p]])),
        lambda_real = lr, lambda_null = lam(null[[tt$p]]),
        lambda_null_rep_min = rng(br$lambda)[1], lambda_null_rep_max = rng(br$lambda)[2],
        n_rep_ge_real = sum(br$lambda >= lr, na.rm = TRUE), emp_p_lambda = (1 + sum(br$lambda >= lr, na.rm = TRUE)) / (1 + B),
        lambda_adj_real = la, lambda_adj_null = lam(null[[adj]]),
        lambda_adj_null_rep_min = rng(br$lambda_adj)[1], lambda_adj_null_rep_max = rng(br$lambda_adj)[2],
        emp_p_lambda_adj = (1 + sum(br$lambda_adj >= la, na.rm = TRUE)) / (1 + B),
        share05_real = shr(real[[tt$p]], 0.05), share05_null = shr(null[[tt$p]], 0.05),
        share05_null_rep_min = rng(br$share05)[1], share05_null_rep_max = rng(br$share05)[2])
    if (!is.null(rl) && tt$p %in% names(rl))
        r[, `:=`(n_real_all = sum(!is.na(rl[[tt$p]])), lambda_real_all = lam(rl[[tt$p]]), lambda_adj_real_all = lam(rl[[adj]]),
                 share05_real_all = shr(rl[[tt$p]], 0.05))]
    r
}), fill = TRUE)
## Westfall-Young type p-value: is the real lambda of this domain larger than the largest lambda of the 7 domains in a permutation?
dm <- byrep[domain %in% doms]
mx <- dm[, .(mx = max(lambda, na.rm = TRUE), mxa = max(lambda_adj, na.rm = TRUE)), by = rep]
tab[domain %in% doms, wy_p_lambda := sapply(lambda_real, function(l) (1 + sum(mx$mx >= l)) / (1 + B))]
tab[domain %in% doms, wy_p_lambda_adj := sapply(lambda_adj_real, function(l) (1 + sum(mx$mxa >= l)) / (1 + B))]
fwrite(tab, paste0(prefix, "_r4_domain_table.txt"), sep = "\t", quote = FALSE)
print(tab[, .(label, lambda_real = round(lambda_real, 2), lambda_null = round(lambda_null, 2), emp_p = round(emp_p_lambda, 3),
              lambda_adj_real = round(lambda_adj_real, 2), lambda_adj_null = round(lambda_adj_null, 2), wy_p = round(wy_p_lambda, 3))])

## ---- by df tertile ----
bydf <- rbindlist(lapply(tests, function(tt) {
    cp <- cps[[tt$name]]; nl <- copy(null); nl[, tn := cp$tn]
    nl <- nl[!is.na(tn)]
    a <- nl[, .(data = "null", n = sum(!is.na(get(tt$p))), median_df = median(get(tt$df), na.rm = TRUE), share05 = shr(get(tt$p), 0.05),
                share001 = shr(get(tt$p), 0.001), lambda = lam(get(tt$p))), by = tn]
    tr <- findInterval(real[[tt$df]], cp$edges) + 1L
    b <- real[!is.na(tr), .(data = "real", n = sum(!is.na(get(tt$p))), median_df = median(get(tt$df), na.rm = TRUE), share05 = shr(get(tt$p), 0.05),
                           share001 = shr(get(tt$p), 0.001), lambda = lam(get(tt$p))), by = .(tn = tr[!is.na(tr)])]
    cbind(domain = tt$name, rbind(a, b))
}), fill = TRUE)
setnames(bydf, "tn", "df_tertile")
fwrite(bydf, paste0(prefix, "_r4_domain_by_df.txt"), sep = "\t", quote = FALSE)

## ---- direction of the effects ----
tstat <- function(dt, d) dt[[paste0("beta_", d)]] / dt[[paste0("se_", d)]]
sg <- rbindlist(lapply(doms, function(d) {
    f <- function(dt, w) { t <- tstat(dt, d); p <- dt[[paste0("pdom_", d)]]; ok <- !is.na(t)
        data.table(domain = d, data = w, n = sum(ok), mean_t = mean(t[ok]), median_t = median(t[ok]), sd_t = sd(t[ok]),
                   frac_pos = mean(t[ok] > 0), n_p05 = sum(p < 0.05, na.rm = TRUE), frac_pos_among_p05 = mean(t[which(p < 0.05)] > 0, na.rm = TRUE)) }
    r <- list(f(real, "real_500"), f(null, "null_pooled"))
    if (!is.null(realall)) r <- c(r, list(f(realall, "real_all")))
    x <- rbindlist(r); x[data == "null_pooled", mean_t_rep_sd := null[, mean(tstat(.SD, d), na.rm = TRUE), by = rep][, sd(V1)]]; x
}), fill = TRUE)
fwrite(sg, paste0(prefix, "_r4_domain_signed.txt"), sep = "\t", quote = FALSE)

## ---- correlation of the t statistics between domains, across proteins ----
tm <- function(dt) sapply(doms, function(d) tstat(dt, d))
cr_real <- cor(tm(real), use = "pairwise.complete.obs")
cr_null <- lapply(reps, function(b) cor(tm(null[rep == b]), use = "pairwise.complete.obs"))
cn_mean <- Reduce(`+`, cr_null) / length(cr_null)
cn_sd <- apply(simplify2array(cr_null), c(1, 2), sd, na.rm = TRUE)
cc <- CJ(d1 = seq_along(doms), d2 = seq_along(doms))[d1 < d2]
cc[, `:=`(domain1 = doms[d1], domain2 = doms[d2], cor_real = cr_real[cbind(d1, d2)], cor_null_mean = cn_mean[cbind(d1, d2)], cor_null_sd = cn_sd[cbind(d1, d2)])]
if (!is.null(realall)) { cra <- cor(tm(realall), use = "pairwise.complete.obs"); cc[, cor_real_all := cra[cbind(d1, d2)]] }
fwrite(cc[, !c("d1", "d2")], paste0(prefix, "_r4_domain_corr.txt"), sep = "\t", quote = FALSE)

## ---- case counts (from the text column domain_counts of the real data) ----
if ("domain_counts" %in% names(real)) {
    m <- regmatches(real$domain_counts, gregexpr("[A-Za-z0-9_]+:n=[0-9]+,cases=[0-9]+", real$domain_counts))
    ct <- rbindlist(lapply(m, function(z) if (length(z)) data.table(domain = sub(":.*", "", z), n = as.integer(sub(".*:n=([0-9]+),.*", "\\1", z)),
                                                                   cases = as.integer(sub(".*cases=", "", z)))))
    if (nrow(ct)) fwrite(ct[, .(median_n = median(n), median_cases = median(cases), min_cases = min(cases), max_cases = max(cases), proteins_lt5_cases = sum(cases < 5)), by = domain],
                         paste0(prefix, "_r4_domain_cases.txt"), sep = "\t", quote = FALSE)
}

## ---- per protein and domain: adjusted p-values and BH q ----
base <- if (!is.null(realall)) realall else real
pc <- rbindlist(lapply(doms, function(d) {
    x <- data.table(protein = base$protein, domain = d, beta = base[[paste0("beta_", d)]], se = base[[paste0("se_", d)]],
                    df = base[[paste0("df_", d)]], p = base[[paste0("pdom_", d)]], p_perm = base[[paste0("adj_", d)]])
    x[, t := beta / se]; x[, q_BH := p.adjust(p_perm, "BH")]; x
}))
setcolorder(pc, c("protein", "domain", "beta", "se", "df", "t", "p", "p_perm", "q_BH"))
fwrite(pc, paste0(prefix, "_r4_domain_pcal.txt"), sep = "\t", quote = FALSE)
top <- pc[!is.na(p_perm)][order(p_perm), head(.SD, 10), by = domain]
fwrite(top, paste0(prefix, "_r4_domain_top.txt"), sep = "\t", quote = FALSE)

## ---- figure ----
obs_qq <- function(p) { p <- sort(p[!is.na(p)]); p[p == 0] <- 1e-300; m <- length(p); list(e = -log10((seq_len(m) - 0.5) / (m + 1)), o = -log10(p)) }
png(paste0(prefix, "_r4_domain_qq_2x4.png"), height = 8.4, width = 16, units = "in", res = 200)
par(mfrow = c(2, 4), mar = c(4.2, 4.2, 3, 1), pty = "s")
for (d in doms) {
    v <- paste0("pdom_", d); rr <- if (!is.null(realall)) realall[[v]] else real[[v]]
    cr <- obs_qq(rr); cn <- obs_qq(null[[v]]); mx2 <- max(cr$e, cr$o, cn$e, cn$o, 1)
    plot(NA, xlim = c(0, mx2), ylim = c(0, mx2), xlab = expression(Expected~-log[10](p)), ylab = expression(Observed~-log[10](p)),
         main = sprintf("%s: lambda real %.2f, null %.2f", lab(d), lam(real[[v]]), lam(null[[v]])), cex.main = 1.0)
    abline(0, 1, col = "grey40")
    for (b in reps) { cb <- obs_qq(null[rep == b][[v]]); points(cb$e, cb$o, col = "grey75", pch = 16, cex = 0.25) }
    points(cn$e, cn$o, col = "black", pch = 16, cex = 0.35); points(cr$e, cr$o, col = "red", pch = 16, cex = 0.45)
    if (d == doms[1]) legend("topleft", bty = "n", pch = 16, col = c("red", "black", "grey75"), cex = 0.9,
                             legend = c(ifelse(is.null(realall), "real data (permutation-run proteins)", "real data (all proteins)"), "pooled permutations", "single permutations"))
}
for (k in seq_len(8 - length(doms) - 1)) plot.new()
lr_all <- tab[domain %in% doms, lambda_real]
hist(mx$mx, breaks = 20, col = "grey80", border = "white", main = "Largest lambda over the domains", xlab = "max over domains of lambda_GC",
     xlim = range(c(mx$mx, lr_all), na.rm = TRUE) * c(0.95, 1.05))
abline(v = lr_all, col = "red", lty = 2); abline(v = max(lr_all), col = "red", lwd = 2)
dev.off()
cat("Written: ", paste0(prefix, c("_r4_domain_table.txt", "_r4_domain_by_df.txt", "_r4_domain_by_replicate.txt", "_r4_domain_signed.txt",
    "_r4_domain_corr.txt", "_r4_domain_cases.txt", "_r4_domain_pcal.txt", "_r4_domain_top.txt", "_r4_domain_qq_2x4.png"), collapse = "\n         "), "\n", sep = "")
