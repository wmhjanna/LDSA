#!/usr/bin/env Rscript
# Reproduce overall and station-class random-effects synthesis.
# Usage: Rscript scripts/02_reml_synthesis.R INPUT.csv OUTPUT_DIRECTORY [classic|HK]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Supply a site-trend CSV and an output directory.")
suppressPackageStartupMessages(library(meta))
ci_method <- if (length(args) >= 3L) args[3] else "classic"
if (!ci_method %in% c("classic", "HK")) stop("Use classic or HK.")
d <- read.csv(args[1], stringsAsFactors = FALSE)
required <- c("Site", "Station_class", "Slope_Percent", "CI_Lower_Percent", "CI_Upper_Percent")
stopifnot(all(required %in% names(d)), nrow(d) >= 2L, !anyDuplicated(d$Site))
if ("Method" %in% names(d) && length(unique(d$Method)) != 1L) stop("Pool only one bootstrap setting at a time.")
valid <- with(d, is.finite(Slope_Percent) & is.finite(CI_Lower_Percent) & is.finite(CI_Upper_Percent) &
  CI_Upper_Percent > CI_Lower_Percent & Slope_Percent >= CI_Lower_Percent & Slope_Percent <= CI_Upper_Percent)
if (!all(valid)) stop("Invalid CI or point estimate outside CI; inspect before pooling.")
# CI-width approximation used in the current synthesis, not the SD of bootstrap draws.
d$SE <- (d$CI_Upper_Percent - d$CI_Lower_Percent) / (2 * qnorm(0.975))
results <- list(); weights <- NULL
for (cl in c("Overall", "UB", "TR", "SUB", "RB")) {
  h <- if (cl == "Overall") d else d[d$Station_class == cl, ]
  if (nrow(h) < 2L) next
  m <- metagen(TE = h$Slope_Percent, seTE = h$SE, studlab = h$Site,
    common = FALSE, random = TRUE, method.tau = "REML", method.random.ci = ci_method)
  results[[cl]] <- data.frame(Station_class = cl, k = m$k, meta_CI_method = ci_method,
    effect = m$TE.random, lower = m$lower.random, upper = m$upper.random,
    p = m$pval.random, tau2 = m$tau2, I2 = m$I2, Q = m$Q)
  if (cl == "Overall") weights <- data.frame(Site = h$Site, SE = h$SE,
    weight_percent = 100 * m$w.random / sum(m$w.random))
}
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, results), file.path(args[2], "REML_summary.csv"), row.names = FALSE)
write.csv(weights, file.path(args[2], "REML_site_weights.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(args[2], "R_session_info.txt"))
