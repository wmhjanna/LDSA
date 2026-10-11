#!/usr/bin/env Rscript
# Theil-Sen trends with the openair moving-block bootstrap.
# Usage: Rscript scripts/01_theilsen_block_bootstrap.R INPUT.csv OUTPUT.csv [1999] [auto|iid|3|6|12]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Supply an input CSV and an output CSV.")
suppressPackageStartupMessages(library(openair))
if (as.character(packageVersion("openair")) != "2.19.0") {
  stop("This wrapper targets openair 2.19.0. Its internal API must be revalidated for other versions.")
}
nboot <- if (length(args) >= 3L) as.integer(args[3]) else 1999L
block_setting <- if (length(args) >= 4L) args[4] else "auto"
if (is.na(nboot) || nboot < 100L) stop("At least 100 bootstrap replicates are required.")
if (!block_setting %in% c("auto", "iid", "3", "6", "12")) stop("Unknown block setting.")
autocor <- block_setting != "iid"
fixed_block <- if (block_setting %in% c("auto", "iid")) NA_integer_ else as.integer(block_setting)

d <- read.csv(args[1], stringsAsFactors = FALSE)
if (!"Station_class" %in% names(d) && "Type" %in% names(d)) d$Station_class <- d$Type
required <- c("group_id", "Site", "Station_class", "Metric", "Phase", "date", "value")
if (!all(required %in% names(d))) stop("Missing input columns: ", paste(setdiff(required, names(d)), collapse = ", "))
d$date <- as.POSIXct(d$date, tz = "UTC")
if (anyNA(d$date) || any(!is.finite(d$value))) stop("Input dates and values must be valid; supply prepared monthly data.")
groups <- split(d, d$group_id)

fit_one <- function(g) {
  g <- g[order(g$date), ]
  if (nrow(g) < 3L || anyDuplicated(format(g$date, "%Y-%m", tz = "UTC"))) {
    stop("Each series needs at least three unique monthly observations.")
  }
  for (field in c("Site", "Station_class", "Metric", "Phase")) {
    if (length(unique(g[[field]])) != 1L) stop("Inconsistent metadata within a group.")
  }
  base_regci <- getFromNamespace("regci", "openair")
  base_sample <- getFromNamespace("samp.boot.block", "openair")
  # Restore the package functions even if fitting fails.
  on.exit({
    assignInNamespace("regci", base_regci, "openair")
    assignInNamespace("samp.boot.block", base_sample, "openair")
  }, add = TRUE)
  # nboot is not a public TheilSen argument in this version.
  reg <- base_regci
  formals(reg)$nboot <- nboot
  environment(reg) <- environment(base_regci)
  assignInNamespace("regci", reg, "openair")
  if (is.finite(fixed_block)) {
    sampler <- function(n, B, blockLength = 20) base_sample(n, B, min(fixed_block, n))
    assignInNamespace("samp.boot.block", sampler, "openair")
  }
  fit <- TheilSen(
    data.frame(date = g$date, value = g$value), pollutant = "value",
    avg.time = "month", deseason = FALSE, slope.percent = TRUE,
    autocor = autocor, plot = FALSE, silent = TRUE
  )
  r <- fit$data$res2[1, ]
  value <- function(name) if (name %in% names(r)) as.numeric(r[[name]]) else NA_real_
  method <- if (block_setting == "auto") "Auto" else if (block_setting == "iid") "IID" else paste0("Block", block_setting)
  lo <- value("lower.percent"); hi <- value("upper.percent")
  point <- value("slope.percent")
  valid <- all(is.finite(c(lo, hi, point))) && lo < hi && lo <= point && point <= hi
  classification <- if (!valid) "CI unavailable" else if (hi < 0) "Declining" else if (lo > 0) "Increasing" else "Uncertain"
  data.frame(
    group_id = g$group_id[1], Site = g$Site[1], Station_class = g$Station_class[1],
    Metric = g$Metric[1], Phase = g$Phase[1], n_months = nrow(g),
    autocor = autocor, nboot = nboot,
    block_length = if (is.finite(fixed_block)) min(fixed_block, nrow(g)) else if (autocor) round(nrow(g)^(1/3)) else 1L,
    Slope_Unit = value("slope"), CI_Lower_Unit = value("lower"), CI_Upper_Unit = value("upper"),
    Slope_Percent = value("slope.percent"), CI_Lower_Percent = lo, CI_Upper_Percent = hi,
    P_value = value("p"), Method = paste(method, nboot, sep = "_"),
    Trend_class = classification
  )
}
results <- lapply(names(groups), function(id) {
  message("Fitting ", id)
  fit_one(groups[[id]])
})
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, results), args[2], row.names = FALSE)
writeLines(capture.output(sessionInfo()), paste0(args[2], ".session.txt"))
