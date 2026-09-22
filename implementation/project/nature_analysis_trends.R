args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: nature_analysis_trends.R input.csv output.csv")

suppressPackageStartupMessages({
  library(openair)
  library(dplyr)
  library(readr)
})

input <- read_csv(args[1], show_col_types = FALSE)
input$date <- as.Date(input$date)
set.seed(20260816)

groups <- split(input, interaction(input$pair_id, input$role, drop = TRUE))
rows <- vector("list", length(groups))
i <- 0L

for (g in groups) {
  i <- i + 1L
  g <- g[order(g$date), ]
  fit <- tryCatch(
    TheilSen(
      data.frame(date = g$date, value = g$value),
      pollutant = "value",
      deseason = FALSE,
      avg.time = "month",
      slope.percent = TRUE,
      plot = FALSE,
      silent = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  r <- fit$data$res2[1, ]
  rows[[i]] <- data.frame(
    pair_id = g$pair_id[1],
    site = g$site[1],
    country = g$country[1],
    site_type = g$site_type[1],
    pollutant = g$pollutant[1],
    role = g$role[1],
    overlap_start = g$overlap_start[1],
    overlap_end = g$overlap_end[1],
    n_months = g$n_months[1],
    span_months = g$span_months[1],
    coverage_fraction = g$coverage_fraction[1],
    max_gap_months = g$max_gap_months[1],
    eligible_5y = g$eligible_5y[1],
    primary_included = g$primary_included[1],
    qc_note = g$qc_note[1],
    slope_percent = r$slope.percent,
    lower_percent = r$lower.percent,
    upper_percent = r$upper.percent,
    p_value = r$p,
    slope_absolute = r$slope,
    lower_absolute = r$lower,
    upper_absolute = r$upper,
    stringsAsFactors = FALSE
  )
}

result <- bind_rows(rows)
write_csv(result, args[2], na = "")
