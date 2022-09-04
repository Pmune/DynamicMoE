# Run all autoregressive timeseries mixture of experts simulations.
rm(list=ls())
graphics.off()
cat("\014")

source("R/static_timeseries_dgp.R")
source("R/dynamic_timeseries_dgp.R")

# generate Figure 5.3 in the paper.

source("plots/lps_diff_boxplot_autoregressive_models.R")
