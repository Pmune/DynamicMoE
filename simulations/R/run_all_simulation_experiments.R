
# This script runs automatically all simulation experiments in section 5.
# Each experiment can be run independently. Therefore if needed each sourced file
# can be run in its separate window to gain computation time. However, be aware
# that some experiments use the parallel computing packages `doParallel` and
# `foreach` which use a lot of memory and CPU. Running several files at once might
# slow down your computer.

# clean the workspace
rm(list=ls())
graphics.off()
cat("\014")

# install the dynamic mixture of experts package - dmoe
source("install_dependencies.R") # install all required libraries
devtools::install_github("Pmune/DynamicMoE/dmoe")

#------ Simulation experiments with dynamic parameter DGP-----------
# The experiments included in this code section are:
#       1.static Poisson DGP,
#       2.dynamic Poisson DGP,
#       3.dynamic mixture of Poisson DGP
# Note that the experiments take long time to execute.
# To reduce the running time,You may need to change the iteration value in each individual file.
source("R/static_poisson_dgp.R") # run the static Poisson simulation DGP
source("R/dynamic_poisson_dgp.R") # run the dynamic Poisson simulation DGP
source("R/dynamic_poisson_mixture_dgp.R") # run the dynamic Poisson mixture simulation DGP

# generate plots for Fig 5.1 and 5.2 in the paper.
source("plots/model_selection_heatmap_plot.R") # Fig 5.1
#source("plots/lps_diff_boxplot_dynamic_models.R") # Fig 5.2


# ---------- Autoregressive timeseries DGP experiments----------------------------
# Experiments included in this code section are:
#       1. Autoregressive timeserie DGP with static parameter.
#       2.Autoregressive timeserie DGP with dynamic parameter.

source("R/static_timeseries_dgp.R") # run the static mixture of Poisson autoregressive DGP
source("R/dynamic_timeseries_dgp.R") # run the dynamic mixture of Poisson autoregressive DGP

# generate Figure 5.3 in the paper.
source("plots/lps_diff_boxplot_autoregressive_models.R") # Fig 5.3

#----------run the effective sample size rate experiment-------------------------
source("R/ess_comparison.R") # Fig 5.4
