
# This script runs automatically all simulation experiments in section 5.
# Each experiment can be run independently. Therefore if needed each sourced file
# can be run in its separate window to gain computation time. However, be aware
# that some experiments use the parallel computing packages `doParallel` and
# `foreach` which use a lot of memory and CPU. Running several files at once might
# slow down your computer.

source("install_dependencies.R") # install all required libraries

#------ Simulation experiments with dynamic parameter DGP-----------
# The experiments included in this code section are:
#       1.static Poisson DGP,
#       2.dynamic Poisson DGP,
#       3.dynamic mixture of Poisson DGP
# The experiments take long time to execute.To reduce the running time,
# I have reduced the iterations in each experiment to 10 iterations.
# With 10 iterations each experiment bellow takes roughly 2 hours.
source("R/static_poisson_dgp.R") # run the static Poisson simulation DGP
source("R/dynamic_poisson_dgp.R") # run the dynamic Poisson simulation DGP
source("R/dynamic_poisson_mixture_dgp.R") # run the dynamic Poisson mixture simulation DGP

# generate plots for Fig 5.1 and 5.2 in the paper.
source("plots/model_selection_heatmap_plot.R")
source("plots/lps_diff_boxplot_dynamic_models.R")


# ---------- Autoregressive timeseries DGP experiments----------------------------
# Experiments included in this code section are:
#       1. Autoregressive timeserie DGP with static parameter.
#       2.  Autoregressive timeserie DGP with dynamic parameter.

source("R/static_timeseries_dgp.R")
source("R/dynamic_timeseries_dgp.R")

# generate Figure 5.3 in the paper.
source("plots/lps_diff_boxplot_autoregressive_models.R")

#----------run the effective sample size rate experiment-------------------------
source("R/ess_comparison.R")
