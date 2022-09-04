# Run all dynamic mixture of experts simulations.
#
source("R/static_poisson_dgp.R") # run the static Poisson simulation DGP
source("R/dynamic_poisson_dgp.R") # run the dynamic Poisson simulation DGP
source("R/dynamic_poisson_mixture_dgp.R") # run the dynamic Poisson mixture simulation DGP

# generate plots for Fig 5.1 and 5.2 in the paper.
source("plots/model_selection_heatmap_plots.R")
source("plots/lps_diff_boxplot_dynamic_models.R")
