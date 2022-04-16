rm(list=ls())
graphics.off()
cat("\014")

source("ess_comparison.R")
n_particles <- 2000 # number of particles (recommended number of particles: 10000)
iterations <- 25
alpha <- 0.6
ess_comp_3vars <- ess_comparizon(paramdim=4, n_particles=n_particles,
                                 iterations=iterations, alpha=alpha)
save(ess_comp_3vars,file="results/ess_comp_3vars.R")
