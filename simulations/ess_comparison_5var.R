rm(list=ls())
graphics.off()
cat("\014")

source("ess_comparison.R")
n_particles <- 2000 # number of particles (recommended number of particles: 10000)
iterations <- 25
alpha <- 0.5
ess_comp_5vars <- ess_comparizon(paramdim=6, n_particles=n_particles,
                                 iterations=iterations, alpha=alpha)
save(ess_comp_5vars,file="results/ess_comp_5vars.R")

plot.ts(ess_comp_5vars$lb)
lines(ess_comp_5vars$ll, col="red")
