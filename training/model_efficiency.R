rm(list=ls())
graphics.off()
cat("\014")


library(dmoe)

# load the data
source("data_cleaner.R")


# ------set the model parameters  ---------

set.seed(500)
y.max<-8
exp.col<-c("ChangedModule","FileComplexity")
mix.col<-c("Commits")
L<-length(time.intervals)-1 # length of the time intervals
N.particles<-1000 # number of particles (recommended number of particles: 10000)
alpha<-0.5
R<-2
output_efficiency <- list()
# dynamic model with 2 components
for (iter in seq_len(100)){
system.time(output<-FAPF(Data, time.intervals, N.particles, mix.col, exp.col,
                         n_comp=2, time_ind=1, response_ind=2, alpha, R, T, h=1))

  output_efficiency [[iter]]<- output

}


save(output_efficiency,file="outputs/posterior_for_model_efficiency.R")

