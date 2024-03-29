# plot the heatmap of the model selection frequency. Figure 5.1 in the paper.
rm(list=ls())

library(ggplot2)
library(ggthemes)
library(viridis)
library(reshape2)

# files are saved with a tag appended to the file name to avoid overwriting files.
# if you have changed the tag please change the file_tag to the appropriate value.
file_tag <- "1" # file tag

load(paste0("results/static_poisson_dgp_results","_", file_tag, ".R")) # Static Poisson
static_pois<- expand.grid(Discount = as.character(c(.3,.4,.5,.6,.7,.8,.9,.95,.99))
                ,Components = c("1", "2", "3"))


static_pois$Value<-as.vector(static_poisson_dgp$freq)
static_pois$Model<-"M1"


load(paste0("results/dynamic_poisson_dgp_results","_", file_tag, ".R")) # Dynamic Poisson
dyn_pois<-  expand.grid(Discount = as.character(c(.3,.4,.5,.6,.7,.8,.9,.95,.99))
                    ,Components = c("1", "2", "3"))
dyn_pois$Value<-as.vector(dynamic_poisson_dgp$freq)
dyn_pois$Model<-"M2"


load(paste0("results/dynamic_poissonmix_dgp_results","_", file_tag, ".R")) # Dynamic mixture of Poisson
dyn_pois_mix<- expand.grid(Discount = as.character(c(.3,.4,.5,.6,.7,.8,.9,.95,.99))
                    ,Components = c("1", "2", "3"))
dyn_pois_mix$Value<-as.vector(dynamic_poissonmix_dgp$freq)
dyn_pois_mix$Model<-"M3"



selection_freqs <- rbind(static_pois, dyn_pois, dyn_pois_mix) # merge all frequency matrices


# plot of the model selection frequency heatmap

freq_plot<-ggplot(data = selection_freqs, aes(x = Components , y =Discount)) +
  geom_tile(aes(fill = Value))+ labs(y = "Discount factor",x = "Number of components") +
  facet_wrap(~Model, ncol=3)+
  labs(fill="Frequency")+
  theme_minimal()+
  theme(axis.ticks=element_blank(),
        plot.title = element_text(size=11, hjust = 0.5),
        strip.background = element_rect(fill = "white", linetype = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right")+

scale_fill_viridis(alpha=0.8)

print(freq_plot)
png("plots/model_selection_frequency_heatmap.png", width=600, height=300)
print(freq_plot)
dev.off()

