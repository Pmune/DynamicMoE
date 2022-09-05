# Libraries
library("tidyverse")
library("hrbrthemes")
library("viridis")
require(gridExtra)
hrbrthemes::import_roboto_condensed()

# Load parfaits model LPS
#
# files are saved with a tag appended to the file name to avoid overwriting files.
# if you have changed the tag please change the file_tag to the appropriate value.
file_tag <- "1" # file tag

load(paste0("results/static_poisson_dgp_results","_", file_tag, ".R")) # Static Poisson
pois_lps_diff <- static_poisson_dgp$LPS.dif

load(paste0("results/dynamic_poisson_dgp_results","_", file_tag, ".R")) # Dynamic Poisson
dynpois_lps_diff <- dynamic_poisson_dgp$LPS.dif

load(paste0("results/dynamic_poissonmix_dgp_results","_", file_tag, ".R")) # Dynamic mixture of Poisson
mix_dynpois_lps_diff <- dynamic_poissonmix_dgp$LPS.dif
iterations <- length(mix_dynpois_lps_diff)

# create the lps diff dataset
lps_diff_data <- data.frame(
  dgp=c(rep("M1",iterations), rep("M2",iterations), rep("M3",iterations)),
  lps=c( pois_lps_diff , dynpois_lps_diff, mix_dynpois_lps_diff)
)

# boxplot of the lps diff
lpsplot <- lps_diff_data %>%
  ggplot( aes(x=reorder(dgp,lps), y=lps, fill=dgp )) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.4) +
  theme_minimal() +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=11, hjust = 0.5),
    strip.background = element_rect(fill = "white", linetype = 0))+
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("LPS diff") +
  xlab(" ")

print(lpsplot)
png("plots/lps_diff_boxplot_m1_m3.png", width=600, height=300)
print(lpsplot)
dev.off()
