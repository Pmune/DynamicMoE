rm(list=ls())
graphics.off()
cat("\014")

library(ggplot2)
library(ggthemes)
library(viridis)
library(reshape2)
load("./simulations/results/static_poison_results.R") # Static Poisson
d.Pois<- expand.grid(Discount = as.character(c(.4,.5,.6,.7,.8,.9,.95,.99))
                ,Components = c("1", "2", "3"))


d.Pois$Value<-as.vector(results$freq)
d.Pois$Model<-"M1"
KL.Po<-results$LPS.dif
LPS.Po<-results$LPS.Meanmat
load("./simulations/results/dynamic_poisson_results.R") # Dynamic Poisson
d.NB<-  expand.grid(Discount = as.character(c(.4,.5,.6,.7,.8,.9,.95,.99))
                    ,Components = c("1", "2", "3"))

d.NB$Value<-as.vector(results$freq)

d.NB$Model<-"M2"
KL.NB<-results$LPS.dif
LPS.NB<-results$LPS.Meanmat
load("./simulations/results/dynamic_poisson_mixture_results.R") # Dynamic mixture of Poisson
d.DNB<- expand.grid(Discount = as.character(c(.4,.5,.6,.7,.8,.9,.95,.99))
                    ,Components = c("1", "2", "3"))

d.DNB$Value<-as.vector(results$freq)
d.DNB$Model<-"M3"
KL.DNB<-results$LPS.dif
LPS.DNB<-results$LPS.Meanmat
d<-rbind(d.Pois,d.NB,d.DNB) # merge all frequency matrices

# LPS

M1.LPS<-expand.grid(Discount = as.character(c(.4,.5,.6,.7,.8,.9,.95,.99))
                    ,Components = c("1", "2", "3"))
M1.LPS$value<-as.vector(LPS.Po)
M1.LPS$Model<-"M1"

M2.LPS<-expand.grid(Discount = as.character(c(.4,.5,.6,.7,.8,.9,.95,.99))
                    ,Components = c("1", "2", "3"))
M2.LPS$value<-as.vector(LPS.NB)
M2.LPS$Model<-"M2"

M3.LPS<-expand.grid(Discount = as.character(c(.4,.5,.6,.7,.8,.9,.95,.99))
                    ,Components = c("1", "2", "3"))
M3.LPS$value<-as.vector(LPS.DNB)
M3.LPS$Model<-"M3"

M.LPS<-rbind(M1.LPS,M2.LPS,M3.LPS) # merge all LPS matrices

# plot of the heatmap facet

p<-ggplot(data = d, aes(x = Components , y =Discount)) +
  geom_tile(aes(fill = Value))+ labs(y = "Discount factor",x = "Number of components") +
  #scale_x_discrete(expand = c(-1, 0)) +
  #scale_y_discrete(expand = c(0, 0))+
  facet_wrap(~Model, ncol=3)+
  #scale_fill_gradient(low ="darkred", high = "yellow", guide="colorbar")+
  labs(fill="Frequency")+
  coord_equal(ratio=1/1.5)+
  theme_tufte(base_family="Helvetica")+
  theme(axis.ticks=element_blank())+
  theme(legend.title=element_text(size=10))+
  theme(legend.title.align=0)+
  theme(legend.position="right")+
  theme(legend.key.size=unit(.8, "cm"))+
  theme(legend.key.width=unit(.3, "cm"))+
  scale_fill_viridis()


p2<-ggplot(data = M.LPS, aes(x = Components , y =Discount)) +
  geom_tile(aes(fill = value))+ labs(y = "Discount factor",x = "Number of components") +
  #scale_x_discrete(expand = c(-1, 0)) +
  #scale_y_discrete(expand = c(0, 0))+
  facet_wrap(~Model, ncol=3)+
  #scale_fill_gradient(low ="darkred", high = "yellow", guide="colorbar")+
  labs(fill="Frequency")+
  coord_equal(ratio=1/1.5)+
  theme_tufte(base_family="Helvetica")+
  theme(axis.ticks=element_blank())+
  theme(legend.title=element_text(size=10))+
  theme(legend.title.align=0)+
  theme(legend.position="right")+
  theme(legend.key.size=unit(.8, "cm"))+
  theme(legend.key.width=unit(.3, "cm"))+
  scale_fill_viridis()

p


#plot of the boxplot
KL.ratios<-cbind("M1"=KL.Po,"M2"=KL.NB,"M3"=KL.DNB)
KL.ratios<-melt(KL.ratios)[,-1]

ggplot(KL.ratios, aes(x=Var2 ,y=value)) +geom_boxplot()+
         labs(x="",y="LPS difference")+scale_y_continuous(limits = c(-10, 30))+
  theme_tufte(base_family="Helvetica")+
  theme(legend.position="none")+
  theme_classic()


library(MASS)
write.matrix(LPS.DNB,file="mean_lps_dpois_mixture",sep=" ")
write.matrix(LPS.NB,file="mean_lps_dpois",sep=" ")
write.matrix(LPS.Po,file="mean_lps_pois",sep=" ")
