rm(list=ls())
graphics.off()
cat("\014")

# load the data
source("data_cleaner.R")
source("interval_predictions.R")

# load posterior output
load("outputs/selected_model_posterior.R")

# setting plot input
exp.col<- c("ChangedModule","FileComplexity")
mix.col<- c("Commits")
n_comp <- 2
pred_intervals <- c(1, 9, 20)
y.max <- 10
L<-length(time.intervals) - 1 # length of the time intervals
predictions <- interval_predictions(Data, Test.data, time.intervals, model_posterior,
                                    pred_intervals, exp.col, mix.col, n_comp, y.max)

observed_pred <- predictions$observed
predicted_pred <- predictions$predicted


## plot the predictive distributions for different intervals
windows(12,4)
par(mfrow=c(1,length(pred_intervals)), font=2,family="sans",
    font.axis=1, cex.lab=.8, cex.axis=.8, cex.main=1)
i <- 1

for (j in pred_intervals){
plot(c(0:y.max),predicted_pred[[j]]$poisson, type = "n", col = "dark red",
     lty = 1,lwd = 2, ylim = range(predicted_pred[[j]]), ylab = " ", xlab =" ", tck = -0.01,
     xaxp = c(0,y.max,y.max))

lines(c(0:y.max),predicted_pred[[j]]$poisson,type="l",col="dark red",lty=1,lwd=2,
      ylab=" ", xlab=" ", tck=-0.01, xaxp=c(0,y.max,y.max))
lines(c(0:y.max),predicted_pred[[j]]$poisson_mix,type="l",col="dark blue",lty=1,lwd=2,
      ylab=" ", xlab=" ", tck=-0.01, xaxp=c(0,y.max,y.max))
lines(c(0:y.max),observed_pred[,i], col = "dark grey",lty = 1,lwd = 3)
mtext(paste0("Interval: ", i+1), 3, cex = 0.8, line = 1, font = 2)
mtext("Number of faults", 2, cex= 0.8, line = 1, font = 2)
mtext("Density", 2, cex= 0.8, line = 1, font = 2)
legend("topright",c("Data","1 Comp", "2 Comp"),col = c("dark grey", "dark red", "dark blue"),
       lty = c(1, 1, 1),lwd = c(3, 2, 2), bty = "n",cex = .8)
i <- i + 1
}

# ----- plot of the regression coeficient posterior 2 components --------

mean_trajectory <- t(apply(as.matrix(1:L),1,function(j)
                          colMeans(model_posterior$poisson_mix$Particles[[j]])))

hpd<-t(apply(as.matrix(1:L),1,function(j)
          HPDinterval(as.mcmc(model_posterior$poisson_mix$Particles[[j]]),.95)))
upper <- hpd[,9:16]
lower <- hpd[,1:8]

sr1<-c() #screens for the first column
sr1[1]<-7
for(i in 2:3){
  sr1[i]<-sr1[i-1]+5
}
sr2<-sr1+2 # screens for the second column
#sr3<-sr2+2 # screens for the third column
xaxis<-c(1:20)
library(shape)
color<-shadepalette(n = 10, "white", "black")
colNames<-c("LogInt","CM","FC","NC","NFC","JF","CF")
# plotting the regression coefficients in the component  models
windows(12,8)
split.screen(c(5,5))
j=1
for(i in sr1){
  screen(i)
  par(mar = c(0, 0, 0.5, 0),mgp = c(3, 0.5, 0))
  plot(xaxis,mean_trajectory[,j],type="n",col="grey",ylab="",xlab="",axes=F,
       ylim=range(mean_trajectory[,j],lower[,j],upper[,j]))

  polygon(c(rev(rep(xaxis,each=2)[-c(1,2*L)]), rep(xaxis,each=2)[-c(1,2*L)]),
          c(rev(rep(lower[-L,j],each=2)),rep(upper[-L,j],each=2)), col ="light grey",border = NA)
  if(j==1){mtext("Component 1",3,cex=0.8,line=1)}
  lines(xaxis,mean_trajectory[,j],col="dark red",type="s")
  ygrid<-range(mean_trajectory[,j],lower[,j],upper[,j])
  ylablist<-round(seq(ygrid[1],ygrid[2],length.out = 4),2)
  if(j<3){abline(h=min(mean_trajectory[,j],lower[,j],upper[,j]))}
  if(j==3){
    axis(1, at=seq(0,20,by=5),labels = seq(0,20,by=5),las=1,tck=-0.03,cex.axis=.7)
    mtext("Time",1,cex=.8,line=2)
  }
  axis(2,at=ylablist,labels=ylablist,las=2,tck=-0.03,cex.axis=.7)
  mtext(colNames[j],2,cex=0.8,line=2)

  j<-j+1
}

for(i in sr2){
  screen(i)
  par(mar = c(0, 0, 0.5, 0),mgp = c(3, 0.5, 0))
  plot(xaxis,mean_trajectory[,j],type="n",col="grey",ylab="",xlab="",axes=F,
       ylim=range(mean_trajectory[,j],lower[,j],upper[,j]))

  polygon(c(rev(rep(xaxis,each=2)[-c(1,2*L)]), rep(xaxis,each=2)[-c(1,2*L)]),
          c(rev(rep(lower[-L,j],each=2)),rep(upper[-L,j],each=2)), col ="light grey",border = NA)
  if(j==4){mtext("Component 2",3,cex=0.8,line=1)}
  lines(xaxis,mean_trajectory[,j],col="dark red",type="s")
  ygrid<-range(mean_trajectory[,j],lower[,j],upper[,j])
  ylablist<-round(seq(ygrid[1],ygrid[2],length.out = 4),2)
  if(j<6){abline(h=min(mean_trajectory[,j],lower[,j],upper[,j]))}
  if(j==6){
    axis(1, at=seq(0,20,by=5),labels = seq(0,20,by=5),las=1,tck=-0.03,cex.axis=.7)
    mtext("Time",1,cex=.8,line=2)
  }
  axis(2,at=ylablist,labels=ylablist,las=2,tck=-0.03,cex.axis=.7)
  #mtext(colNames[j],2,cex=0.8,line=2)
  j<-j+1
}

# plotting the regression coefficients in the mixture weights

windows(8,4)
par(mfrow=c(1,2),mgp = c(3, 0.5, 0))
for(j in 7:8){
  plot(xaxis,mean_trajectory[,j],type="n",col="grey",ylab="",xlab="",axes=F,
       ylim=range(mean_trajectory[,j],lower[,j],upper[,j]))
  polygon(c(rev(rep(xaxis,each=2)[-c(1,2*L)]), rep(xaxis,each=2)[-c(1,2*L)]),
          c(rev(rep(lower[-L,j],each=2)),rep(upper[-L,j],each=2)), col ="light grey",border = NA)
  lines(xaxis,mean_trajectory[,j],col="dark blue",type="s")
  if (j==7){mtext("LogInt",2,cex=.8,line=2)
  }
  else{ mtext("NC",2,cex=.8,line=2)
  }
  ygrid<-range(mean_trajectory[,j],lower[,j],upper[,j])
  ylablist<-round(seq(ygrid[1],ygrid[2],length.out = 4),2)
  axis(1, at=seq(0,20,by=5), labels = seq(0,20,by=5),las=1,tck=-0.03,cex.axis=.7)
  mtext("Training periods",1,cex=.8,line=2)
  axis(2,at=ylablist,labels=ylablist,las=2,tck=-0.03,cex.axis=.7)
}

