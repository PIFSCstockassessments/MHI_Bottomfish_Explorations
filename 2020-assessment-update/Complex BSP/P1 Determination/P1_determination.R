##############################################################
# Results of initial proportion of carrying capacity sensitivity. 
# Used to determine the best prior mean value for P1
# Read in each PX sensitivity from 0.1 to 0.9, find the one
# which minimized MSE of the CPUE series. 
# Ideally, both series are minimized by the same P1, but
# selectively choose CPUE_1 given proximity to 1949.

# Brian Langseth, PIFSC, June 2017

# Updated June 27, 2017
#############################################################


rm(list=ls())
library(R2WinBUGS)

vals<-seq(0.1,1.0,0.1)

rmse=matrix(NA,length(vals),3,dimnames=list(c(vals),c("CPUE1","CPUE2","Survey")))

for(i in 1:length(vals)){
  addname <- paste0("P1 determinationproper_rad_NA_P1",vals[i])      #d7_2020_q3_P1
  src.dir <- paste('C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\',addname,"\\",sep="")
  setwd(src.dir)
  
  logtxt <-  bugs.log(paste(src.dir,"log.txt",sep=""))
  stats <- logtxt[["stats"]]
  DICresults <- logtxt[["DIC"]]
  
  rmse1=stats["RMSE1",1] 
  rmse2=stats["RMSE2",1] 
  survey=stats["RMSE3",1] 
  
  rmse[i,]=c(rmse1,rmse2,survey)
}

par(mar=c(5.1, 5.1, 4.1, 2.1))

plot(rmse[,"CPUE1"]+rmse[,"CPUE2"],x=vals,ylim=c(23,25),xlim=c(0,1.1),xlab="Mean prior for P1",ylab="Root mean squared error of CPUE fit",
     cex=2, type="o",col="black",pch=16,cex.lab=1.75, cex.axis=1.5, xaxs="i", yaxs="i", lwd=2)


#
plot(rmse[,"CPUE1"]+rmse[,"CPUE2"]+rmse[,"Survey"],x=vals,ylim=c(26,27.5),xlim=c(0,1.1),xlab=expression("Mean prior for "~italic("P"[italic("1")])),ylab="",
     cex=2, type="o",col="black",pch=16,cex.lab=1.75, cex.axis=1.5, xaxs="i", yaxs="i", lwd=2)
title(ylab="Root mean squared error of CPUE fit", line=2.45, cex.lab=1.75)


rmse[,"CPUE1"]+rmse[,"CPUE2"]+rmse[,"Survey"]


######Here are the original RMSEs for P1 values

P1<-seq(0.1,0.9,0.1)
RMES<-c(26.691,26.594,26.581,26.495,26.532,26.661,26.814,27.010,27.269,27.506)

points(RMES~P1)

legend(0.2,27.7, c("r = 0.1", "r = 0.2"), 
        pch=c(1,16),cex=1.3,col=c("black","red"))

