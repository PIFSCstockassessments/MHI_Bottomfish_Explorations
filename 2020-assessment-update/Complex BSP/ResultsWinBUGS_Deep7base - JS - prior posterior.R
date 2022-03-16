##############################################################################
## Import results from WinBUGS txt log, create figures, stats tests
## Created by Annie Yau, annie.yau@noaa.gov, 2014.01.24
## 
## Modified for Main Hawaiian Islands Deep 7 Bottomfish 2017 assessment
## by Brian Langseth: 
## Last updated: April 18, 2017
##############################################################################

rm(list=ls())
library(R2WinBUGS)
library(moments)
library(lawstat)

## import data used for assessment

#DATA = read.csv("C:\\Users\\John.Syslo\\Documents\\Deep 7 BSP\\Base\\log_sd\\Final Base\\Resid\\d7_2017_data_REMLF.csv",header=T)
# = read.csv("C:\\Users\\John.Syslo\\Desktop\\base model CIE\\q3 base\\CV 0.2\\d7_2017_data_REMLF.csv",header=T)
DATA = read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\data_NA.csv",header=T)

DATA=DATA[-1,]

#Need to change:
#surv=4604640/1000000*2.20462
#surv_se=(891127.6/1000000*2.20462)

Survey_data = read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Survey_data.csv",header=T)


## Set working directory, read in data file, and read in MCMC runs
addname <- "d7_2020_base_NA_555" ##<-----------------------name of model-------Change accordingly
#src.dir <- paste('D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\model code\\SPM code\\',addname,"\\",sep="") # Change accordingly
#src.dir<-paste('C:\\Users\\John.Syslo\\Documents\\Deep 7 BSP\\Base\\log_sd\\Final Base\\Resid\\')

src.dir <- paste("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\",addname,"\\",sep="")

setwd(src.dir)

logtxt <-  bugs.log(paste(src.dir,"log.txt",sep=""))
stats <- logtxt[["stats"]]
DICresults <- logtxt[["DIC"]]

coda1 = read.bugs(paste(src.dir,"coda1.txt",sep=""))
coda2 = read.bugs(paste(src.dir,"coda2.txt",sep=""))
coda3 = read.bugs(paste(src.dir,"coda3.txt",sep=""))


#######################################################################################
################ The following sections will need editing #######################

#### CHANGE THIS SECTION ACCORDING TO DATA STRUCTURE ##########################
ncountry = (ncol(DATA)-3)/2 ## number of CPUE series
count <- rep("Main Hawaiian Islands Deep 7 Bottomfish",2) ##<------CPUE countries------ Change accordingly
yr <- DATA$Year
pyr <- 4   ## number of projected years
rnames <- rownames(stats)
cutoff <- 2003 ## year when CMLs were assigned to same individual annually, start of new CPUE series
N <- nrow(DATA)
n1 <- sum(length(DATA$CPUE_1)) #Total number of records
n2 <- n1 - length(which(yr<cutoff)) #Length of later portion of timeseries
YR <- list(yr[1:(n1-n2+1)],yr[(n1-n2+1):n1])

CPUE1 <- cbind(DATA$CPUE_1) ##<---CPUE data---Change accordingly
CPUE1 <- CPUE1[!is.na(CPUE1)]
CPUE2 <- cbind(DATA$CPUE_2) ##<---CPUE data---Change accordingly
CPUE2 <- CPUE2[!is.na(CPUE2)]
CPUE <- list(CPUE1,CPUE2)

predCPUE1 <- stats[which(rnames=="PRED_CPUE[1]"):which(rnames==paste("PRED_CPUE[",n1-n2+1,"]",sep="")),1] ## <--- Change accordingly
predCPUE2 <- stats[which(rnames==paste0("PRED_CPUE2[",n1-n2+1,"]")):which(rnames==paste0("PRED_CPUE2[",n1,"]",sep="")),1] ## <--- Change accordingly
predCPUE <- list(predCPUE1,predCPUE2)

slogres1 <- stats[which(rnames=="STD_LOG_RESID1[1]"):which(rnames==paste("STD_LOG_RESID1[",n1-n2+1,"]",sep="")),1] ## <--- Change accordingly
slogres2 <- stats[which(rnames==paste0("STD_LOG_RESID2[",(n1-n2+1),"]")):which(rnames==paste0("STD_LOG_RESID2[",n1,"]")),1] ## <--- Change accordingly
slogres <- list(slogres1,slogres2)

B <- stats[which(rnames=="B[1]"):which(rnames==paste("B[",length(yr),"]",sep="")),1] ## <--- Change accordingly
Blo <- stats[which(rnames=="B[1]"):which(rnames==paste("B[",length(yr),"]",sep="")),4] ## <--- Change accordingly
Bhi <- stats[which(rnames=="B[1]"):which(rnames==paste("B[",length(yr),"]",sep="")),6] ## <--- Change accordingly
Bstat <- stats[which(rnames=="BSTATUS[1]"):which(rnames==paste("BSTATUS[",length(yr),"]",sep="")),1] ## <--- Change accordingly 
BMSY <- stats["BMSY",1]

H <- stats[which(rnames=="H[1]"):which(rnames==paste("H[",length(yr),"]",sep="")),1] ## <--- Change accordingly
Hlo <- stats[which(rnames=="H[1]"):which(rnames==paste("H[",length(yr),"]",sep="")),4] ## <--- Change accordingly
Hhi <- stats[which(rnames=="H[1]"):which(rnames==paste("H[",length(yr),"]",sep="")),6] ## <--- Change accordingly
Hstat <- stats[which(rnames=="HSTATUS[1]"):which(rnames==paste("HSTATUS[",length(yr),"]",sep="")),5] ## USE MEDIAN <--- Change accordingly 
HMSY <- stats["HMSY",1]

predS <- stats[which(rnames=="PRED_Bio2017"),1]
Slo <- stats[which(rnames=="PRED_Bio2017"),4]
Shi <- stats[which(rnames=="PRED_Bio2017"),6]

# Bproj <- stats[which(rnames=="B_PROJ[1]"):which(rnames==paste("B_PROJ[",pyr,"]",sep="")),1] ## <--- Change accordingly
# Bprojlo <- stats[which(rnames=="B_PROJ[1]"):which(rnames==paste("B_PROJ[",pyr,"]",sep="")),4] ## <--- Change accordingly
# Bprojhi <- stats[which(rnames=="B_PROJ[1]"):which(rnames==paste("B_PROJ[",pyr,"]",sep="")),6] ## <--- Change accordingly
# Bprojstat <- stats[which(rnames=="B_PROJ_STATUS[1]"):which(rnames==paste("B_PROJ_STATUS[",pyr,"]",sep="")),1] ## <--- Change accordingly
# 
# Hproj <- stats[which(rnames=="H_PROJ[1]"):which(rnames==paste("H_PROJ[",pyr,"]",sep="")),1] ## <--- Change accordingly
# Hprojlo <- stats[which(rnames=="H_PROJ[1]"):which(rnames==paste("H_PROJ[",pyr,"]",sep="")),4] ## <--- Change accordingly
# Hprojhi <- stats[which(rnames=="H_PROJ[1]"):which(rnames==paste("H_PROJ[",pyr,"]",sep="")),6] ## <--- Change accordingly
# Hprojstat <- stats[which(rnames=="H_PROJ_STATUS[1]"):which(rnames==paste("H_PROJ_STATUS[",pyr,"]",sep="")),1] ## <--- Change accordingly
# 
# Cproj <- stats[which(rnames=="CATCH_PROJ[1]"):which(rnames==paste("CATCH_PROJ[",pyr,"]",sep="")),1] ## <--- Change accordingly

M <- stats["M",1] ## <---Change accordingly
r <- stats["r",1] ## <---Change accordingly
K <- stats["K",1] ## <---Change accordingly
MSY <- stats["MSY",1] ## <---Change accordingly


plot(coda1[[1]][,"s_lambda"],trace=F,density=T,show.obs=F,xlab="slambda",main="",ylab="Posterior density") 
plot(coda2[[1]][,"s_lambda"],trace=F,density=T,show.obs=F,xlab="slambda",main="",ylab="Posterior density") 
plot(coda3[[1]][,"s_lambda"],trace=F,density=T,show.obs=F,xlab="slambda",main="",ylab="Posterior density") 


#######################################################################################
################ The following sections should not need editing #######################


### PLOTS #############################################################################

pdf(paste(src.dir,addname,"_summaryfigs.pdf",sep="")) ## save output to pdf

## Observed vs predicted CPUE
for (i in 1:(ncountry)) {
  par(mar=c(4.5,4.5,3.5,2)) 
  plot(YR[[i]], CPUE[[i]], type="o", col="black", lty=1, lwd=2,
       main=count[i], xlab="Year", ylab="CPUE (lbs/trip)", 
       ylim=c(0,round(1.2*max(c(CPUE[[i]],predCPUE[[i]]),na.rm = TRUE),digits=2)),xlim=c(min(yr),max(yr)),
       # xlim=c(round(yr[1],-1), round(yr[length(yr)],-1)),
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  lines(YR[[i]],predCPUE[[i]], type = "o", col="black", lty = 2, lwd=1.5, pch = 0)
  legend("top",c("Observed CPUE", "Predicted CPUE"), lty=c(1,2), pch=c(1,0), lwd=c(2,1.5), cex=1.3,bty="n")
}

## Standard log residuals
for (i in 1:(ncountry)) {
  par(mar=c(4.5,4.5,3.5,2)) 
  barplot(slogres[[i]], names.arg=YR[[i]], axis.lty = 1,
      main=count[i], xlab="Year", ylab="Standardized log residuals", xpd=F,
      ylim = c(floor(min(slogres[[i]], na.rm = TRUE))-1,ceiling(max(slogres[[i]], na.rm = TRUE))+1),
      cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h=0)
}

## Biomass
par(mar=c(4.5,4.5,1,1)) 
plot(yr, B, type="o", col="black", lty=1, lwd=2,
     xlab="Year", ylab="Exploitable biomass (million lbs)",
     ylim=c(0, ceiling(1.2*max(Bhi))), xlim=c(yr[1], yr[length(yr)]+2),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
lines(yr,Blo, type = "l", col="black", lty = 2, lwd=1.5) ## lower 95% CI
lines(yr,Bhi, type = "l", col="black", lty = 2, lwd=1.5) ## upper 95% CI
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Bproj, type="o", col="gray50", lty=1, lwd=2) ## projected biomass
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Bprojlo, type="l", col="gray50", lty=2, lwd=1.5 )
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Bprojhi, type="l", col="gray50", lty=2, lwd=1.5 )
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Cproj, type="o", col="black", pch=0,lty=1,lwd=1.5 ) ## projected catch
abline(h=BMSY, lty=1, col="gray60", lwd=2)
text(yr[2], BMSY-0.15*BMSY, "BMSY", col="gray60")
legend('topright', c("Biomass","95% Confidence Interval"), 
       col=c("black","black"),lty=c(1,2), lwd=c(2,1.5), pch=c(1,-1), cex=1.3, bty="n")
points(2017,predS,type="o",col="blue",pch=19)
segments(x0=2017,y0=Slo,y1=Shi,col="blue",lwd=2)
points(2017-0.5,surv)
segments(x0=2017-0.5,y0=surv+2*surv_se,y1=surv-2*surv_se,col="black",lwd=2)

#add previous assessment biomass trend
addname <- "d7_2017_baseWPSAR" ##<-----------------------name of model-------Change accordingly
#src.dir <- paste('D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\model code\\SPM code\\',addname,"\\",sep="") # Change accordingly
#src.dir<-paste('C:\\Users\\John.Syslo\\Documents\\Deep 7 BSP\\Base\\log_sd\\Final Base\\Resid\\')

src.dir <- paste("C:\\Users\\John.Syslo\\Documents\\2017 Deep 7 BSP\\Complex BSP\\",addname,"\\",sep="")

setwd(src.dir)


logtxt_old <-  bugs.log(paste(src.dir,"log.txt",sep=""))
stats_old  <- logtxt_old[["stats"]]
rnames_old <- rownames(stats_old)
yr_old<-yr[1:67]

B_old  <- stats_old [which(rnames_old=="B[1]"):which(rnames_old==paste("B[",length(yr_old),"]",sep="")),1] ## <--- Change accordingly
Blo_old  <- stats_old [which(rnames_old=="B[1]"):which(rnames_old==paste("B[",length(yr_old),"]",sep="")),4] ## <--- Change accordingly
Bhi_old  <- stats_old [which(rnames_old=="B[1]"):which(rnames_old==paste("B[",length(yr_old),"]",sep="")),6] ## <--- Change accordingly
Bstat_old  <- stats_old [which(rnames_old=="BSTATUS[1]"):which(rnames_old==paste("BSTATUS[",length(yr_old),"]",sep="")),1] ## <--- Change accordingly 
BMSY_old  <- stats_old ["BMSY",1]

lines(B_old~yr[1:67],col="red",lwd=2)
lines(Blo_old~yr[1:67],col="red",lwd=2,lty=2)
lines(Bhi_old~yr[1:67],col="red",lwd=2,lty=2)


## Harvest rate
par(mar=c(4.5,4.5,1,1)) 
plot(yr, H, type="o", col="black", lty=1, lwd=2,
     xlab="Year", ylab="Harvest rate",
     ylim = c(0,round(1.2*max(Hhi),digits=2)), yaxp = c(0,0.3,6), xlim=c(yr[1], yr[length(yr)]),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
lines(yr,Hlo, type = "l", col="black", lty = 2, lwd=1.5) # 95% CI
lines(yr,Hhi, type = "l", col="black", lty = 2, lwd=1.5) # 95% CI
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Hproj, type="o", col="gray50", lty=1, lwd=2) ## projected biomass
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Hprojlo, type="l", col="gray50", lty=2, lwd=1.5 )
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Hprojhi, type="l", col="gray50", lty=2, lwd=1.5 )
abline(h=HMSY, lty=1, col="gray60", lwd=2) # add HMSY line
text(yr[2], HMSY+0.1*HMSY, "HMSY", col="gray60")
legend('topleft', c("Harvest rate","95% Confidence Interval"), 
       col=c("black","black"),lty=c(1,2), lwd=c(2,1.5), pch=c(1,-1), cex=1.3 , bty="n" )

## KOBE plot
par(mar=c(4.5,5,1,1)) 
plot(Bstat, Hstat, type="o", col="black", lty=1, lwd=2,
     xlab=expression(B/B[MSY]), ylab=expression(H/H[MSY]), 
     xlim=c(0,round(1.2*max(Bstat,na.rm=T), digits=1)), 
     ylim=c(0,round(1.2*max(Hstat,na.rm=T), digits=1)),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(Bstat[1], Hstat[1],  pch=22, col="blue", bg="gray50", lwd=2)
points(Bstat[length(Bstat)], Hstat[length(Hstat)], pch=22,  col="green", bg="gray50", lwd=2)
# points(Bprojstat, Hprojstat, type="o", col="gray50", lty=2, lwd=1.5)
abline(h=1, lty=1, col="gray50", lwd=2)
abline(v=1, lty=1, col="gray50", lwd=2)
text(0.95*Bstat[1], 0.8*Hstat[1], paste(yr[1]), cex=1.2,col="blue")
text(1.05*Bstat[length(yr)], 1.3*Hstat[length(yr)], paste(yr[length(yr)]), cex=1.2, col="green" )
# text(Bprojstat[pyr], 1.08*Hprojstat[pyr], paste(yr[length(yr)]+pyr), col="gray50" )
# legend('topleft', c("Modeled status","Projected status"), 
#       col=c("black","gray50"),lty=c(1,2), lwd=c(2,1.5), pch=c(1,1), cex=1.3 )

## KOBE plot in color
par(mar=c(4.5,5,1,1)) 
plot(Bstat, Hstat, type="o", col="black", lty=1, lwd=2,
     xlab=expression(B/B[MSY]), ylab=expression(H/H[MSY]), 
     xlim=c(0,round(1.2*max(Bstat,na.rm=T), digits=1)), 
     ylim=c(0,round(1.2*max(Hstat,na.rm=T), digits=1)),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
lim <- par("usr")
rect(lim[1], 1, 1, lim[4], border = "red", col = "red")
rect(1, lim[3], lim[2], 1, border = "green", col = "green")
rect(lim[1], lim[3], 1, 1, border = "yellow", col = "yellow")
rect(1, 1, lim[2], lim[4], border = "yellow", col = "yellow")
lines(Bstat, Hstat, type="o", col="black", lty=1, lwd=2,
     xlab="B/BMSY", ylab="H/HMSY", 
     xlim=c(0,round(1.2*max(Bstat,na.rm=T), digits=1)), 
     ylim=c(0,round(1.2*max(Hstat,na.rm=T), digits=1)),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
box(col="black", lwd=2)
points(Bstat[1], Hstat[1], pch=22, col="black", bg="white", lwd=2)
points(Bstat[length(Bstat)], Hstat[length(Hstat)], pch=22, col="black", bg="white", lwd=2)
# points(Bprojstat, Hprojstat, type="o", col="gray50", lty=2, lwd=1.5)
abline(h=1, lty=1, col="gray50", lwd=2)
abline(v=1, lty=1, col="gray50", lwd=2)
text(0.95*Bstat[1], 0.8*Hstat[1], paste(yr[1]), cex=1.2)
text(1.05*Bstat[length(yr)], 1.3*Hstat[length(yr)], paste(yr[length(yr)]), cex=1.2 )
# text(Bprojstat[pyr], 1.08*Hprojstat[pyr], paste(yr[length(yr)]+pyr), col="gray50" )
# legend('topleft', c("Modeled status","Projected status"), 
#       col=c("black","gray50"),lty=c(1,2), lwd=c(2,1.5), pch=c(1,1), cex=1.3 )

## Production vs biomass
Btest = seq(0, 1.1*K, 0.1)
Ytest = r*Btest*(1-((Btest/K)^M)) # yield, aka surplus production
par(mar=c(4.5,5,1,1)) 
plot(Btest, Ytest, type="l", col="black", lty=1, lwd=2,
  xlab="Biomass (1000 mt)", ylab="Surplus production (million lbs)", 
  xlim = c(0, max(Btest[Ytest>=0])), ylim=c(0,1.1*max(Ytest)),
  cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(Btest[which.max(Ytest)], max(Ytest), pch=1, col="Black", lwd=2, cex=1.5 )
abline(v=0.5*max(Btest[Ytest>=0]), lty=3, lwd=2, col="gray50")
# abline(v=Btest[which.max(Ytest)], lty=1, lwd=2, col="gray50")
# abline(h=max(Ytest), lty=2, lwd=2, col="gray50")
text(0.5*max(Btest[Ytest>=0]),0, paste0("K/2 = ",round(0.5*max(Btest[Ytest>=0]),1), " million lbs"), pos=4)
text(Btest[which.max(Ytest)],max(Ytest), paste0("MSY = ",round(max(Ytest),1)," million lbs"), pos=4)
text(Btest[which.max(Ytest)],0.95*max(Ytest), paste0("BMSY = ",round(Btest[which.max(Ytest)],1)," million lbs"), pos=4)
legend('topleft', c(paste0("r = ",round(r,2)), paste0("K = ",round(K,2)), paste0("M = ",round(M,2))),bty="n") 
       
       
       20.03+0.111*20.03*(1-(20.03/27.55)^2.263)-0.639
       
       r*(1-(1/(M+1)))*K*(M+1)^(-1/M)
  

 #3/26/21
  # production vs biomass FIGURE 7       
       
       Btest = seq(0, 1.1*K, 0.5)
       Ytest1 = r*Btest*(1-((Btest/K)^0.25)) # yield, aka surplus production  #0.25 instead of M
       Ytest2 = r*Btest*(1-((Btest/K)^1)) # yield, aka surplus production  #1 instead of M
       Ytest3 = r*Btest*(1-((Btest/K)^4)) # yield, aka surplus production  #4 instead of M
       
       
       Ytest1 = r*Btest*(1-((Btest/K)^0.25))/max(r*Btest*(1-((Btest/K)^0.25))) # yield, aka surplus production  #0.25 instead of M
       Ytest2 = r*Btest*(1-((Btest/K)^1))/max(r*Btest*(1-((Btest/K)^1))) # yield, aka surplus production  #1 instead of M
       Ytest3 = r*Btest*(1-((Btest/K)^4))/max(r*Btest*(1-((Btest/K)^4))) # yield, aka surplus production  #4 instead of M
       
       
       par(mar=c(4.5,5,1,1)) 
       plot(Btest/K, Ytest2, type="l", col="black", lty=1, lwd=2,
            xlab="Biomass in proportion to carrying capacity", ylab="Proportion of maximum surplus production", 
            xlim = c(0, 1), ylim=c(0,1.1*max(Ytest2)),
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xaxs = "i",yaxs = "i")
       lines(Btest/K, Ytest1, col="black", lty=2, lwd=2)
       lines(Btest/K, Ytest3, col="black", lty=3, lwd=2)
       
       legend(0.4,0.4, c(expression(~italic("M")~"= 0.25"), expression(~italic("M")~"= 1"), expression(~italic("M")~"= 4")), 
              lty=c(2,1,3), lwd=c(2,1.5), cex=1.5, col=c("black","black"))
       
       

       
       
dev.off()

##Posteriors and priors overlap of model parameters

#Establish prior functions - Set at 10000 draws
lognormal_prior <- function(Target_VAR_Prior_Avg,CV_VAR){ 
  VAR_Prior_Precision <- 1.0/log(1+CV_VAR*CV_VAR)
  VAR_Prior_Avg  <- log(Target_VAR_Prior_Avg) - (0.5/VAR_Prior_Precision)
  VAR_Prior_sd = sqrt(1/VAR_Prior_Precision)
  VAR = rlnorm(10000,VAR_Prior_Avg,VAR_Prior_sd)
  return(VAR)
}
gamma_prior <- function(VAR_Shape,VAR_Scale){ 
  VAR = rgamma(10000,shape=VAR_Shape,rate=VAR_Scale)
  return(VAR)
}
invgamma_prior <- function(VAR_Shape,VAR_Scale){ 
  VAR_Precision = rgamma(1000,shape=VAR_Shape,rate=VAR_Scale)
  VAR = 1/VAR_Precision
  return(VAR)
}
unif_prior <- function(min,max){
  VAR = runif(100000,min,max)
  return(VAR)
}

### Now plot posteriors overlayed with their priors


##NEED TO COMBINE POSTERIORS FROM ALL CHAINS (coda1, 2, and 3)
coda_big<-rbind(as.matrix(coda1),as.matrix(coda2),as.matrix(coda3))


newdat<-cbind(coda_big[,"K"],coda_big[,"r"],coda_big[,"P[1]"],coda_big[,"M"], coda_big[,"q1"],coda_big[,"q2"],coda_big[,"rad"],coda_big[,"tau2_1"],
              coda_big[,"tau2_2"], coda_big[,"sigma2"])
pairs(newdat, labels=c("K","R", expression(P[1]), "M", expression(q[1]), expression(q[2]), "rad", expression(paste(tau[1]^2)), 
                       expression(paste(tau[2]^2)),expression(paste(sigma^2))),lower.panel=NULL)


#pairs plot with italics - 3-29-21
pairs(newdat, labels=c(expression(~italic("K")),expression(~italic("R")), expression(~italic("P"[italic("1")])),expression(~italic("M")), 
                       expression(~italic("q"[italic("1")])),expression(~italic("q"[italic("2")])),expression(~italic("rad")), expression(paste(italic(tau[1]^2))), 
                       expression(paste(italic(tau[2]^2))),expression(paste(italic(sigma^2)))),lower.panel=NULL)




correlations<-cor(newdat)
colnames(correlations)<-c("K","R", "P1", "M", "q1", "q2", "rad","tau1_2", "tau2_2", "sigma2")
rownames(correlations)<-c("K","R", "P1", "M", "q1", "q2", "rad","tau1_2", "tau2_2", "sigma2")

correlations

write.table(correlations, "clipboard", sep="\t", row.names=FALSE)
dev.new(width=10,height=8)
par(mfrow=c(4,3),mai=c(0.6,0.6,0.3,0.3))

#par(mfrow=c(2,2),mai=c(0.8,1,0.3,0.3))

#plot(coda1[[1]][,"K"],trace=F,density=T,show.obs=F,xlab="K",main="",ylab="Posterior density") 
plot(density(coda_big[,"K"]),xlab=expression(~italic("K")),main="",ylab="Density",ylim=c(0,1.1*max(density(coda_big[,"K"])$y)), 
     xlim=c(0,1.1*max(density(coda_big[,"K"])$x)),,xaxs="i",yaxs="i",cex.lab=2, cex.axis=1.5, lwd=2)
lines(density(lognormal_prior(29,0.5)),lty=2,col=2, lwd = 2)

#plot(coda1[[1]][,"r"],trace=F,density=T,show.obs=F,xlab="r",main="",ylab="Posterior density") 

plot(density(coda_big[,"r"]),xlab=expression(~italic("R")),main="",ylab="Density",ylim=c(0,1.2*max(density(coda_big[,"r"])$y)), 
     xlim=c(0,1.1*max(density(coda_big[,"r"])$x)),,xaxs="i",yaxs="i",cex.lab=2, cex.axis=1.5, lwd=2)
lines(density(lognormal_prior(0.1,0.25)),lty=2,col=2,lwd = 2)


#plot(coda1[[1]][,"P[1]"],trace=F,density=T,show.obs=F,xlab="P1",main="",ylab="Posterior density") 

plot(density(coda_big[,"P[1]"]),xlab=expression(~italic("P"[italic("1")])),main="",ylab="Density",ylim=c(0,1.2*max(density(coda_big[,"P[1]"])$y)), 
     xlim=c(0,1.1*max(density(coda_big[,"P[1]"])$x)),,xaxs="i",yaxs="i",cex.lab=2, cex.axis=1.5, lwd=2)
lines(density(lognormal_prior(0.53,0.2)),lty=2,col=2,lwd = 2)

#plot(coda1[[1]][,"M"],trace=F,density=T,show.obs=F,xlab="M",main="",ylab="Posterior density") 
plot(density(coda_big[,"M"]),xlab=expression(~italic("M")),main="",ylab="Density",ylim=c(-0.005,3*max(density(coda_big[,"M"])$y)), 
     xlim=c(0,1.1*max(density(coda_big[,"M"])$x)),,xaxs="i",yaxs="i",cex.lab=2, cex.axis=1.5, lwd=2)
lines(density(gamma_prior(0.5,0.5)),lty=2,col=2,lwd = 2) #shape and scale

#plot(coda1[[1]][,"q1"],trace=F,density=T,show.obs=F,xlab="q1",main="",ylab="Posterior density") 
plot(density(coda_big[,"q1"]),xlab=expression(~italic("q"[italic("1")])),main="",ylab="Density",ylim=c(-0.005,1.1*max(density(coda_big[,"q1"])$y)), 
     xlim=c(0,1.1*max(density(coda_big[,"q1"])$x)),,xaxs="i",yaxs="i",cex.lab=2, cex.axis=1.5, lwd=2)
lines(density(unif_prior(0.00001,100000)),lty=2,col=2,lwd = 2)
# q1col=which(colnames(coda1[[1]])=="q1[1]")
# qmeans=as.vector(NA)
# for(i in q1col:(q1col+65-1)){
#   varname=colnames(coda1[[1]])[i]
#   plot(coda1[[1]][,i],trace=F,density=T,show.obs=F,xlab=varname,main="",ylab="Posterior density") 
#   lines(density(unif_prior(0.00001,1000)),lty=2,col=2)
#   qmeans[i-q1col+1]=mean(coda1[[1]][,i])
# }
# plot(qmeans,type="l")

#plot(coda1[[1]][,"q2"],trace=F,density=T,show.obs=F,xlab="q2",main="",ylab="Posterior density") 

plot(density(coda_big[,"q2"]),xlab=expression(~italic("q"[italic("2")])),main="",ylab="Density",ylim=c(-0.005,1.1*max(density(coda_big[,"q2"])$y)), 
     xlim=c(0,1.1*max(density(coda_big[,"q2"])$x)),,xaxs="i",yaxs="i",cex.lab=2, cex.axis=1.5, lwd=2)
lines(density(unif_prior(0.00001,100000)),lty=2,col=2,lwd = 2)

plot(density(coda_big[,"rad"]),xlab=expression(~italic("rad")),main="",ylab="Density",ylim=c(0,1.2*max(density(coda_big[,"rad"])$y)), 
     xlim=c(0,1.1*max(density(coda_big[,"rad"])$x)),,xaxs="i",yaxs="i",cex.lab=2, cex.axis=1.5, lwd=2)
lines(density(lognormal_prior(27.6,0.5)),lty=2,col=2,lwd = 2)

#plot(coda1[[1]][,"sigma2"],trace=F,density=T,show.obs=F,xlab=expression(sigma^2),main="",ylab="Posterior density") 
plot(density(coda_big[,"sigma2"]),xlab=expression(paste(italic(sigma^2))),main="",ylab="Density",ylim=c(-0.05,1.1*max(density(coda_big[,"sigma2"])$y)), 
     xlim=c(0,1.1*max(density(coda_big[,"sigma2"])$x)),,xaxs="i",yaxs="i",cex.lab=2, cex.axis=1.5, lwd=2)
lines(density(invgamma_prior(0.2,0.1)),lty=2,col=2,lwd = 2)

#plot(coda1[[1]][,"tau2_1"],trace=F,density=T,show.obs=F,xlab=expression(tau[1]^2),main="",ylab="Posterior density") 
plot(density(coda_big[,"tau2_1"]),xlab=expression(paste(italic(tau[1]^2))),main="",ylab="Density",ylim=c(-0.05,1.1*max(density(coda_big[,"tau2_1"])$y)), 
     xlim=c(0,1.1*max(density(coda_big[,"tau2_1"])$x)),,xaxs="i",yaxs="i",cex.lab=2, cex.axis=1.5, lwd=2)
lines(density(invgamma_prior(0.2,1)),lty=2,col=2,lwd = 2)

#plot(density(invgamma_prior(0.2,1)),lty=2,col=2, xlim = c(0,2000000000000000000000000000))

#plot(coda1[[1]][,"tau2_2"],trace=F,density=T,show.obs=F,xlab=expression(tau[2]^2),main="",ylab="Posterior density") 
plot(density(coda_big[,"tau2_2"]),xlab=expression(paste(italic(tau[2]^2))),main="",ylab="Density",ylim=c(-0.05,1.1*max(density(coda_big[,"tau2_2"])$y)), 
     xlim=c(0,1.1*max(density(coda_big[,"tau2_2"])$x)),,xaxs="i",yaxs="i",cex.lab=2, cex.axis=1.5, lwd=2)
lines(density(invgamma_prior(0.2,1)),lty=2,col=2,lwd = 2)


## Trace plots of model parameters
##NEED TO COMBINE POSTERIORS FROM ALL CHAINS (coda1, 2, and 3)
plot(coda1[[1]][,"K"],trace=T,density=F,xlab="MCMC iteration",main="",ylab="K") 
plot(coda1[[1]][,"r"],trace=T,density=F,xlab="MCMC iteration",main="",ylab="r") 
plot(coda1[[1]][,"P[1]"],trace=T,density=F,xlab="MCMC iteration",main="",ylab="P1") 
plot(coda1[[1]][,"M"],trace=T,density=F,xlab="MCMC iteration",main="",ylab="M") 
plot(coda1[[1]][,"q1"],trace=T,density=F,xlab="MCMC iteration",main="",ylab="q1") 
plot(coda1[[1]][,"q2"],trace=T,density=F,xlab="MCMC iteration",main="",ylab="q2") 
plot(coda1[[1]][,"sigma2"],trace=T,density=F,xlab="MCMC iteration",main="",ylab="sigma2") 
plot(coda1[[1]][,"tau2_1"],trace=T,density=F,xlab="MCMC iteration",main="",ylab="tau2_1") 
plot(coda1[[1]][,"tau2_2"],trace=T,density=F,xlab="MCMC iteration",main="",ylab="tau2_2") 



#look at prior/posterior for catch in 2015

par(mfrow=c(1,1))
par(mfrow=c(3,3))

catch<-coda_big[,"B[70]"]*coda_big[,"H[70]"]
plot(density(catch),main="", xlab="Catch" ,ylab="Density",ylim=c(0,1.2*max(density(catch)$y)), 
     xlim=c(0.3,1.1*max(density(catch)$x)),xaxs="i",yaxs="i",cex.lab=1.5, cex.axis=1.5, lwd=2)
lines(density(unif_prior(0.240+0.240*0.6*1.06,0.240+0.240*1.4*1.06)),lty=2,col=2,lwd = 2)


#2010   - these match
catch<-coda_big[,"B[62]"]*coda_big[,"H[62]"]
plot(density(catch),main="", xlab="Catch" ,ylab="Density",ylim=c(0,1.2*max(density(catch)$y)), 
     xlim=c(0.0,1.1*max(density(catch)$x)),xaxs="i",yaxs="i",cex.lab=1.5, cex.axis=1.5, lwd=2)
obs_catch_2010<-0.212206
lines(density(unif_prior(obs_catch_2010+obs_catch_2010*0.6*1.06,obs_catch_2010+obs_catch_2010*1.4*1.06)),lty=2,col=2,lwd = 2)



#2000


obs_catch_2000<-0.308749

#1980



obs_catch_1980<-0.251514


#1960
catch<-coda_big[,"B[12]"]*coda_big[,"H[12]"]
plot(density(catch),main="", xlab="Catch" ,ylab="Density",ylim=c(0,1.2*max(density(catch)$y)), 
     xlim=c(0.0,1.1*max(density(catch)$x)),xaxs="i",yaxs="i",cex.lab=1.5, cex.axis=1.5, lwd=2)
obs_catch_1960<-0.163391
lines(density(unif_prior(obs_catch_1960+obs_catch_1960*0.6*1.629131,obs_catch_1960+obs_catch_1960*1.4*1.629131)),lty=2,col=2,lwd = 2)







plot(density(coda_big[,"tau2_1"]),xlab=expression(tau[1]^2),main="",ylab="Density",ylim=c(-0.05,1.1*max(density(coda_big[,"tau2_1"])$y)), 
     xlim=c(0,1.1*max(density(coda_big[,"tau2_1"])$x)),,xaxs="i",yaxs="i",cex.lab=1.5, cex.axis=1.5, lwd=2)








#look at prior/posterior for MSY - related quantities
#BMSY <- K*pow(M+1.0,(-1.0/M))
#MSY <- r*BMSY*(1.0-(1.0/(M+1.0)))
#HMSY <- min(r*(1.0-(1.0/(M+1.0))),0.999)
#PMSY <- BMSY/K
#FMSY <- -log(1-HMSY)
#CPUE_MSY <- q2*BMSY

K<-lognormal_prior(29,0.5)
M<-gamma_prior(0.5,0.5)
r<-lognormal_prior(0.1,0.25)

prior_Bmsy<-K*(M+1)^(-1/M)
prior_MSY<-r*prior_Bmsy*(1.0-(1.0/(M+1.0)))
 prior_Hmsy<-r*(1.0-(1.0/(M+1.0)))
prior_Pmsy<-prior_Bmsy/K


par(mfrow=c(2,2),mai=c(0.7,0.7,0.3,0.3))
plot(density(coda_big[,"MSY"]),xlab="MSY",main="",ylab="Density",ylim=c(0,2), 
     xlim=c(0,1.1*max(density(coda_big[,"MSY"])$x)),,xaxs="i",yaxs="i",cex.lab=1.5, cex.axis=1.5, lwd=2)
lines(density(prior_MSY),lty=2,col=2, lwd = 2)


plot(density(coda_big[,"BMSY"]),xlab=expression(B[MSY]),main="",ylab="Density",ylim=c(0,0.15), 
     xlim=c(0,1.1*max(density(coda_big[,"BMSY"])$x)),,xaxs="i",yaxs="i",cex.lab=1.5, cex.axis=1.5, lwd=2)
lines(density(prior_Bmsy),lty=2,col=2, lwd = 2)


plot(density(coda_big[,"HMSY"]),xlab=expression(H[MSY]),main="",ylab="Density",ylim=c(0,25), 
     xlim=c(0,1.1*max(density(coda_big[,"HMSY"])$x)),,xaxs="i",yaxs="i",cex.lab=1.5, cex.axis=1.5, lwd=2)
lines(density(prior_Hmsy),lty=2,col=2, lwd = 2)

plot(density(coda_big[,"PMSY"]),xlab=expression(P[MSY]),main="",ylab="Density",ylim=c(0,9), 
     xlim=c(0,1.1*max(density(coda_big[,"PMSY"])$x)),,xaxs="i",yaxs="i",cex.lab=1.5, cex.axis=1.5, lwd=2)
lines(density(prior_Pmsy),lty=2,col=2, lwd = 2)


mean(prior_MSY)
mean(coda_big[,"MSY"])

mean(coda_big[,"MSY"])/mean(prior_MSY)

mean(prior_Bmsy)
mean(coda_big[,"BMSY"])

mean(coda_big[,"BMSY"])/mean(prior_Bmsy)

mean(prior_Hmsy)
mean(coda_big[,"HMSY"])

mean(coda_big[,"HMSY"])/mean(prior_Hmsy)

mean(prior_Pmsy)
mean(coda_big[,"PMSY"])

mean(coda_big[,"PMSY"])/mean(prior_Pmsy)

sd(coda_big[,"PMSY"])


### STATS ###################################################################

sink(paste(src.dir,addname,"_stats.txt",sep=""))
## Display major parameters
 writeLines("************** Parameter posteriors w/SD and 95% CI *****************
            ")
 writeLines(paste("r = ",stats["r",1]," +/- ",stats["r",2],", [",stats["r",4],", ",stats["r",6],"]",sep=""))
 writeLines(paste("K = ",stats["K",1]," +/- ",stats["K",2],", [",stats["K",4],", ",stats["K",6],"]",sep=""))
 writeLines(paste("M = ",stats["M",1]," +/- ",stats["M",2],", [",stats["M",4],", ",stats["M",6],"]",sep=""))
 writeLines(paste("P1 = ",stats["P[1]",1]," +/- ",stats["P[1]",2],", [",stats["P[1]",4],", ",stats["P[1]",6],"]",sep=""))
 writeLines(paste("BMSY = ",stats["BMSY",1]," +/- ",stats["BMSY",2],", [",stats["BMSY",4],", ",stats["BMSY",6],"]",sep=""))
 writeLines(paste("HMSY = ",stats["HMSY",1]," +/- ",stats["HMSY",2],", [",stats["HMSY",4],", ",stats["HMSY",6],"]",sep=""))
 writeLines(paste("MSY = ",stats["MSY",1]," +/- ",stats["MSY",2],", [",stats["MSY",4],", ",stats["MSY",6],"]",sep=""))
 writeLines("
            ")

## Calculate correlations between CPUE series
writeLines("************** Pearson correlation coefficient between CPUE series **************
           ")
for (i in 1:(ncountry)) {
    writeLines("Observed CPUE correlated with predicted CPUE")
    print(cor(cbind(CPUE[[i]],predCPUE[[i]]), use="pairwise.complete.obs", method="pearson"))   
    writeLines("
               ")
}
writeLines("
            ")
  
## Analyze residuals for predicted CPUE
for (i in 1:ncountry) {
  writeLines(paste("************",count[i],"tests of residuals for fitted CPUE ************"))

  writeLines("
---------------------------------------------------------------------
LINEAR REGRESSION against time, no time trend if slope pval>0.05")
  fit1 <- lm(slogres[[i]] ~ YR[[i]])
  print(summary(fit1))
  print(anova(fit1)) # anova table

  writeLines("
---------------------------------------------------------------------
BARTLETT's test, variances are homogeneous if pval>0.05, split stdlogres into 2 time series")
  test1 <- slogres[[i]][!is.na(slogres[[i]])]
  if (length(test1) %% 2 == 0 ) {
    list1 <- test1[1:(length(test1)/2)]
    list2 <- test1[(length(test1)/2 +1):length(test1)]
  }else{
      list1 <- test1[1:( (length(test1)+1)/2 )]
      list2 <- test1[( (length(test1)+1)/2 +1):length(test1)]
      }
  print( bartlett.test( list(list1, list2) ) )
  
  
  writeLines(paste("   
---------------------------------------------------------------------
SHAPIRO-WILK Normality test, pval>0.05 indicates normality"))
  print(shapiro.test(slogres[[i]][!is.na(slogres[[i]])]))
  
  writeLines(paste("
---------------------------------------------------------------------
SKEWNESS ~= 0?"))
  print(skewness(slogres[[i]], na.rm=T))  

  writeLines(paste("
---------------------------------------------------------------------
KURTOSIS ~= 3?"))
  print(kurtosis(slogres[[i]], na.rm=T))

  writeLines("
             ")
  }


## Print out DIC results
writeLines(paste("************ DIC results **************************************
                 "))
DICresults

sink()  ## end last sink to txt

############################################################################

  
