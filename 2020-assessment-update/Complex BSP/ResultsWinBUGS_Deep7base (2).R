##############################################################################
## Import results from WinBUGS txt log, create figures, stats tests
## Created by Annie Yau, annie.yau@noaa.gov, 2014.01.24
## Last modified 2014.04.25 for Main Hawaiian Islands Deep 7 Bottomfish
## Catch is in million pounds
## CPUE is in lbs/trip
##############################################################################

rm(list=ls())
library(R2WinBUGS)
library(moments)
library(lawstat)

## import data used for assessment
#DATA = read.table("X:/AYau/Bottomfish/Deep7Assess_2014/deep7_data.csv",header=T,sep=",")
#DATA = read.csv("C:\\Users\\John.Syslo\\Documents\\Deep 7 BSP\\Base\\log_sd\\Final Base\\d7_2017_data_REMLF.csv",header=T)
DATA = read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\data_NA.csv",header=T)

Survey_data = read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Survey_data.csv",header=T)

DATA=DATA[-1,]
head(DATA)

## Set working directory and read in data file
addname <- "d7_2020_base_NA_555" ##<-----------------------name of model-------Change accordingly
#src.dir = paste("X:/AYau/Bottomfish/Deep7Assess_2014/",addname,"/",sep="")  #<---Change accordingly
#src.dir<-paste('C:\\Users\\John.Syslo\\Documents\\Deep 7 BSP\\Base\\log_sd\\Final Base\\')

src.dir<-paste("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\",addname,"\\",sep="")
setwd(src.dir)

logtxt <-  bugs.log(paste(src.dir,"log.txt",sep=""))
stats <- logtxt[["stats"]]
DICresults <- logtxt[["DIC"]]

#### CHANGE THIS SECTION ACCORDING TO DATA STRUCTURE ##########################
ncountry = (ncol(DATA)-3)/2 ## number of CPUE series
count <- c("Main Hawaiian Islands Deep 7 Bottomfish") ##<------CPUE countries------ Change accordingly
yr <- DATA$Year
pyr <- 4   ## number of projected years
rnames <- rownames(stats)
cutoff <- 2004 ## year when CMLs were assigned to same individual annually, start of new CPUE series 
n1 <- sum(length(DATA$CPUE_1))-15
n2 <- n1 - length(which(yr<cutoff))


CPUE_1<-  (DATA$CPUE_1) ##<---CPUE data---Change accordingly
CPUE_2 <- (DATA$CPUE_2) ##<---CPUE data---Change accordingly

predCPUE1 <- stats[which(rnames=="PRED_CPUE[1]"):which(rnames==paste("PRED_CPUE[",n1,"]",sep="")),1] ## <--- Change accordingly
predCPUE_2 <- stats[which(rnames=="PRED_CPUE2[55]"):which(rnames==paste("PRED_CPUE2[",70,"]",sep="")),1] ## <--- Change accordingly

predCPUE <- CPUE
predCPUE[which(!is.na(predCPUE[,1])),1] <- predCPUE1
# predCPUE[which(!is.na(predCPUE[,2])),2] <- predCPUE2

slogres1 <- stats[which(rnames=="STD_LOG_RESID1[1]"):which(rnames==paste("STD_LOG_RESID1[",55,"]",sep="")),1] ## <--- Change accordingly
slogres2 <- stats[which(rnames==paste0("STD_LOG_RESID2[",(n1-n2),"]")):which(rnames==paste0("STD_LOG_RESID2[",n1+15,"]")),1] ## <--- Change accordingly
slogresAll <- rbind(as.matrix(slogres1), as.matrix(slogres2))

slogres <- CPUE
slogres[which(!is.na(slogres[,1])),1] <- slogresAll
# slogres[which(!is.na(slogres[,2])),2] <- slogres2

B <- stats[which(rnames=="B[1]"):which(rnames==paste("B[",length(yr),"]",sep="")),1] ## <--- Change accordingly
Blo <- stats[which(rnames=="B[1]"):which(rnames==paste("B[",length(yr),"]",sep="")),4] ## <--- Change accordingly
Bhi <- stats[which(rnames=="B[1]"):which(rnames==paste("B[",length(yr),"]",sep="")),6] ## <--- Change accordingly
Bstat <- stats[which(rnames=="BSTATUS[1]"):which(rnames==paste("BSTATUS[",length(yr),"]",sep="")),1] ## <--- Change accordingly 
BMSY <- stats["BMSY",1]
poflB <- stats[which(rnames=="pOFL_B[1]"):which(rnames==paste("pOFL_B[",length(yr),"]",sep="")),1] ##

H <- stats[which(rnames=="H[1]"):which(rnames==paste("H[",length(yr),"]",sep="")),1] ## <--- Change accordingly
Hlo <- stats[which(rnames=="H[1]"):which(rnames==paste("H[",length(yr),"]",sep="")),4] ## <--- Change accordingly
Hhi <- stats[which(rnames=="H[1]"):which(rnames==paste("H[",length(yr),"]",sep="")),6] ## <--- Change accordingly
Hstat <- stats[which(rnames=="HSTATUS[1]"):which(rnames==paste("HSTATUS[",length(yr),"]",sep="")),5] ## USE MEDIAN <--- Change accordingly 
HMSY <- stats["HMSY",1]
poflH <- stats[which(rnames=="pOFL_H[1]"):which(rnames==paste("pOFL_H[",length(yr),"]",sep="")),1] ##

#make big table

x<-cbind(B,Bstat,poflB,H,Hstat,poflH)
write.table(x,"clipboard-128", sep = '\t')

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

M <- stats["M",1] ## <-------------------------------------------------------------------------------------- Change accordingly
r <- stats["r",1] ## <-------------------------------------------------------------------------------------- Change accordingly
K <- stats["K",1] ## <-------------------------------------------------------------------------------------- Change accordingly
MSY <- stats["MSY",1] ## <--------------------------------------------------------------------------------- Change accordingly

#######################################################################################
#######################################



#7-20-20
###########################################JS had some difficulty running above code, and used the code below for pred vs observed CPUE###########################################################.5,4.5,3.5,2)) 
    plot( CPUE_1[1:55]~seq(1949,2003,1), type="o", col="black", lty=1, lwd=2,
         xlab="Year", ylab="Standardized CPUE (lbs/single reporting day)", 
         ylim=c(0,200),xlim=c(1948,2010),
         # xlim=c(round(yr[1],-1), round(yr[length(yr)],-1)),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,cex=1.5, xaxs="i", yaxs="i")
    lines(predCPUE1~seq(1949,2003,1), type = "o", col="black", lty = 2, lwd=1.5, pch = 0,cex=1.5)
    legend("topright",c("Observed CPUE", "Predicted CPUE"), lty=c(1,2), pch=c(1,0), lwd=c(2,1.5), cex=1.75,bty="n")
    axis(side = 2, at = seq(0, 200, by = 10), labels = FALSE, tcl = -0.3) 
    axis(side = 1, at = seq(1949, 2010, by = 1), labels = FALSE, tcl = -0.3) 
    
    

    par(mar=c(4.5,4.5,3.5,2)) 
    plot(CPUE_2[55:70]~seq(2003,2018,1), type="o", col="black", lty=1, lwd=2,
         xlab="Year", ylab="Standardized CPUE (lbs/hour)", 
         ylim=c(0,20),xlim=c(2000,2020),
         # xlim=c(round(yr[1],-1), round(yr[length(yr)],-1)),
         cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, cex=1.5, xaxs="i", yaxs="i")
    lines(predCPUE_2~seq(2003,2018), type = "o", col="black", lty = 2, lwd=1.5, pch = 0, cex=1.5)
    legend("top",c("Observed CPUE", "Predicted CPUE"), lty=c(1,2), pch=c(1,0), lwd=c(2,1.5), cex=1.75,bty="n")
    axis(side = 2, at = seq(0, 20, by = 2), labels = FALSE, tcl = -0.3) 
    axis(side = 1, at = seq(2000, 2020, by = 1), labels = FALSE, tcl = -0.3) 
    


##11-15-17 - not sure if the chunk below is useful########################
#second cpue time series
par(mar=c(4.5,4.5,3.5,2)) 
plot(yr[which(yr>=(cutoff-1))], CPUE[(yr>=(cutoff-1)),], type="o", lty=1, lwd=2,
     xlab="Year", ylab="Standardized CPUE (lbs/trip)", 
     ylim=c(0,20),xlim=c(2000,2020),cex=1.5,
     # ylim=c(0,round(1.1*max(c(CPUE[,i],predCPUE[,i]),na.rm = TRUE),digits=2))
     # xlim=c(round(yr[1],-1), round(yr[length(yr)],-1)),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xaxs = "i", yaxs = "i")
lines(yr[which(yr>=(cutoff-1))],CPUE[which(yr>=(cutoff-1)),i], type = "o", lty = 1, lwd=2, pch = 1,cex=1.5)
#lines(yr[which(yr<cutoff)],predCPUE[which(yr<cutoff),i], type = "o", lty = 2, lwd=1.5, pch = 0,cex=1.5)
lines(yr[which(yr>=(cutoff-1))],predCPUE_2, type = "o", lty = 2, lwd=1.5, pch = 0,cex=1.5)
legend('topright', c( "Observed CPUE", "Predicted CPUE"), 
       lty=c(1,1), pch=c(1,0), lwd=c(2,1.5), cex=1.5, col=c("black","black"))
#axis(side = 2, at = seq(40, 200, by = 10), labels = FALSE, tcl = -0.3) 
axis(side = 1, at = seq(2000, 2020, by = 1), labels = FALSE, tcl = -0.3) 




## Standard log residuals#################
for (i in 1:(ncountry)) {
  par(mar=c(4.5,4.5,3.5,2)) 
  barplot(slogres[,i], names.arg=yr, axis.lty = 1,
      main=count[i], xlab="Year", ylab="Standardized log residuals", xpd=F,
      ylim = c(floor(min(slogres[,i], na.rm = TRUE))-1,ceiling(max(slogres[,i], na.rm = TRUE))+1),
      cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  abline(h=0)
}



#USE THIS 2020
## Biomass
  yr<-seq(1949,2018,1)
par(mar=c(4.5,4.5,1,1)) 
plot(yr, B, type="o", col="black", lty=1, lwd=2,
     xlab="Year", ylab="Exploitable biomass (million lb)",
     ylim=c(0, 65), xlim=c(yr[1], yr[length(yr)]+1),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xaxs = "i",yaxs = "i")
lines(yr,Blo, type = "l", col="black", lty = 2, lwd=1.75) ## lower 95% CI
lines(yr,Bhi, type = "l", col="black", lty = 2, lwd=1.75) ## upper 95% CI
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Bproj, type="o", col="gray50", lty=1, lwd=2) ## projected biomass
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Bprojlo, type="l", col="gray50", lty=2, lwd=1.5 )
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Bprojhi, type="l", col="gray50", lty=2, lwd=1.5 )
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Cproj, type="o", col="black", pch=0,lty=1,lwd=1.5 ) ## projected catch
abline(h=BMSY, lty=1, col="gray60", lwd=2)
text(yr[5], BMSY-0.06*BMSY, expression(B[MSY]), col="gray60",cex=1.25)
legend('topright', c("Biomass","95% Confidence Interval"), 
       col=c("black","black"),lty=c(1,2), lwd=c(2,1.5), pch=c(1,-1), cex=1.3 )
axis(side = 2, at = seq(0, 65, by = 2), labels = FALSE, tcl = -0.3) 
axis(side = 1, at = seq(1949, 2018, by = 1), labels = FALSE, tcl = -0.3) 
abline(h=BMSY*0.844, lty=1, lwd=2)
text(yr[7], 11.5, expression("0.844*"~B[MSY]),cex=1.25)

#Add survey data points - match spm years
surv_fac<-Survey_data$Biomass_kg[1:2]*2.20462/1000000
surv_err<-2*Survey_data$SE_Biomass_kg[1:2]*2.20462/1000000
surv_hi<-surv_fac+surv_err
surv_lo<-surv_fac-surv_err
points(surv_fac~yr[69:70],col=2)
arrows(yr[69:70],surv_fac,yr[69:70],surv_hi,col=2, angle=90, length=0.1)
arrows(yr[69:70],surv_fac,yr[69:70],surv_lo,col=2, angle=90, length=0.1)

#Add survey data points - all survey years
## Biomass
par(mar=c(4.5,4.5,1,1)) 
plot(yr, B, type="o", col="black", lty=1, lwd=2,
     xlab="Year", ylab="Exploitable biomass (million lbs)",
     ylim=c(0, 60), xlim=c(yr[1], yr[length(yr)]+3),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xaxs = "i",yaxs = "i")
lines(yr,Blo, type = "l", col="black", lty = 2, lwd=1.75) ## lower 95% CI
lines(yr,Bhi, type = "l", col="black", lty = 2, lwd=1.75) ## upper 95% CI
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Bproj, type="o", col="gray50", lty=1, lwd=2) ## projected biomass
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Bprojlo, type="l", col="gray50", lty=2, lwd=1.5 )
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Bprojhi, type="l", col="gray50", lty=2, lwd=1.5 )
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Cproj, type="o", col="black", pch=0,lty=1,lwd=1.5 ) ## projected catch
abline(h=BMSY, lty=1, col="gray60", lwd=2)
text(yr[5], BMSY-0.035*BMSY, "BMSY", col="gray60",cex=1.5)
legend('topright', c("Biomass","95% Confidence Interval"), 
       col=c("black","black"),lty=c(1,2), lwd=c(2,1.5), pch=c(1,-1), cex=1.3 )
axis(side = 2, at = seq(0, 50, by = 2), labels = FALSE, tcl = -0.3) 
axis(side = 1, at = seq(1949, 2018, by = 1), labels = FALSE, tcl = -0.3) 
abline(h=BMSY*0.844, lty=1, lwd=2)
text(yr[5], 11, "0.844*BMSY",cex=1.5)

#Add survey data points - match spm years
surv_fac<-Survey_data$Biomass_kg*2.20462/1000000
surv_err<-2*Survey_data$SE_Biomass_kg*2.20462/1000000
surv_hi<-surv_fac+surv_err
surv_lo<-surv_fac-surv_err

sv_yr<-seq(2017,2020,1)
points(surv_fac~sv_yr,col=2)
arrows(sv_yr,surv_fac,sv_yr,surv_hi,col=2, angle=90, length=0.1)
arrows(sv_yr,surv_fac,sv_yr,surv_lo,col=2, angle=90, length=0.1)








##USE THIS 2020
## Harvest rate
par(mar=c(4.5,4.5,1,1)) 
plot(yr, H, type="o", col="black", lty=1, lwd=2,
     xlab="Year", ylab="Harvest rate",
     ylim = c(0,0.28), yaxp = c(0,0.3,6), xlim=c(yr[1], yr[length(yr)]),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xaxs = "i",yaxs = "i")
lines(yr,Hlo, type = "l", col="black", lty = 2, lwd=1.75) # 95% CI
lines(yr,Hhi, type = "l", col="black", lty = 2, lwd=1.75) # 95% CI
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Hproj, type="o", col="gray50", lty=1, lwd=2) ## projected biomass
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Hprojlo, type="l", col="gray50", lty=2, lwd=1.5 )
#lines((yr[length(yr)]+1):(yr[length(yr)]+pyr), Hprojhi, type="l", col="gray50", lty=2, lwd=1.5 )
abline(h=HMSY, lty=1, col="gray60", lwd=2) # add HMSY line
text(yr[15], HMSY+0.1*HMSY, expression(H[MSY]), col="gray60", cex=1.25)
legend('topleft', c("Harvest rate","95% Confidence Interval"), 
       col=c("black","black"),lty=c(1,2), lwd=c(2,1.5), pch=c(1,-1), cex=1.3 )
axis(side = 2, at = seq(0, 0.28, by = 0.01), labels = FALSE, tcl = -0.3) 
axis(side = 1, at = seq(1949, 2018, by = 1), labels = FALSE, tcl = -0.3) 


## KOBE plot
par(mar=c(4.5,5,1,1)) 
plot(Bstat, Hstat, type="o", col="black", lty=1, lwd=2,
     xlab=expression(B/B[MSY]), ylab=expression(H/H[MSY]), 
     xlim=c(0,2.0), 
     ylim=c(0,2.0),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(Bstat[1], Hstat[1],  pch=17, col="black", bg="gray50", lwd=2, cex = 2.5)
points(Bstat[length(Bstat)], Hstat[length(Hstat)], pch=17,  col="black", bg="gray50", lwd=2, cex=2.5)
# points(Bprojstat, Hprojstat, type="o", col="gray50", lty=2, lwd=1.5)
abline(h=1, lty=1, col="gray50", lwd=2)
abline(v=0.844, lty=1, col="gray50", lwd=2)
text(0.9*Bstat[1], Hstat[1], paste(yr[1]), cex=1.5)
text(0.9*Bstat[length(yr)], Hstat[length(yr)], paste(yr[length(yr)]), cex=1.5 )
axis(side = 1, at = seq(0, 2, by = 0.1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 2, by = 0.1), labels = FALSE, tcl = -0.3) 


# text(Bprojstat[pyr], 1.08*Hprojstat[pyr], paste(yr[length(yr)]+pyr), col="gray50" )
# legend('topleft', c("Modeled status","Projected status"), 
#       col=c("black","gray50"),lty=c(1,2), lwd=c(2,1.5), pch=c(1,1), cex=1.3 )


##USE THIS 2020
## KOBE plot in color################################ 11-15-17 - THIS PLOT IS AS IN REPORT#################################
par(mar=c(4.5,5,1,1)) 
plot(Bstat, Hstat, type="o", col="black",lty=1, lwd=2,
     xlab=expression(B/B[MSY]), ylab=expression(H/H[MSY]), 
     xlim=c(0,round(1.2*max(Bstat,na.rm=T), digits=1)), 
     ylim=c(0,round(1.2*max(Hstat,na.rm=T), digits=1)),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
lim <- par("usr")
rect(lim[1], 1, 0.844, lim[4], border = "firebrick1", col = "firebrick1")
rect(0.844, lim[3], lim[2], 1, border = "palegreen", col = "palegreen")
rect(lim[1], lim[3], 0.844, 1, border = "yellow", col = "yellow")
rect(0.844, 1, lim[2], lim[4], border = "orange", col = "orange")
lines(Bstat, Hstat, type="o", col="black", lty=1, lwd=2,
     xlab="B/BMSY", ylab="H/HMSY", 
     xlim=c(0,round(1.2*max(Bstat,na.rm=T), digits=1)), 
     ylim=c(0,round(1.2*max(Hstat,na.rm=T), digits=1)),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
box(col="black", lwd=2)
points(Bstat[1], Hstat[1], pch=17, col="black", bg="white", lwd=2,cex=2)
points(Bstat[length(Bstat)], Hstat[length(Hstat)], pch=17, col="black", bg="white", lwd=2,cex=2)
# points(Bprojstat, Hprojstat, type="o", col="gray50", lty=2, lwd=1.5)
abline(h=1, lty=1, col="gray50", lwd=2)
abline(v=0.844, lty=1, col="gray50", lwd=2)
text(0.9*Bstat[1], Hstat[1]-0.08, paste(yr[1]), cex=1)
text(0.9*Bstat[length(yr)], Hstat[length(yr)]-0.1, paste(yr[length(yr)]), cex=1 )
# text(Bprojstat[pyr], 1.08*Hprojstat[pyr], paste(yr[length(yr)]+pyr), col="gray50" )
# legend('topleft', c("Modeled status","Projected status"), 
#       col=c("black","gray50"),lty=c(1,2), lwd=c(2,1.5), pch=c(1,1), cex=1.3 )

## Production vs biomass
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

legend(0.4,0.4, c( "M = 0.25", "M = 1", "M = 4"), 
       lty=c(2,1,3), lwd=c(2,1.5), cex=1.5, col=c("black","black"))



###try production curves with B/K

Btest<-seq(0,1,0.1)
Ytest1 = r*Btest*(1-((Btest)^0.25)) # yield, aka surplus production  #0.25 instead of M
Ytest2 = r*Btest*(1-((Btest)^1)) # yield, aka surplus production  #1 instead of M
Ytest3 = r*Btest*(1-((Btest)^4)) # yield, aka surplus production  #4 instead of M

par(mar=c(4.5,5,1,1)) 
plot(Btest, Ytest2, type="l", col="black", lty=1, lwd=2,
     xlab="Biomass (million lbs)", ylab="Proportion of maximum surplus production", 
     xlim = c(0, 2), ylim=c(0,1.1*max(Ytest2)),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xaxs = "i",yaxs = "i")
lines(Btest, Ytest1, col="black", lty=2, lwd=2)
lines(Btest, Ytest3, col="black", lty=3, lwd=2)


#points(Btest[which.max(Ytest)], max(Ytest), pch=1, col="Black", lwd=2, cex=1.5 )
#abline(v=0.5*max(Btest[Ytest>=0]), lty=3, lwd=2, col="gray50")
# abline(v=Btest[which.max(Ytest)], lty=1, lwd=2, col="gray50")
# abline(h=max(Ytest), lty=2, lwd=2, col="gray50")
#abline(v=13.69, col="gray50")
#text(0.5*max(Btest[Ytest>=0]),0, paste0("K/2 = ",round(0.5*max(Btest[Ytest>=0]),1), " million lbs"), pos=4)
#text(Btest[which.max(Ytest)],max(Ytest), paste0("MSY = ",round(max(Ytest),1)," million lbs"), pos=4)
#text(Btest[which.max(Ytest)],0.95*max(Ytest), paste0("BMSY = ",round(Btest[which.max(Ytest)],1)," million lbs"), pos=4)
legend('topright', c(paste0("r = ",round(r,2)), paste0("K = ",round(K,2)), paste0("M = ",round(M,2))),cex=1.5)


dev.off()



### STATS ###################################################################

sink(paste(src.dir,addname,"_stats.txt",sep=""))
## Display major parameters
 writeLines("************** Parameter posteriors w/SE and 95% CI *****************
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
  for (j in 1:(ncountry)) {
    writeLines(paste0(count[i]," CPUE correlated with ",count[j], " CPUE"  ))
    print(cor(cbind(CPUE[,i],CPUE[,j]), use="pairwise.complete.obs", method="pearson"))   
    writeLines("
               ")
}}
writeLines("
            ")
  
## Analyze residuals for predicted CPUE
for (i in 1:ncountry) {
  writeLines(paste("************",count[i],"tests of residuals for fitted CPUE ************"))

  writeLines("
---------------------------------------------------------------------
LINEAR REGRESSION against time, no time trend if slope pval>0.05")
  fit1 <- lm(slogres[,i] ~ DATA$Year)
  print(summary(fit1))
  print(anova(fit1)) # anova table

  writeLines("
---------------------------------------------------------------------
BARTLETT's test, variances are homogeneous if pval>0.05, split stdlogres into 2 time series")
  test1 <- slogres1[!is.na(slogres[,i]),i]
  if (length(test1) %% 2 == 0 ) {
    list1 <- test1[1:(length(test1)/2)]
    list2 <- test1[(length(test1)/2 +1):length(test1)]}
    else {
      list1 <- test1[1:( (length(test1)+1)/2 )]
      list2 <- test1[( (length(test1)+1)/2 +1):length(test1)]
      }
  print( bartlett.test( list(list1, list2) ) )
  
  
  writeLines(paste("   
---------------------------------------------------------------------
SHAPIRO-WILK Normality test, pval>0.05 indicates normality"))
  print(shapiro.test(slogres[!is.na(slogres[,i]),i]))
  
  writeLines(paste("
---------------------------------------------------------------------
SKEWNESS ~= 0?"))
  print(skewness(slogres[,i], na.rm=T))  

  writeLines(paste("
---------------------------------------------------------------------
KURTOSIS ~= 0?"))
  print(kurtosis(slogres[,i], na.rm=T))

  writeLines("
             ")
  }


## Print out DIC results
writeLines(paste("************ DIC results **************************************
                 "))
DICresults

sink()  ## end last sink to txt

############################################################################

  