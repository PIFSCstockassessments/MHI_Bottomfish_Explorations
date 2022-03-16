#Model diagnostics as separate file from results



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

slogres3 <- stats[which(rnames==paste0("STD_LOG_RESID3[",69,"]")):which(rnames==paste0("STD_LOG_RESID3[",72,"]")),1] ## <--- Change accordingly


##################################
#plots of CPUE fits###############
##################################

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



#Survey CPUE


predSurv<- stats[which(rnames=="PRED_Bio[69]"):which(rnames==paste("PRED_Bio[",72,"]",sep="")),1]*1000000 
Survey_B<-Survey_data$Biomass_kg*2.20462/(25892*104.4653) #removed factoring out millions


par(mar=c(4.5,4.5,3.5,2)) 
plot(Survey_B~seq(2017,2020,1), type="o", col="black", lty=1, lwd=2,
     xlab="Year", ylab="Survey CPUE (lb/camera-swept area)", 
     ylim=c(2,6),xlim=c(2016,2021),
     # xlim=c(round(yr[1],-1), round(yr[length(yr)],-1)),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, cex=1.5, xaxs="i", yaxs="i")
lines(predSurv~seq(2017,2020,1), type = "o", col="black", lty = 2, lwd=1.5, pch = 0, cex=1.5)
legend("topright",c("Observed Survey", "Predicted Survey"), lty=c(1,2), pch=c(1,0), lwd=c(2,1.5), cex=1.75,bty="n")
axis(side = 2, at = seq(0, 20, by = 2), labels = FALSE, tcl = -0.3) 
axis(side = 1, at = seq(2016, 2021, by = 1), labels = FALSE, tcl = -0.3) 





################################
#Residual diagnostics###########
################################
#trend:

#summary(lm(slogres1~yr))
# summary(lm(slogres1~c(seq(1,6,1),c(1,9,10),c(12,13),c(seq(15,55,1)))))  #removing large values - 2017
summary(lm(slogres1~seq(1,55,1)))  #PASS -no n eed to remove large values 2020 -not significant


#normality:
shapiro.test(slogres1)   #2020 this is significantly different from normal with all values
shapiro.test(slogres1[-c(7,11,14)]) #PASSES without 1955,1959,1962
# constant variance:
library(car)  
ncvTest(lm(slogres1~seq(1,55,1)))

#BARTLETT's test, variances are homogeneous if pval>0.05, split stdlogres into 2 time series")
test1 <- slogres1[!is.na(slogres1)]
if (length(test1) %% 2 == 0 ) {
  list1 <- test1[1:(length(test1)/2)]
  list2 <- test1[(length(test1)/2 +1):length(test1)]
}else{
  list1 <- test1[1:( (length(test1)+1)/2 )]
  list2 <- test1[( (length(test1)+1)/2 +1):length(test1)]
}
print( bartlett.test( list(list1, list2) ) ) #FAIL



yr<-seq(1949,2003,1)
barplot(slogres1, axis.lty = 1,
        names.arg=yr, xlab="Year", ylab="Standardized log residuals", xpd=F,
        ylim = c(floor(min(slogres1, na.rm = TRUE))-1,ceiling(max(slogres1, na.rm = TRUE))+1),
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xaxs = "i",yaxs = "i")



abline(h=0)

# add residual diagnostics
text(20, 2.5, expression("Trend: No, "~italic("P")~"=0.85"), cex=1, pos=4)
text(20, 2.0, expression("Normality: No, "~italic("P")~"<0.01"), cex=1, pos=4)
text(20, 1.5, expression("Constant variance: No, "~italic("P")~"<0.01"), cex=1, pos=4)

#######################################
#########Time period 2#################
#######################################

#trend:
yr2<-seq(2003,2018,1)

summary(lm(slogres2~yr2)) #PASS
#normality:
shapiro.test(slogres2) #PASS
# constant variance:
library(car)
ncvTest(lm(slogres2~yr2)) 


#BARTLETT's test, variances are homogeneous if pval>0.05, split stdlogres into 2 time series")

test2 <- slogres2[!is.na(slogres2)]
if (length(test2) %% 2 == 0 ) {
  list1 <- test2[1:(length(test2)/2)]
  list2 <- test2[(length(test2)/2 +1):length(test2)]
}else{
  list1 <- test2[1:( (length(test2)+1)/2 )]
  list2 <- test2[( (length(test2)+1)/2 +1):length(test2)]
}
print( bartlett.test( list(list1, list2) ) )  #FAIL

#Try Bartlett's test removing large residuals
slogres2_red<-slogres2[-c(9,10,12,13,16)]   #-c(65,68) #need to remove several years still no pass
test2 <- slogres2_red[!is.na(slogres2_red)]
if (length(test2) %% 2 == 0 ) {
  list1 <- test2[1:(length(test2)/2)]
  list2 <- test2[(length(test2)/2 +1):length(test2)]
}else{
  list1 <- test2[1:( (length(test2)+1)/2 )]
  list2 <- test2[( (length(test2)+1)/2 +1):length(test2)]
}
print( bartlett.test( list(list1, list2) ) )  #FAIL



barplot(slogres2, axis.lty = 1,
        names.arg=yr2, xlab="Year", ylab="Standardized log residuals", xpd=F,
        ylim = c(floor(min(slogres1, na.rm = TRUE))-1,ceiling(max(slogres1, na.rm = TRUE))+1),
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xaxs = "i",yaxs = "i")
abline(h=0)

text(8, 2.5, expression("Trend: No, "~italic("P")~"=0.10"), cex=1, pos=4)
text(8, 2.0, expression("Normality: Yes, "~italic("P")~"=0.60"), cex=1, pos=4)
text(8, 1.5, expression("Constant variance: No, "~italic("P")~"=0.02"), cex=1, pos=4)


#BARTLETT's test, variances are homogeneous if pval>0.05, split stdlogres into 2 time series")

  test2 <- slogres2[!is.na(slogres2)]
if (length(test2) %% 2 == 0 ) {
list1 <- test2[1:(length(test2)/2)]
list2 <- test2[(length(test2)/2 +1):length(test2)]
}else{
list1 <- test2[1:( (length(test2)+1)/2 )]
list2 <- test2[( (length(test2)+1)/2 +1):length(test2)]
}
print( bartlett.test( list(list1, list2) ) )


##survey
barplot(slogres3) #not worth running tests
abline(h=0)

yr3<-seq(2017,2020,1)

barplot(slogres3, axis.lty = 1,
        names.arg=yr3, xlab="Year", ylab="Standardized log residuals", xpd=F,
        ylim = c(floor(min(slogres1, na.rm = TRUE))-1,ceiling(max(slogres1, na.rm = TRUE))+1),
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xaxs = "i",yaxs = "i")
