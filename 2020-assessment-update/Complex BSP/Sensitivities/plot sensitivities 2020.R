

library(R2WinBUGS)

#set up matrix to hold parameter estimates for table
# R,K,M,P1,MSY,BMSY,HMSY,B2015,H2015,poflH,poflB

par_mat<-array(NA,dim=c(40,13))

colnames(par_mat)<-c("R","K","M","P1","MSY","BMSY","HMSY","B2018","H2018","poflH","poflB","H/HMSY","B/BMSY")


#####Plot sensitivity of biomass to K#########################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_K14.5"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\origq_sens_k14.5"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

#pull out parameter estimates of interest
par_mat[1,"R"]<-stats[(rnames="r"),1]
par_mat[1,"K"]<-stats[(rnames="K"),1]
par_mat[1,"M"]<-stats[(rnames="M"),1]
par_mat[1,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[1,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[1,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[1,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[1,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[1,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[1,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[1,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[1,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[1,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_K21.75"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\origq_sens_k21.75"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

#pull out parameter estimates of interest
par_mat[2,"R"]<-stats[(rnames="r"),1]
par_mat[2,"K"]<-stats[(rnames="K"),1]
par_mat[2,"M"]<-stats[(rnames="M"),1]
par_mat[2,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[2,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[2,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[2,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[2,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[2,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[2,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[2,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[2,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[2,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_K36.25"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\origq_sens_k36.25"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

#pull out parameter estimates of interest
par_mat[3,"R"]<-stats[(rnames="r"),1]
par_mat[3,"K"]<-stats[(rnames="K"),1]
par_mat[3,"M"]<-stats[(rnames="M"),1]
par_mat[3,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[3,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[3,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[3,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[3,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[3,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[3,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[3,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[3,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[3,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_K43.5"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\origq_sens_k43.5"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

#pull out parameter estimates of interest
par_mat[4,"R"]<-stats[(rnames="r"),1]
par_mat[4,"K"]<-stats[(rnames="K"),1]
par_mat[4,"M"]<-stats[(rnames="M"),1]
par_mat[4,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[4,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[4,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[4,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[4,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[4,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[4,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[4,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[4,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[4,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

dev.new(height=4,width=6)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,35), main=expression(paste("Sensitivity to alternative prior mean for carrying capacity (",italic("K"),")")),
     xlab="Year", ylab="Mean Biomass (million lb)",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.5, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.5, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.5, pch=18)
lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")
legend(2000,10, c(expression(~italic("K")~"= 14.50"), expression(~italic("K")~"= 21.75"), expression(~italic("K")~"= 36.25"), expression(~italic("K")~"= 43.50")), 
       lty=c(1,1,1,1), pch=c(15,16,17,18), lwd=c(2,2,2,2), cex=1.25)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 35, by = 1), labels = FALSE, tcl = -0.3) 

expression(~italic("K")~"= 43.50")

#####Plot sensitivity of harvest rate to K#########################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_K14.5"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\origq_sens_k14.5"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_K21.75"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\origq_sens_k21.75"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_K36.25"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\origq_sens_k36.25"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_K43.5"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\origq_sens_k43.5"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly


dev.new(height=4,width=6)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,0.15), main=expression(paste("Sensitivity to alternative prior mean for carrying capacity (",italic("K"),")")),
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.5, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.5, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.5, pch=18)
lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")
legend(2000,0.14, c(expression(~italic("K")~"= 14.50"), expression(~italic("K")~"= 21.75"), expression(~italic("K")~"= 36.25"), expression(~italic("K")~"= 43.50")), 
       lty=c(1,1,1,1), pch=c(15,16,17,18), lwd=c(2,2,2,2), cex=1.25)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.15, by = 0.01), labels = FALSE, tcl = -0.3) 


#####Plot sensitivity of biomass to M#################################################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_M_0.5"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_M0.33" #mean M = 1.5
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


#pull out parameter estimates of interest
par_mat[5,"R"]<-stats[(rnames="r"),1]
par_mat[5,"K"]<-stats[(rnames="K"),1]
par_mat[5,"M"]<-stats[(rnames="M"),1]
par_mat[5,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[5,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[5,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[5,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[5,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[5,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[5,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[5,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[5,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[5,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_M_0.75"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_M0.4" #mean M = 1.25

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

#pull out parameter estimates of interest
par_mat[6,"R"]<-stats[(rnames="r"),1]
par_mat[6,"K"]<-stats[(rnames="K"),1]
par_mat[6,"M"]<-stats[(rnames="M"),1]
par_mat[6,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[6,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[6,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[6,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[6,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[6,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[6,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[6,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[6,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[6,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_M_1.25"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_M0.67"  #mean M = 0.75
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

#pull out parameter estimates of interest
par_mat[7,"R"]<-stats[(rnames="r"),1]
par_mat[7,"K"]<-stats[(rnames="K"),1]
par_mat[7,"M"]<-stats[(rnames="M"),1]
par_mat[7,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[7,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[7,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[7,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[7,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[7,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[7,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[7,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[7,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[7,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_M_1.5"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_M1"  #mean M = 0.5
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

#pull out parameter estimates of interest
par_mat[8,"R"]<-stats[(rnames="r"),1]
par_mat[8,"K"]<-stats[(rnames="K"),1]
par_mat[8,"M"]<-stats[(rnames="M"),1]
par_mat[8,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[8,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[8,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[8,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[8,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[8,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[8,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[8,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[8,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[8,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

dev.new(height=4,width=6)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,35), main=expression(paste("Sensitivity to alternative prior mean for shape parameter (",italic("M"),")")),
     xlab="Year", ylab="Mean Biomass (million lb)",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.5, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.5, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.5, pch=18)
lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")
legend(2000,34, c(expression(~italic("M")~"= 1.50"), expression(~italic("M")~"= 1.25"), expression(~italic("M")~"= 0.75"), expression(~italic("M")~"= 0.50")), 
       lty=c(1,1,1,1), pch=c(15,16,17,18), lwd=c(2,2,2,2), cex=1.25)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 35, by = 1), labels = FALSE, tcl = -0.3) 



#####Plot sensitivity of harvest rate to M#################################################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_M_0.5"

src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_M0.33" #mean M = 1.5
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_M_0.75_itau2init"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_M0.4" #mean M = 1.25
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_M_1.25"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_M0.67" #mean M = 0.75
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_M_1.5"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_M1"  #mean M = 0.5
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly



##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

dev.new(height=4,width=6)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,0.15), main=expression(paste("Sensitivity to alternative prior mean for shape parameter (",italic("M"),")")),
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.5, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.5, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.5, pch=18)
lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")
legend(2000,0.14, c(expression(~italic("M")~"= 1.50"), expression(~italic("M")~"= 1.25"), expression(~italic("M")~"= 0.75"), expression(~italic("M")~"= 0.50")), 
       lty=c(1,1,1,1), pch=c(15,16,17,18), lwd=c(2,2,2,2), cex=1.25)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.15, by = 0.01), labels = FALSE, tcl = -0.3) 





#####Plot sensitivity of biomass to r#########################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_r0.05"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_r0.05"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

#pull out parameter estimates of interest
par_mat[9,"R"]<-stats[(rnames="r"),1]
par_mat[9,"K"]<-stats[(rnames="K"),1]
par_mat[9,"M"]<-stats[(rnames="M"),1]
par_mat[9,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[9,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[9,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[9,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[9,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[9,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[9,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[9,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[9,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[9,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_r0.15"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_r0.15"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


par_mat[10,"R"]<-stats[(rnames="r"),1]
par_mat[10,"K"]<-stats[(rnames="K"),1]
par_mat[10,"M"]<-stats[(rnames="M"),1]
par_mat[10,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[10,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[10,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[10,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[10,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[10,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[10,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[10,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[10,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[10,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_r0.25"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_r0.25"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


par_mat[11,"R"]<-stats[(rnames="r"),1]
par_mat[11,"K"]<-stats[(rnames="K"),1]
par_mat[11,"M"]<-stats[(rnames="M"),1]
par_mat[11,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[11,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[11,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[11,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[11,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[11,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[11,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[11,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[11,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[11,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]



##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

dev.new(height=4,width=6)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,35), main=expression(paste("Sensitivity to alternative prior mean for intrinsic growth rate (",italic("R"),")")),
     xlab="Year", ylab="Mean Biomass (million lb)",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.5, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.5, pch=17)

lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")
legend(2000,34, c(expression(~italic("R")~"= 0.05"), expression(~italic("R")~"= 0.15"), expression(~italic("R")~"= 0.25")), 
       lty=c(1,1,1), pch=c(15,16,17), lwd=c(2,2,2), cex=1.25)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 35, by = 1), labels = FALSE, tcl = -0.3) 


#####Plot sensitivity of harvest rate to r#########################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_r0.05"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_r0.05"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_r0.15"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_r0.15"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_r0.25"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_r0.25"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly





##grab the default biomass estimates
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

dev.new(height=4,width=6)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,0.15), main=expression(paste("Sensitivity to alternative prior mean for intrinsic growth rate (",italic("R"),")")),
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.5, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.5, pch=17)

lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")
legend(2000,0.14, c(expression(~italic("R")~"= 0.05"), expression(~italic("R")~"= 0.15"), expression(~italic("R")~"= 0.25")), 
       lty=c(1,1,1), pch=c(15,16,17), lwd=c(2,2,2), cex=1.25)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.14, by = 0.01), labels = FALSE, tcl = -0.3) 



#####Plot sensitivity of biomass to P1#########################################################
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_P1_0.265"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_p1_0.265"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

par_mat[12,"R"]<-stats[(rnames="r"),1]
par_mat[12,"K"]<-stats[(rnames="K"),1]
par_mat[12,"M"]<-stats[(rnames="M"),1]
par_mat[12,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[12,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[12,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[12,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[12,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[12,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[12,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[12,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[12,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[12,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_P1_0.3975"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_p1_0.4"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


par_mat[13,"R"]<-stats[(rnames="r"),1]
par_mat[13,"K"]<-stats[(rnames="K"),1]
par_mat[13,"M"]<-stats[(rnames="M"),1]
par_mat[13,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[13,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[13,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[13,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[13,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[13,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[13,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[13,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[13,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[13,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_P1_0.6625"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_p1_0.66"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


par_mat[14,"R"]<-stats[(rnames="r"),1]
par_mat[14,"K"]<-stats[(rnames="K"),1]
par_mat[14,"M"]<-stats[(rnames="M"),1]
par_mat[14,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[14,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[14,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[14,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[14,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[14,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[14,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[14,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[14,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[14,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_P1_0.795"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_p1_0.795"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


par_mat[15,"R"]<-stats[(rnames="r"),1]
par_mat[15,"K"]<-stats[(rnames="K"),1]
par_mat[15,"M"]<-stats[(rnames="M"),1]
par_mat[15,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[15,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[15,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[15,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[15,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[15,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[15,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[15,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[15,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[15,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

dev.new(height=4,width=6)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,35), main=expression(paste("Sensitivity to alternative prior mean for initial proportion of carrying capacity (",italic("P"[italic("1")]),")")),
     xlab="Year", ylab="Mean Biomass (million lb)",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.5, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.5, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.5, pch=18)

lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")
legend(2000,34, c(expression(~italic("P"[italic("1")])~"= 0.27"), expression(~italic("P"[italic("1")])~"= 0.40"), expression(~italic("P"[italic("1")])~"= 0.66"),expression(~italic("P"[italic("1")])~"= 0.80")), 
       lty=c(1,1,1,1), pch=c(15,16,17,18), lwd=c(2,2,2,2), cex=1.25)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 35, by = 1), labels = FALSE, tcl = -0.3) 


#####Plot sensitivity of Harvest Rate to P1#########################################################
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_P1_0.265"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_p1_0.265"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_P1_0.3975"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_p1_0.4"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_P1_0.6625"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_p1_0.66"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_P1_0.795"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_p1_0.795"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

dev.new(height=4,width=6)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,0.15), main=expression(paste("Sensitivity to alternative prior mean for initial proportion of carrying capacity (",italic("P"[italic("1")]),")")),
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.5, pch=15, cex.lab=1.25, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.25, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.25, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.25, pch=18)

lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")
legend(2000,0.14, c(expression(~italic("P"[italic("1")])~"= 0.27"), expression(~italic("P"[italic("1")])~"= 0.40"), expression(~italic("P"[italic("1")])~"= 0.66"),expression(~italic("P"[italic("1")])~"= 0.80")), 
       lty=c(1,1,1,1), pch=c(15,16,17,18), lwd=c(2,2,2,2), cex=1.25)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.15, by = 0.01), labels = FALSE, tcl = -0.3) 




#####################Sensitivity of Biomass to inclusion of the survey##################################


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_nosurvey"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Retrospective\\d7_2020_retro_NOSURVEY\\d7_retro0"
setwd(src.dir)
logtxt <-  bugs.log("d7_retro0_log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly
bio_1_low<- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),4] ## <--- Change accordingly
bio_1_hi<- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),6] ## <--- Change accordingly
bio_1_sd<-stats[which(rnames=="B[1]"):which(rnames=="B[70]"),2] ## <--- Change accordingly

par_mat[16,"R"]<-stats[(rnames="r"),1]
par_mat[16,"K"]<-stats[(rnames="K"),1]
par_mat[16,"M"]<-stats[(rnames="M"),1]
par_mat[16,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[16,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[16,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[16,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[16,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[16,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[16,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]  #2020 - need to monitor!
par_mat[16,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]  #2020 - need to monitor!
par_mat[16,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[16,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly
bio_base_low<- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),4] ## <--- Change accordingly
bio_base_hi<- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),6] ## <--- Change accordingly
bio_base_sd <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),2] ## <--- Change accordingly

dev.new(height=4,width=6)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,75), main="Sensitivity to inclusion of fishery-independent survey",
     xlab="Year", ylab="Mean Biomass (million lb)",
     type="l", lwd=2,col="red",cex=1.5,  cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
polygon(c(seq(1949,2018,1),rev(seq(1949,2018,1))),c(bio_1_low,rev(bio_1_hi)),col=rgb(1,0,0,0.2), border = "NA")
lines(bio_1~seq(1949,2018,1),col="red") #,type="o"
lines(seq(1949,2018,1),bio_1_low,lty=5,col=rgb(1,0,0,0.5),lwd=3)
lines(seq(1949,2018,1),bio_1_hi,lty=5,col=rgb(1,0,0,0.5),lwd=3)

lines(bio_base~seq(1949,2018,1),lwd=4,col="blue")
polygon(c(seq(1949,2018,1),rev(seq(1949,2018,1))),c(bio_base_low,rev(bio_base_hi)),col=rgb(0,0,1,0.2), border = "NA")
lines(bio_base~seq(1949,2018,1),lwd=4, col="blue")
lines(seq(1949,2018,1),bio_base_low,lty=3,col=rgb(0,0,1,0.5),lwd=3)
lines(seq(1949,2018,1),bio_base_hi,lty=3,col=rgb(0,0,1,0.5),lwd=3)

legend(1985,74, c("Survey excluded - mean", "Survey excluded - 95% CI", "Survey included - mean","Survey included - 95% CI"), 
       lty=c(1,5,1,3), pch=c(15,NA,NA,NA), lwd=c(2,2), cex=1.25,col=c("red","red","blue","blue"))


axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 75, by = 1), labels = FALSE, tcl = -0.3) 


abline(h=16.46*0.844 , lwd=3, lty=3, col="red")
text(2014,16,expression("0.844*"~B[MSY]), cex=1.25, col="red")
abline(h=15.46*0.844,lwd=3, lty=5,col="blue")
text(2014,11,expression("0.844*"~B[MSY]),cex=1.25,col="blue")

#calculate CI width differences to mention in text?
ci1<-bio_1_hi-bio_1_low
cibase<-bio_base_hi-bio_base_low

ci1/bio_1

cibase/bio_base


a<-bio_base_sd/bio_base

b<-bio_1_sd/bio_1


#####################Sensitivity of Harvest to inclusion of the survey##################################


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_nosurvey"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Retrospective\\d7_2020_retro_NOSURVEY\\d7_retro0"
setwd(src.dir)
logtxt <-  bugs.log("d7_retro0_log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly
bio_1_low<- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),4] ## <--- Change accordingly
bio_1_hi<- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),6] ## <--- Change accordingly


##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly
bio_base_low<- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),4] ## <--- Change accordingly
bio_base_hi<- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),6] ## <--- Change accordingly


dev.new(height=4,width=6)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,0.25), main="Sensitivity to inclusion of fishery-independent survey",
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i",col="red")
polygon(c(seq(1949,2018,1),rev(seq(1949,2018,1))),c(bio_1_low,rev(bio_1_hi)),col=rgb(1,0,0,0.2), border = "NA")
lines(bio_1~seq(1949,2018,1),type="o",col="red")
lines(seq(1949,2018,1),bio_1_low,lty=5,col=rgb(1,0,0,0.5),lwd=3)
lines(seq(1949,2018,1),bio_1_hi,lty=5,col=rgb(1,0,0,0.5),lwd=3)

lines(bio_base~seq(1949,2018,1),lwd=4, col="blue")
polygon(c(seq(1949,2018,1),rev(seq(1949,2018,1))),c(bio_base_low,rev(bio_base_hi)),col=rgb(0,0,1,0.2), border = "NA")
lines(bio_base~seq(1949,2018,1),lwd=4, col="blue")
lines(seq(1949,2018,1),bio_base_low,lty=3,col=rgb(0,0,1,0.5),lwd=3)
lines(seq(1949,2018,1),bio_base_hi,lty=3,col=rgb(0,0,1,0.5),lwd=3)

legend(1990,0.24, c("Survey excluded - mean", "Survey excluded - 95% CI", "Survey included - mean","Survey included - 95% CI"), 
       lty=c(1,5,1,3), pch=c(15,NA,NA,NA), lwd=c(2,2), cex=1.25,col=c("red","red","blue","blue"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.15, by = 0.01), labels = FALSE, tcl = -0.3) 


abline(h=0.068,lwd=3, lty=5,col="blue")
text(1962,0.0633,expression(H[MSY]), cex=1.5, col="blue")
abline(h=0.069 , lwd=3, lty=3, col="red")
text(1962,0.074,expression(H[MSY]), cex=1.5, col="red")




#####################Sensitivity of Biomass to reduced uncertainty in survey radius##################################


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_survey"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_survey"


setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly
bio_1_low<- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),4] ## <--- Change accordingly
bio_1_hi<- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),6] ## <--- Change accordingly
bio_1_sd<-stats[which(rnames=="B[1]"):which(rnames=="B[70]"),2] ## <--- Change accordingly

par_mat[17,"R"]<-stats[(rnames="r"),1]
par_mat[17,"K"]<-stats[(rnames="K"),1]
par_mat[17,"M"]<-stats[(rnames="M"),1]
par_mat[17,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[17,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[17,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[17,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[17,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[17,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[17,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[17,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[17,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[17,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly
bio_base_low<- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),4] ## <--- Change accordingly
bio_base_hi<- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),6] ## <--- Change accordingly
bio_base_sd <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),2] ## <--- Change accordingly

dev.new(height=4,width=6)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,65), main="Sensitivity to decreased uncertainty for fishery-independent survey",
     xlab="Year", ylab="Mean Biomass (million lb)",
     type="o", col="red",cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
polygon(c(seq(1949,2018,1),rev(seq(1949,2018,1))),c(bio_1_low,rev(bio_1_hi)),col=rgb(1,0,0,0.2), border = "NA")
lines(bio_1~seq(1949,2018,1),type="o",col="red")
lines(seq(1949,2018,1),bio_1_low,lty=5,col=rgb(1,0,0,0.5),lwd=3)
lines(seq(1949,2018,1),bio_1_hi,lty=5,col=rgb(1,0,0,0.5),lwd=3)

lines(bio_base~seq(1949,2018,1),lwd=4,col="blue")
polygon(c(seq(1949,2018,1),rev(seq(1949,2018,1))),c(bio_base_low,rev(bio_base_hi)),col=rgb(0,0,1,0.2), border = "NA")
lines(bio_base~seq(1949,2018,1),lwd=4, col="blue")
lines(seq(1949,2018,1),bio_base_low,lty=3,col=rgb(0,0,1,0.5),lwd=3)
lines(seq(1949,2018,1),bio_base_hi,lty=3,col=rgb(0,0,1,0.5),lwd=3)

legend(1980,64, c("Survey CV reduced - mean", "Survey CV reduced - 95% CI", "Survey base - mean","Survey base - 95% CI"), 
       lty=c(1,5,1,3), pch=c(15,NA,NA,NA), lwd=c(2,2), cex=1.25,col=c("red","red","blue","blue"))


axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 65, by = 1), labels = FALSE, tcl = -0.3) 

abline(h=12.72*0.844,lwd=3, lty=5,col="red")
text(1960,10,expression("0.844*"~B[MSY]),cex=1.5,col="red")
abline(h=15.46*0.844 , lwd=3, lty=3, col="blue")
text(1960,14,expression("0.844*"~B[MSY]), cex=1.5, col="blue")


#calculate CI width differences to mention in text?
ci1<-bio_1_hi-bio_1_low
cibase<-bio_base_hi-bio_base_low

ci1/bio_1

cibase/bio_base


a<-bio_base_sd/bio_base

b<-bio_1_sd/bio_1


#####################Sensitivity of Harvest to reduced uncertainty in survey radius##################################


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_survey"

src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_survey"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly
bio_1_low<- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),4] ## <--- Change accordingly
bio_1_hi<- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),6] ## <--- Change accordingly


##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly
bio_base_low<- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),4] ## <--- Change accordingly
bio_base_hi<- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),6] ## <--- Change accordingly


dev.new(height=4,width=6)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,0.25), main="Sensitivity to decreased uncertainty for fishery-independent survey",
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i",col="red")
polygon(c(seq(1949,2018,1),rev(seq(1949,2018,1))),c(bio_1_low,rev(bio_1_hi)),col=rgb(1,0,0,0.2), border = "NA")
lines(bio_1~seq(1949,2018,1),type="o",col="red")
lines(seq(1949,2018,1),bio_1_low,lty=5,col=rgb(1,0,0,0.5),lwd=3)
lines(seq(1949,2018,1),bio_1_hi,lty=5,col=rgb(1,0,0,0.5),lwd=3)

lines(bio_base~seq(1949,2018,1),lwd=4, col="blue")
polygon(c(seq(1949,2018,1),rev(seq(1949,2018,1))),c(bio_base_low,rev(bio_base_hi)),col=rgb(0,0,1,0.2), border = "NA")
lines(bio_base~seq(1949,2018,1),lwd=4, col="blue")
lines(seq(1949,2018,1),bio_base_low,lty=3,col=rgb(0,0,1,0.5),lwd=3)
lines(seq(1949,2018,1),bio_base_hi,lty=3,col=rgb(0,0,1,0.5),lwd=3)

legend(1946,0.24, c("Survey CV reduced - mean", "Survey CV reduced - 95% CI", "Survey base - mean","Survey base - 95% CI"), 
       lty=c(1,5,1,3), pch=c(15,NA,NA,NA), lwd=c(2,2), cex=1.25,col=c("red","red","blue","blue"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.25, by = 0.01), labels = FALSE, tcl = -0.3) 


abline(h=0.068,lwd=3, lty=5,col="blue")
text(1962,0.074,expression(H[MSY]), cex=1.5, col="blue")
abline(h=0.06591 , lwd=3, lty=3, col="red")
text(1961,0.06,expression(H[MSY]), cex=1.5, col="red")



#####################Sensitivity of Biomass to uniform versus inverse guassian prior on error##################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_uniferr"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_uniferr"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

par_mat[18,"R"]<-stats[(rnames="r"),1]
par_mat[18,"K"]<-stats[(rnames="K"),1]
par_mat[18,"M"]<-stats[(rnames="M"),1]
par_mat[18,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[18,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[18,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[18,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[18,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[18,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[18,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[18,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[18,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[18,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]


##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


dev.new(wicth=6,height=4)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,35), main="Sensitivity to uniform versus inverse gamma prior distributions 
     for observation and process errors",
     xlab="Year", ylab="Mean Biomass (million lb)",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")

lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")
legend(1985,34, c("Uniform prior"), 
       lty=c(1,1,1,1), pch=c(15), lwd=c(2), cex=1.25)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 35, by = 1), labels = FALSE, tcl = -0.3) 


#####################Sensitivity of Harvest Rate to uniform versus inverse guassian prior on error##################################
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_uniferr"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_uniferr"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly


##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

dev.new(width=6,height=4)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,0.12), main="Sensitivity to uniform versus inverse gamma prior distributions 
     for observation and process errors",
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")

lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")
legend(1950,0.115, c("Uniform prior"), 
       lty=c(1,1,1,1), pch=c(15), lwd=c(2), cex=1.25)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.15, by = 0.01), labels = FALSE, tcl = -0.3) 


#####Plot sensitivity of biomass to process error mode#########################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_sigma2_0.001"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_sig0.001"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


par_mat[19,"R"]<-stats[(rnames="r"),1]
par_mat[19,"K"]<-stats[(rnames="K"),1]
par_mat[19,"M"]<-stats[(rnames="M"),1]
par_mat[19,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[19,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[19,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[19,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[19,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[19,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[19,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[19,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[19,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[19,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_sigma2_0.01"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_sig0.01"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


par_mat[20,"R"]<-stats[(rnames="r"),1]
par_mat[20,"K"]<-stats[(rnames="K"),1]
par_mat[20,"M"]<-stats[(rnames="M"),1]
par_mat[20,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[20,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[20,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[20,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[20,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[20,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[20,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[20,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[20,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[20,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_sigma2_1"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_sig1"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


par_mat[21,"R"]<-stats[(rnames="r"),1]
par_mat[21,"K"]<-stats[(rnames="K"),1]
par_mat[21,"M"]<-stats[(rnames="M"),1]
par_mat[21,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[21,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[21,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[21,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[21,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[21,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[21,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[21,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[21,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[21,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_sigma2_10_q bounds reduced"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_sig10"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

par_mat[22,"R"]<-stats[(rnames="r"),1]
par_mat[22,"K"]<-stats[(rnames="K"),1]
par_mat[22,"M"]<-stats[(rnames="M"),1]
par_mat[22,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[22,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[22,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[22,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[22,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[22,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[22,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[22,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[22,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[22,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

dev.new(width=6,height=4)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,40), main="Sensitivity to alternative prior mode for process error variance",
     xlab="Year", ylab="Mean Biomass (million lb)",
     type="o", cex=1.25, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.25, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.25, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.25, pch=18)
lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")


a<-expression(paste(sigma^2,' = 0.83 x ', 10^-3))
b<-expression(paste(sigma^2,' = 0.83 x ', 10^-2))
c<-expression(paste(sigma^2,' = 0.83 x ', 10^-0))
d<-expression(paste(sigma^2,' = 0.83 x ', 10^1))

legend(1990,35, c(a,b,c,d), lty=c(1,1,1,1), pch=c(15,16,17,18), lwd=c(2,2,2,2), cex=1.25, bty="n") 


axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 40, by = 1), labels = FALSE, tcl = -0.3) 


#####Plot sensitivity of harvest rate to process error mode#########################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_sigma2_0.001"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_sig0.001"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_sigma2_0.01"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_sig0.01"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_sigma2_1"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_sig1"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_sigma2_10_q bounds reduced"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_sig10"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly


dev.new(width=6,height=4)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,0.15), main="Sensitivity to alternative prior mode for process error variance",
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.25, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.25, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.25, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.25, pch=18)
lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")
a<-expression(paste(sigma^2,' = 0.83 x ', 10^-3))
b<-expression(paste(sigma^2,' = 0.83 x ', 10^-2))
c<-expression(paste(sigma^2,' = 0.83 x ', 10^-0))
d<-expression(paste(sigma^2,' = 0.83 x ', 10^1))

legend(1950,0.13, c(a,b,c,d), lty=c(1,1,1,1), pch=c(15,16,17,18), lwd=c(2,2,2,2), cex=1.25, bty="n") 
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.15, by = 0.01), labels = FALSE, tcl = -0.3) 





#####Plot sensitivity of biomass to observation error mode#########################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_tau2_0.01"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_tau0.01"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

par_mat[23,"R"]<-stats[(rnames="r"),1]
par_mat[23,"K"]<-stats[(rnames="K"),1]
par_mat[23,"M"]<-stats[(rnames="M"),1]
par_mat[23,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[23,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[23,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[23,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[23,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[23,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[23,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[23,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[23,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[23,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_tau2_0.1"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_tau0.1"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

par_mat[24,"R"]<-stats[(rnames="r"),1]
par_mat[24,"K"]<-stats[(rnames="K"),1]
par_mat[24,"M"]<-stats[(rnames="M"),1]
par_mat[24,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[24,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[24,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[24,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[24,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[24,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[24,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[24,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[24,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[24,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_tau2_10"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_tau10"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

par_mat[25,"R"]<-stats[(rnames="r"),1]
par_mat[25,"K"]<-stats[(rnames="K"),1]
par_mat[25,"M"]<-stats[(rnames="M"),1]
par_mat[25,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[25,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[25,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[25,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[25,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[25,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[25,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[25,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[25,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[25,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_tau2_100"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_tau100"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

par_mat[26,"R"]<-stats[(rnames="r"),1]
par_mat[26,"K"]<-stats[(rnames="K"),1]
par_mat[26,"M"]<-stats[(rnames="M"),1]
par_mat[26,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[26,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[26,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[26,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[26,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[26,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[26,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[26,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[26,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[26,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]


##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

dev.new(width=6,height=4)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,35), main="Sensitivity to alternative prior mode for observation error variance",
     xlab="Year", ylab="Mean Biomass (million lb)",
     type="o", cex=1.25, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.25, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.25, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.25, pch=18)
lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")


a<-expression(paste(tau^2,' = 0.83 x ', 10^-2))
b<-expression(paste(tau^2,' = 0.83 x ', 10^-1))
c<-expression(paste(tau^2,' = 0.83 x ', 10^1))
d<-expression(paste(tau^2,' = 0.83 x ', 10^2))

legend(1991,34, c(a,b,c,d), lty=c(1,1,1,1), pch=c(15,16,17,18), lwd=c(2,2,2,2), cex=1.25, bty="n") 


axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 35, by = 1), labels = FALSE, tcl = -0.3) 


#####Plot sensitivity of harvest rate to observation error mode#########################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_tau2_0.01"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_tau0.01"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_tau2_0.1"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_tau0.1"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_tau2_10"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_tau10"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_tau2_100"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_tau100"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly


dev.new(width=6,height=4)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,0.15), main="Sensitivity to alternative prior mode for observation error variance",
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.5, pch=15, cex.lab=1.25, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.25, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.25, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.25, pch=18)
lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")
a<-expression(paste(tau^2,' = 0.83 x ', 10^-2))
b<-expression(paste(tau^2,' = 0.83 x ', 10^-1))
c<-expression(paste(tau^2,' = 0.83 x ', 10^1))
d<-expression(paste(tau^2,' = 0.83 x ', 10^2))

legend(1950,0.14, c(a,b,c,d), lty=c(1,1,1,1), pch=c(15,16,17,18), lwd=c(2,2,2,2), cex=1.25, bty="n") 
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.15, by = 0.01), labels = FALSE, tcl = -0.3) 



#####Plot sensitivity of biomass to alternative unreported catch ratios #########################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC1"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_unrepC1"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

par_mat[27,"R"]<-stats[(rnames="r"),1]
par_mat[27,"K"]<-stats[(rnames="K"),1]
par_mat[27,"M"]<-stats[(rnames="M"),1]
par_mat[27,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[27,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[27,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[27,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[27,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[27,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[27,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[27,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[27,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[27,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC3"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_unrepC3"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

par_mat[28,"R"]<-stats[(rnames="r"),1]
par_mat[28,"K"]<-stats[(rnames="K"),1]
par_mat[28,"M"]<-stats[(rnames="M"),1]
par_mat[28,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[28,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[28,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[28,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[28,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[28,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[28,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[28,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[28,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[28,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC4"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_unrepC4"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

par_mat[29,"R"]<-stats[(rnames="r"),1]
par_mat[29,"K"]<-stats[(rnames="K"),1]
par_mat[29,"M"]<-stats[(rnames="M"),1]
par_mat[29,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[29,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[29,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[29,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[29,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[29,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[29,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[29,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[29,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[29,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC5"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_unrepC5"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

par_mat[30,"R"]<-stats[(rnames="r"),1]
par_mat[30,"K"]<-stats[(rnames="K"),1]
par_mat[30,"M"]<-stats[(rnames="M"),1]
par_mat[30,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[30,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[30,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[30,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[30,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[30,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[30,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[30,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[30,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[30,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]


##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

dev.new(width=6,height=4)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,35), main="Sensitivity to alternative unreported catch ratios",
     xlab="Year", ylab="Mean Biomass (million lb)",
     type="o", cex=1.25, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.25, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.25, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.25, pch=18)
lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")



legend(1990,34, c("Scenario I", "Scenario II", "Scenario III", "Scenario IV"), lty=c(1,1,1,1), pch=c(15,16,17,18), lwd=c(2,2,2,2), cex=1.25) 


axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 35, by = 1), labels = FALSE, tcl = -0.3) 

######REDO for CIE presentation####################################
#BIOMASS
src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC1"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC3"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC4"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC5"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

##grab the default biomass estimates
src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,30), main="Sensitivity to alternative unreported catch ratios",
     xlab="Year", ylab="Mean Biomass (millions of lbs.)",
     type="o", cex=1.5, pch=18, col = 6, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, lwd=2.5, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.5, pch=17, col = 3,lwd=2.5)
points(bio_3~seq(1949,2018,1), type="o",cex=1.5, pch=20, col = 4, lwd=2.5)
points(bio_4~seq(1949,2018,1), type="o",cex=1.5, pch=15, col = 8, lwd=2.5)
lines(bio_base~seq(1949,2018,1), type="o",pch=3, col=1,lwd=2.5)



legend(1990,29, c("Base", "Scen1", "Scen2", "Scen3", "Scen4"), lty=c(1,1,1,1,1), pch=c(3,18,17,20,15), lwd=c(3,3,3,3,3), col=c(1,6,3,4,8),cex=1.7) 


axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 30, by = 1), labels = FALSE, tcl = -0.3) 


abline(h=0.844*16.24,lwd=2.5,col=6)
abline(h=0.844*13.26,lwd=2.5,col=3)
abline(h=0.844*12.9,lwd=2.5,col=4)
abline(h=0.844*10.33,lwd=2.5,col=8)
abline(h=0.844*13.69,lwd=2.5,col=1)

#Harvest rate#######################################################################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC1"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_unrepC1"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC3"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_unrepC3"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC4"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_unrepC4"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC5"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_unrepC5"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly


dev.new(width=6,height=4)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,0.2), main="Sensitivity to alternative unreported catch ratios",
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.25,cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i", pch=18, col = 6, lwd=2.5)
points(bio_2~seq(1949,2018,1), type="o",cex=1.25, pch=17, col = 3, lwd=2.5)
points(bio_3~seq(1949,2018,1), type="o",cex=1.25, pch=20, col = 4, lwd=2.5)
points(bio_4~seq(1949,2018,1), type="o",cex=1.25, pch=15, col = 8, lwd=2.5)
lines(bio_base~seq(1949,2018,1), pch=3, col = 1, lwd=2.5)


legend(1990,29, c("Base", "Scen1", "Scen2", "Scen3", "Scen4"), lty=c(1,1,1,1,1), pch=c(3,18,17,20,15), lwd=c(3,3,3,3,3), col=c(1,6,3,4,8),cex=1.7) 

axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.2, by = 0.01), labels = FALSE, tcl = -0.3) 




abline(h=0.07057,lwd=2.5,col=6)
abline(h=0.06514,lwd=2.5,col=3)
abline(h=0.06394,lwd=2.5,col=4)
abline(h=0.04629,lwd=2.5,col=8)
abline(h=0.066,lwd=2.5,col=1)

#####Plot sensitivity of harvest rate to alternative unreported catch ratios#########################################################
#remove scenario II 10-6-17
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC1"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_unrepC1"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC3"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_unrepC3"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC4"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_unrepC4"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_unrepC5"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sens_unrepC5"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

##grab the default biomass estimates
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly


dev.new(width=6,height=4)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,0.15), main="Sensitivity to alternative unreported catch ratios",
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.25, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.25, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.25, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.25, pch=18)
lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")


legend(1995,0.14, c("Scenario I", "Scenario II", "Scenario III", "Scenario IV"), lty=c(1,1,1,1), pch=c(15,16,17,18), lwd=c(2,2,2,2), cex=1.25) 
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.15, by = 0.01), labels = FALSE, tcl = -0.3) 



#####Plot sensitivity of biomass to alternative unreported catch error  #########################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_UCe_0.01"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_unrepCe0.01"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

par_mat[31,"R"]<-stats[(rnames="r"),1]
par_mat[31,"K"]<-stats[(rnames="K"),1]
par_mat[31,"M"]<-stats[(rnames="M"),1]
par_mat[31,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[31,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[31,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[31,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[31,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[31,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[31,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[31,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[31,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[31,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_UCe_0.2"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_unrepCe0.2"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


par_mat[32,"R"]<-stats[(rnames="r"),1]
par_mat[32,"K"]<-stats[(rnames="K"),1]
par_mat[32,"M"]<-stats[(rnames="M"),1]
par_mat[32,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[32,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[32,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[32,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[32,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[32,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[32,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[32,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[32,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[32,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_UCe_0.6"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_unrepCe0.6"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


par_mat[33,"R"]<-stats[(rnames="r"),1]
par_mat[33,"K"]<-stats[(rnames="K"),1]
par_mat[33,"M"]<-stats[(rnames="M"),1]
par_mat[33,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[33,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[33,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[33,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[33,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[33,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[33,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[33,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[33,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[33,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_UCe_neg25"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\unrepCe_neg25"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


par_mat[34,"R"]<-stats[(rnames="r"),1]
par_mat[34,"K"]<-stats[(rnames="K"),1]
par_mat[34,"M"]<-stats[(rnames="M"),1]
par_mat[34,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[34,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[34,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[34,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[34,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[34,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[34,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[34,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[34,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[34,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_UCe_pos25"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\unrepCe_pos25"

setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_5 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly


par_mat[35,"R"]<-stats[(rnames="r"),1]
par_mat[35,"K"]<-stats[(rnames="K"),1]
par_mat[35,"M"]<-stats[(rnames="M"),1]
par_mat[35,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[35,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[35,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[35,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[35,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[35,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[35,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[35,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[35,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[35,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]


##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

dev.new(width=6,height=4)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,30), main="Sensitivity to alternative error distributions for unreported catch",
     xlab="Year", ylab="Mean Biomass (million lb)",
     type="o", cex=1.25, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.25, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.25, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.25, pch=18)
points(bio_5~seq(1949,2018,1), type="o",cex=1.25, pch=4)
lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")



legend(1990,14, c("[0.9999, 1.0001]", "[0.80,1.20]", "[0.40, 1.60]", "25% Decrease", "25% Increase" ), lty=c(1,1,1,1,1), pch=c(15,16,17,18,4), 
       lwd=c(2,2,2,2,2), cex=1.25) 


axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 30, by = 1), labels = FALSE, tcl = -0.3) 


#####Plot sensitivity of harvest rate to alternative unreported catch error#########################################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_UCe_0.01"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_unrepCe0.01"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_UCe_0.2"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_unrepCe0.2"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_UCe_0.6"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_unrepCe0.6"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_3 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_UCe_neg25"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\unrepCe_neg25"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_4 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly


#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_UCe_pos25"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\unrepCe_pos25"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_5 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly





##grab the default biomass estimates
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly


dev.new(width=6,height=4)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,0.15), main="Sensitivity to alternative error distributions for unreported catch",
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.25, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
points(bio_2~seq(1949,2018,1), type="o",cex=1.25, pch=16)
points(bio_3~seq(1949,2018,1), type="o",cex=1.25, pch=17)
points(bio_4~seq(1949,2018,1), type="o",cex=1.25, pch=18)
points(bio_5~seq(1949,2018,1), type="o",cex=1.25, pch=4)
lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")


legend(1950,0.14, c("[0.9999, 1.0001]", "[0.80,1.20]", "[0.40, 1.60]", "25% Decrease", "25% Increase" ), lty=c(1,1,1,1), pch=c(15,16,17,18), lwd=c(2,2,2,2), cex=1.25) 
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.15, by = 0.01), labels = FALSE, tcl = -0.3) 




#####################Sensitivity of Biomass to time-varying catchability##################################

#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_varq"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_varq"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

par_mat[36,"R"]<-stats[(rnames="r"),1]
par_mat[36,"K"]<-stats[(rnames="K"),1]
par_mat[36,"M"]<-stats[(rnames="M"),1]
par_mat[36,"P1"]<-stats[(rnames="P[1]"),1]
par_mat[36,"MSY"]<-stats[(rnames="MSY"),1]
par_mat[36,"BMSY"]<-stats[(rnames="BMSY"),1]
par_mat[36,"HMSY"]<-stats[(rnames="HMSY"),1]
par_mat[36,"B2018"]<-stats[(rnames="B[70]"),1]
par_mat[36,"H2018"]<-stats[(rnames="H[70]"),1]
par_mat[36,"poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_mat[36,"poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_mat[36,"H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_mat[36,"B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]



#END of par_mat -- Clean Up and Export
par_mat<-par_mat[1:36,]

rownames(par_mat)<-c("K = 14.50", "K = 21.75", "K = 36.25", "K = 43.50", "M = 1.50", "M = 1.25",
                    "M = 0.75","M = 0.50", "R = 0.05", "R = 0.15","R = 0.25","P1 = 0.22","P1 = 0.33",
                    "P1 = 0.55", "P1 = 0.66", "No survey", "Survey CV", "Uniform error", "Sigma x 0.01",
                    "Sigma x 0.1","Sigma x 10","Sigma x 100","Tau x 0.01","Tau x 0.1","Tau x 10",
                    "Tau x 100","Catch I","Catch II","Catch III","Catch IV", "CU 0.01", "CU 0.2", "CU 0.6",
                    "CU neg25", "CU pos25", "q")

#bring in base results to compute % change
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)

par_base<-as.vector(rep(NA,13))
names(par_base)<-c("R","K","M","P1","MSY","BMSY","HMSY","B2018","H2018","poflH","poflB","H/HMSY","B/BMSY")


par_base["R"]<-stats[(rnames="r"),1]
par_base["K"]<-stats[(rnames="K"),1]
par_base["M"]<-stats[(rnames="M"),1]
par_base["P1"]<-stats[(rnames="P[1]"),1]
par_base["MSY"]<-stats[(rnames="MSY"),1]
par_base["BMSY"]<-stats[(rnames="BMSY"),1]
par_base["HMSY"]<-stats[(rnames="HMSY"),1]
par_base["B2018"]<-stats[(rnames="B[70]"),1]
par_base["H2018"]<-stats[(rnames="H[70]"),1]
par_base["poflH"]<-stats[(rnames="pOFL_H[70]"),1]
par_base["poflB"]<-stats[(rnames="pOFL_B[70]"),1]
par_base["H/HMSY"]<-stats[(rnames="HSTATUS[70]"),5]
par_base["B/BMSY"]<-stats[(rnames="BSTATUS[70]"),1]

#compute proportional change
par_mat_rel<-apply(par_mat,1,function(x) (x-par_base)/par_base)

write.csv(t(par_mat_rel),file="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_res_table.csv")
write.csv(t(par_mat),file="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_res_table_absolute.csv")




##grab the default biomass estimates
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2017_baseWPSAR"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[70]"),1] ## <--- Change accordingly

dev.new(width=6, height=4)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,40), main="Sensitivity to inclusion of random-walk catchability",
     xlab="Year", ylab="Mean Biomass (million lb)",
     type="o", cex=1.25, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")

lines(bio_base~seq(1949,2018,1),lwd=4, col="darkgray")
legend(1950,9, c("Random-walk catchability"), 
       lty=c(1,1,1,1), pch=c(15), lwd=c(2), cex=1.75)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 35, by = 1), labels = FALSE, tcl = -0.3) 

#dont include msst
abline(h=0.844*18.46,lwd=2)
abline(h=0.844*15.46 , lwd=2,col="darkgray")
text(1960,15,"0.844*BMSY",cex=1.5)
text(1960,12.5,"0.844*BMSY", cex=1.5, col="darkgray")

#####################Sensitivity of Harvest Rate to time varying catchability##################################
#src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_varq"
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\sens_varq"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

##grab the default biomass estimates
src.dir ="C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[70]"),1] ## <--- Change accordingly

dev.new(width=6,height=4)
plot(bio_1~seq(1949,2018,1),xlim=c(1945,2020), ylim=c(0,0.15), main="Sensitivity to inclusion of random-walk catchability",
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")

lines(bio_base~seq(1949,2018,1), ,lwd=4, col="darkgray")
legend(1950,0.14, c("Random-walk catchability"), 
       lty=c(1,1,1,1), pch=c(15), lwd=c(2), cex=1.75)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.2, by = 0.01), labels = FALSE, tcl = -0.3) 

abline(h=0.05613,lwd=2)
abline(h=0.0693 , lwd=2,col="darkgray")
text(1962,0.054,"HMSY",cex=1.5)
text(1962,0.066,"HMSY", cex=1.5, col="darkgray")

##NOT RUN IN 2020
#####################Sensitivity of Biomass to starting in 1980##################################
src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_1980"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[36]"),1] ## <--- Change accordingly


##grab the default biomass estimates
src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[32]"):which(rnames=="B[67]"),1] ## <--- Change accordingly

plot(bio_1~seq(1980,2015,1),xlim=c(1975,2020), ylim=c(0,30), main="Sensitivity to starting model in 1980",
     xlab="Year", ylab="Mean Biomass (million of lbs.)",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")

lines(bio_base~seq(1980,2015,1), ,lwd=4, col="darkgray")
legend(1995,29, c("Start in 1980"), 
       lty=c(1,1), pch=c(15), lwd=c(2), cex=1.3)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 30, by = 1), labels = FALSE, tcl = -0.3) 

abline(h=13.24*0.844,lwd=3)
#text(2014,15.5,"0.844*BMSY",cex=1.5)
abline(h=13.69*0.844 , lwd=3, col="darkgray")
text(2014,10.5,"0.844*BMSY", cex=1.5)


#####################Sensitivity of Harvest Rate to starting in 1980 ##################################
src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/sensWPSARresults/sens_1980"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[36]"),1] ## <--- Change accordingly

##grab the default biomass estimates
src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[32]"):which(rnames=="H[67]"),1] ## <--- Change accordingly


plot(bio_1~seq(1980,2015,1),xlim=c(1975,2020), ylim=c(0,0.15), main="Sensitivity to starting model in 1980",
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")

lines(bio_base~seq(1980,2015,1), ,lwd=4, col="darkgray")
legend(1995,0.14, c("Start in 1980"), 
       lty=c(1,1,1,1), pch=c(15), lwd=c(2), cex=1.3)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.2, by = 0.01), labels = FALSE, tcl = -0.3) 


abline(h=0.05667,lwd=3)
#text(2014,15.5,"0.844*BMSY",cex=1.5)
abline(h=0.0658 , lwd=3, col="darkgray")
text(2014,0.05,"HMSY", cex=1.5)


##########################Sensitivity to P1 Uniform##############################################################
#Sensitivity of biomass to P1 uniform 0,1 and 0,1.5##############################################################

src.dir ="C:/Users/John.Syslo/Documents/Deep 7 BSP/Base/log_sd/Sens_P1_uniform"
setwd(src.dir)

logtxt <-  bugs.log("log_0to1.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[67]"),1] ## <--- Change accordingly


logtxt <-  bugs.log("log_0to1.5.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="B[1]"):which(rnames=="B[67]"),1] ## <--- Change accordingly



##grab the default biomass estimates
src.dir ="C:/Users/John.Syslo/Documents/Deep 7 BSP/Base/log_sd/Final Base"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[67]"),1] ## <--- Change accordingly



plot(bio_1~seq(1949,2015,1),xlim=c(1945,2020), ylim=c(0,30), main="Sensitivity to P1 prior distribution",
     xlab="Year", ylab="Mean Biomass (million of lbs.)",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
lines(bio_2~seq(1949,2015,1), cex=1.5, type="o", pch=17,lwd=2)
lines(bio_base~seq(1949,2015,1),lwd=4, col="darkgray")


legend(1995,29, c("Uniform [0,1]","Uniform [0,1.5]"), 
       lty=c(1,1), pch=c(15,17), lwd=c(2), cex=2)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 30, by = 1), labels = FALSE, tcl = -0.3) 

abline(h=13.69*0.844,col="darkgray", lwd=2)
abline(h=12.72*0.844,lwd=2) #0,1
abline(h=12.29*0.844,lwd=2)  #0,1.5

text(1962,11,"0.844*BMSY Uniform [0,1]",cex=1.5)
text(1962,10,"0.844*BMSY Uniform [0,1.5]",cex=1.5)
text(1962,12,"0.844*BMSY Base", cex=1.5, col="darkgray")






#Sensitivity of harvest rate to P1 uniform 0,1 and 0,1.5##############################################################

src.dir ="C:/Users/John.Syslo/Documents/Deep 7 BSP/Base/log_sd/Sens_sv_downweight"
setwd(src.dir)

logtxt <-  bugs.log("log_0to1.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[67]"),1] ## <--- Change accordingly


logtxt <-  bugs.log("log_0to1.5.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_2 <- stats[which(rnames=="H[1]"):which(rnames=="H[67]"),1] ## <--- Change accordingly



##grab the default biomass estimates
src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[67]"),1] ## <--- Change accordingly



plot(bio_1~seq(1949,2015,1),xlim=c(1945,2020), ylim=c(0,0.15), main="Sensitivity to P1 prior distribution",
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.75, cex.sub=1.5,xaxs = "i",yaxs = "i")
lines(bio_2~seq(1949,2015,1), cex=1.5, type="o", pch=17,lwd=2)
lines(bio_base~seq(1949,2015,1),lwd=4, col="darkgray")


legend(1995,0.14, c("Uniform [0,1]","Uniform [0,1.5]"), 
       lty=c(1,1), pch=c(15,17), lwd=c(2), cex=2)  #col=c(colors, "black","black"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.15, by = 0.1), labels = FALSE, tcl = -0.3) 

abline(h=0.066,col="darkgray")
abline(h=0.06916) #0,1
abline(h=0.06912)  #0,1.5

text(1962,0.072,"HMSY Uniform",cex=1.5)
text(1962,0.063,"HMSY Base", cex=1.5, col="darkgray")

#####################Sensitivity of Biomass to survey downweighting##################################


src.dir ="C:/Users/John.Syslo/Documents/Deep 7 BSP/Base/log_sd/Sens_sv_downweight"
setwd(src.dir)
logtxt <-  bugs.log("log_sweight0.1.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="B[1]"):which(rnames=="B[67]"),1] ## <--- Change accordingly
bio_1_low<- stats[which(rnames=="B[1]"):which(rnames=="B[67]"),4] ## <--- Change accordingly
bio_1_hi<- stats[which(rnames=="B[1]"):which(rnames=="B[67]"),6] ## <--- Change accordingly
bio_1_sd<-stats[which(rnames=="B[1]"):which(rnames=="B[67]"),2] ## <--- Change accordingly

##grab the default biomass estimates
src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="B[1]"):which(rnames=="B[67]"),1] ## <--- Change accordingly
bio_base_low<- stats[which(rnames=="B[1]"):which(rnames=="B[67]"),4] ## <--- Change accordingly
bio_base_hi<- stats[which(rnames=="B[1]"):which(rnames=="B[67]"),6] ## <--- Change accordingly
bio_base_sd <- stats[which(rnames=="B[1]"):which(rnames=="B[67]"),2] ## <--- Change accordingly


plot(bio_1~seq(1949,2015,1),xlim=c(1945,2020), ylim=c(0,60), main="Sensitivity to downweighting fishery-independent survey",
     xlab="Year", ylab="Mean Biomass (million of lbs.)",
     type="o", col="red",cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xaxs = "i",yaxs = "i")
polygon(c(seq(1949,2015,1),rev(seq(1949,2015,1))),c(bio_1_low,rev(bio_1_hi)),col=rgb(1,0,0,0.2), border = "NA")
lines(bio_1~seq(1949,2015,1),type="o",col="red")
lines(seq(1949,2015,1),bio_1_low,lty=5,col=rgb(1,0,0,0.5),lwd=2)
lines(seq(1949,2015,1),bio_1_hi,lty=5,col=rgb(1,0,0,0.5),lwd=2)

lines(bio_base~seq(1949,2015,1),lwd=4,col="blue")
polygon(c(seq(1949,2015,1),rev(seq(1949,2015,1))),c(bio_base_low,rev(bio_base_hi)),col=rgb(0,0,1,0.2), border = "NA")
lines(bio_base~seq(1949,2015,1),lwd=4, col="blue")
lines(seq(1949,2015,1),bio_base_low,lty=3,col=rgb(0,0,1,0.5),lwd=2)
lines(seq(1949,2015,1),bio_base_hi,lty=3,col=rgb(0,0,1,0.5),lwd=2)

legend(1980,59, c("Survey downweighted - mean", "Survey downweighted - 95% CI", "Survey included - mean","Survey included - 95% CI"), 
       lty=c(1,5,1,3), pch=c(15,NA,NA,NA), lwd=c(2,2), cex=0.9,col=c("red","red","blue","blue"))


axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 60, by = 1), labels = FALSE, tcl = -0.3) 

abline(h=16.37*0.844,lwd=3, lty=5,col="red")
text(2010,15.5,"0.844*BMSY",cex=1.25,col="red")
abline(h=13.69*0.844 , lwd=3, lty=3, col="blue")
text(2010,10,"0.844*BMSY", cex=1.25, col="blue")


#calculate CI width differences to mention in text?
ci1<-bio_1_hi-bio_1_low
cibase<-bio_base_hi-bio_base_low

ci1/bio_1

cibase/bio_base


a<-bio_base_sd/bio_base

b<-bio_1_sd/bio_1


#####################Sensitivity of Harvest to survey downweighting##################################


src.dir =src.dir ="C:/Users/John.Syslo/Documents/Deep 7 BSP/Base/log_sd/Sens_sv_downweight"
setwd(src.dir)
logtxt <-  bugs.log("log_sweight0.1.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[67]"),1] ## <--- Change accordingly
bio_1_low<- stats[which(rnames=="H[1]"):which(rnames=="H[67]"),4] ## <--- Change accordingly
bio_1_hi<- stats[which(rnames=="H[1]"):which(rnames=="H[67]"),6] ## <--- Change accordingly


##grab
setwd(src.dir)
logtxt <-  bugs.log("log_sweight0.1.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_1 <- stats[which(rnames=="H[1]"):which(rnames=="H[67]"),1] ## <--- Change accordingly
bio_1_low<- stats[which(rnames=="H[1]"):which(rnames=="H[67]"),4] ## <--- Change accordingly
bio_1_hi<- stats[which(rnames=="H[1]"):which(rnames=="H[67]"),6] ## <--- Change accordingly


##grab the default biomass estimates
src.dir ="H:/Bottomfish/2017assessmentbenchmark/Complex BSP/d7_2020_base_NA_555"
setwd(src.dir)
logtxt <-  bugs.log("log.txt")
stats <- logtxt[["stats"]]
rnames <- rownames(stats)
bio_base <- stats[which(rnames=="H[1]"):which(rnames=="H[67]"),1] ## <--- Change accordingly
bio_base_low<- stats[which(rnames=="H[1]"):which(rnames=="H[67]"),4] ## <--- Change accordingly
bio_base_hi<- stats[which(rnames=="H[1]"):which(rnames=="H[67]"),6] ## <--- Change accordingly



plot(bio_1~seq(1949,2015,1),xlim=c(1945,2020), ylim=c(0,0.25), main="Sensitivity to downweighting of fishery-independent survey",
     xlab="Year", ylab="Mean Harvest Rate",
     type="o", cex=1.5, pch=15, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xaxs = "i",yaxs = "i",col="red")
polygon(c(seq(1949,2015,1),rev(seq(1949,2015,1))),c(bio_1_low,rev(bio_1_hi)),col=rgb(1,0,0,0.2), border = "NA")
lines(bio_1~seq(1949,2015,1),type="o",col="red")
lines(seq(1949,2015,1),bio_1_low,lty=5,col=rgb(1,0,0,0.5),lwd=2)
lines(seq(1949,2015,1),bio_1_hi,lty=5,col=rgb(1,0,0,0.5),lwd=2)

lines(bio_base~seq(1949,2015,1),lwd=4, col="blue")
polygon(c(seq(1949,2015,1),rev(seq(1949,2015,1))),c(bio_base_low,rev(bio_base_hi)),col=rgb(0,0,1,0.2), border = "NA")
lines(bio_base~seq(1949,2015,1),lwd=4, col="blue")
lines(seq(1949,2015,1),bio_base_low,lty=3,col=rgb(0,0,1,0.5),lwd=2)
lines(seq(1949,2015,1),bio_base_hi,lty=3,col=rgb(0,0,1,0.5),lwd=2)

legend(1950,0.24, c("Survey downweighted - mean", "Survey downweighted - 95% CI", "Survey included - mean","Survey included - 95% CI"), 
       lty=c(1,5,1,3), pch=c(15,NA,NA,NA), lwd=c(2,2), cex=1,col=c("red","red","blue","blue"))
axis(side = 1, at = seq(1945, 2020, by = 1), labels = FALSE, tcl = -0.3) 
axis(side = 2, at = seq(0, 0.25, by = 0.01), labels = FALSE, tcl = -0.3) 


abline(h=0.06614,lwd=3, lty=5,col="blue")
text(1962,0.07,"HMSY")
abline(h=0.0659 , lwd=3, lty=3, col="red")
