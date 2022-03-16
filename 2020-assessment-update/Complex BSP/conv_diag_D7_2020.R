############################################################################
### Test for convergence of MCMC chains #############
### Uses 'coda' library and automates testing using "superdiag' function.
### Created by Annie Yau, annie.yau@noaa.gov, 2014.01.15
### Last modified 2014.05.06

rm(list=ls())
library(R2WinBUGS)  	# Load the R2WinBUGS library
library(coda)

addname <- 'd7_2020_base_NA_555'   ##<---------name of model----------------- ## Change accordingly
src.dir = paste("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\",addname,"\\",sep="")  ##<---##1Change accordingly
setwd(src.dir)


### Coda output files from WinBUGS provides different MCMC chains, which need 
### to be combined into a single mcmc list object for input to superdiag(). 
# read.coda.interactive() ## can also be used to create single mcmc list. 

### Read in each coda files as list and combine into one MCMC list object
coda1 <- read.coda("coda1.txt","codaIndex.txt")
coda2 <- read.coda("coda2.txt","codaIndex.txt")
coda3 <- read.coda("coda3.txt","codaIndex.txt")
codalist <- mcmc.list(coda1,coda2,coda3)

# ### Can also import results from log file, to create figures etc.
# logtxt <-  bugs.log(paste(src.dir,addname,"_log.txt",sep=""))
# stats <- logtxt[["stats"]]
# DICresults <- logtxt[["DIC"]]

#### Convergence tests ###########################################################

### List nodes to determine which node # you want to test
n = length(codalist[1,][[1]])
data.frame(id=c(1:n), var=c(names(codalist[1,][[1]])))
names <- as.character(varnames(codalist))
### Specify the node you want to test 
nodes <- which(names%in%c('BMSY', 'FMSY', 'HMSY', 'K', 'M', 'MSY', 'PMSY', 'P1', 'q1', 'q2', 'q3','r','rad', 'sigma2', 'tau2_1','tau2_2'))
names[nodes]

### Diagnostics!
gR_upCI=rep(NA,length(nodes))
hw_stat_pass=hw_stat_iter=hw_stat_pval=hw_half_pass=matrix(NA,3,length(nodes))
gw_z=matrix(NA,3,length(nodes))

sink(paste(addname,"_convergencediagnostics.txt",sep="")) ## save output as txt file

for (i in 1:length(nodes)) {
  writeLines(paste("**********DIAGNOSTICS FOR [",as.character(nodes[i]),"]",names[nodes[i]],"************"))
  writeLines(paste("Gelman, >>>>> 1?"))
  gR=gelman.diag(codalist[,nodes[i]],autoburnin=FALSE)
  print(gR)
  gR_upCI[i]=gR$psrf[2]
  
  writeLines(paste("Heidelberger & Welch, passed & p-value < 0.05?"))
  hw = heidel.diag(codalist[,nodes[i]])
  print(hw)
  hw_stat_pass[1,i]=hw[[1]][1]
  hw_stat_iter[1,i]=hw[[1]][2]
  hw_stat_pval[1,i]=hw[[1]][3]
  hw_half_pass[1,i]=hw[[1]][4]
  hw_stat_pass[2,i]=hw[[2]][1]
  hw_stat_iter[2,i]=hw[[2]][2]
  hw_stat_pval[2,i]=hw[[2]][3]
  hw_half_pass[2,i]=hw[[2]][4]
  hw_stat_pass[3,i]=hw[[3]][1]
  hw_stat_iter[3,i]=hw[[3]][2]
  hw_stat_pval[3,i]=hw[[3]][3]
  hw_half_pass[3,i]=hw[[3]][4]
  
  writeLines(paste("Geweke, absolute value < 2?"))
  gw = geweke.diag(codalist[,nodes[i]])
  print(gw)
  gw_z[1,i]=gw[[1]]$z
  gw_z[2,i]=gw[[2]]$z 
  gw_z[3,i]=gw[[3]]$z
  
  writeLines(paste("Autocorrelation"))
  print(paste("Node",i,"chain 1"))
  print(autocorr.diag(codalist[[1]][,nodes[i]]))
  print(paste("Node",i,"chain 2"))
  print(autocorr.diag(codalist[[2]][,nodes[i]]))
  print(paste("Node",i,"chain 3"))
  print(autocorr.diag(codalist[[3]][,nodes[i]]))
}
print("GR")
gR_upCI
print("HW iter")
hw_stat_iter
print("HW pass")
hw_stat_pass
print("HW half pass")
hw_half_pass
print("Geweke")
gw_z

sink()

hw_stat_iter
gw_z


#Plots of trace
names[nodes]
plot(coda1[,"tau2_2"],trace=T,density=T,xlab="tau2_2",main="",ylab="Posterior density") 
mean(coda1[,"tau2_2"])
plot(coda2[,"tau2_2"],trace=T,density=T,xlab="tau2_2",main="",ylab="Posterior density") 
mean(coda2[,"tau2_2"])
plot(coda3[,"tau2_2"],trace=T,density=T,xlab="tau2_2",main="",ylab="Posterior density") 
mean(coda3[,"tau2_2"])
plot(coda2[,"sigma2"],trace=T,density=T,xlab="tau2_2",main="",ylab="Posterior density") 
mean(coda2[,"sigma2_1"])


##For creating the table in the assessment report
adjust=c(1,2,3,6,7,10,4,5,8,9,11:15)
names[nodes][adjust]
hw.p=NA
for (i in 1:length(nodes)) {
  hw = heidel.diag(codalist[,nodes[i]])
  hw.p[i]=min(c(hw[[1]][3],hw[[2]][3],hw[[3]][3])) #minimum pvalue across chains
}
ac1=autocorr.diag(codalist[[1]][,nodes])[2:3,]
ac2=autocorr.diag(codalist[[2]][,nodes])[2:3,]
ac3=autocorr.diag(codalist[[3]][,nodes])[2:3,]
acall=pmax(ac1,ac2,ac3)

cbind(gR_upCI,"gw.z"=apply(abs(gw_z),2,max),hw.p,"lag1"=acall[1,],"lag2"=acall[2,])[adjust,]
