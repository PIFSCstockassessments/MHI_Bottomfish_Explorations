
rm(list=ls())

 ##<---------name of model----------------- ## Change accordingly
addname <- 'd7_2020_base_PROJ_NA'
src.dir<-paste("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\",addname,"\\",sep="")

setwd(src.dir)
library(R2WinBUGS)  	# Load the R2WinBUGS library
library(coda)

### Coda output files from WinBUGS provides different MCMC chains, which need 
### to be combined into a single mcmc list object for input to superdiag(). 
# read.coda.interactive() ## can also be used to create single mcmc list. 

### Read in each coda files as list and combine into one MCMC list object
coda1 <- read.coda("coda1.txt", "codaIndex.txt")
coda2 <- read.coda("coda2.txt", "codaIndex.txt")
coda3 <- read.coda("coda1.txt", "codaIndex.txt")
codalist <- mcmc.list(coda1, coda2, coda3)



head(codalist)

x<-summary(codalist)
stats<-(x$statistics)  #first element here is the mean, median will be used for estimate of harvest rate

#only 4 years of B because 2020 not affected by simulated TAC, nothing labeled as B2
x1<-stats[which(rownames(stats)=="proj_pOFL_B[1]"):which(rownames(stats)=="proj_pOFL_B[501]"),1]
x2<-stats[which(rownames(stats)=="proj_pOFL_B3[1]"):which(rownames(stats)=="proj_pOFL_B3[501]"),1]
x3<-stats[which(rownames(stats)=="proj_pOFL_B4[1]"):which(rownames(stats)=="proj_pOFL_B4[501]"),1]
x4<-stats[which(rownames(stats)=="proj_pOFL_B5[1]"):which(rownames(stats)=="proj_pOFL_B5[501]"),1]

plot(x1~seq(0,1000,2),type="l",ylim=c(0,1),xlab= "Reported Catch (1,000 lb)", ylab = "Probability of being overfished")
lines(x2~seq(0,1000,2),lty=2)
lines(x3~seq(0,1000,2),lty=3)
lines(x4~seq(0,1000,2),lty=4)


legend(100,0.9,c("2022","2023","2024", "2025"),lty=c(1,2,3,4))



x6<-stats[which(rownames(stats)=="proj_pOFL_H1[1]"):which(rownames(stats)=="proj_pOFL_H1[501]"),1]
x7<-stats[which(rownames(stats)=="proj_pOFL_H2[1]"):which(rownames(stats)=="proj_pOFL_H2[501]"),1]
x8<-stats[which(rownames(stats)=="proj_pOFL_H3[1]"):which(rownames(stats)=="proj_pOFL_H3[501]"),1]
x9<-stats[which(rownames(stats)=="proj_pOFL_H4[1]"):which(rownames(stats)=="proj_pOFL_H4[501]"),1]
x10<-stats[which(rownames(stats)=="proj_pOFL_H5[1]"):which(rownames(stats)=="proj_pOFL_H5[501]"),1]


plot(x6~seq(0,1000,2),type="l",ylim=c(0,1),xlab= "Reported Catch (1,000 lb)", ylab = "Probability of overfishing")
lines(x7~seq(0,1000,2),lty=2)
lines(x8~seq(0,1000,2),lty=3)
lines(x9~seq(0,1000,2),lty=4)
lines(x10~seq(0,1000,2),lty=5)

legend(100,0.9,c("2021","2022","2023","2024", "2025"),lty=c(1,2,3,4,5))

#what catch corresponds to 0.5 for each year?
tab<-cbind(seq(0,1000,2),x6,x7,x8,x9,x10)


tab[200:300,] #looks like 556 for 2025

cbind(seq(0,1000,2),x10) #looks like 556 for 2025

#2017 report was 558 for 2022

###monitored absolute values of B and H in separate run, so change directory for this

#rm(list=ls())

#addname <- 'd7_2020_base_PROJECT'   ##<---------name of model----------------- ## Change accordingly
addname <- 'd7_2020_base_PROJ_NA_B_H'
#src.dir = paste("X:/AYau/Bottomfish/Deep7Assess_2014/",addname,"/",sep="")  ##<---##1Change accordingly
#src.dir<-paste('C:\\Users\\John.Syslo\\Documents\\Deep 7 BSP\\Base\\log_sd')
src.dir<-paste("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\",addname,"\\",sep="")

setwd(src.dir)
library(R2WinBUGS)  	# Load the R2WinBUGS library
library(coda)

### Coda output files from WinBUGS provides different MCMC chains, which need 
### to be combined into a single mcmc list object for input to superdiag(). 
# read.coda.interactive() ## can also be used to create single mcmc list. 

### Read in each coda files as list and combine into one MCMC list object
coda1 <- read.coda("coda1.txt", "codaIndex.txt")
coda2 <- read.coda("coda2.txt", "codaIndex.txt")
coda3 <- read.coda("coda1.txt", "codaIndex.txt")
codalist <- mcmc.list(coda1, coda2, coda3)



head(codalist)

y<-summary(codalist)
stats2<-(y$statistics)  #first element here is the mean
quants<-(y$quantiles)

#only 4 years of B because 2020 not affected by simulated TAC, nothing labeled as B2
x11<-stats2[which(rownames(stats2)=="proj_B2022[1]"):which(rownames(stats2)=="proj_B2022[501]"),1]
x12<-stats2[which(rownames(stats2)=="proj_B2023[1]"):which(rownames(stats2)=="proj_B2023[501]"),1]
x13<-stats2[which(rownames(stats2)=="proj_B2024[1]"):which(rownames(stats2)=="proj_B2024[501]"),1]
x14<-stats2[which(rownames(stats2)=="proj_B2025[1]"):which(rownames(stats2)=="proj_B2025[501]"),1]
plot(x11~seq(0,1000,2),type="l",ylim=c(10,24),xlab= "Reported Catch (1,000 lb)",ylab="Mean exploitable biomass (million lb)",lty=1)
lines(x12~seq(0,1000,2),lty=2)
lines(x13~seq(0,1000,2),lty=3)
lines(x14~seq(0,1000,2),lty=4)
legend(100,16,c("2022","2023","2024", "2025"),lty=c(1,2,3,4))


#this needs to be median and not mean of posteriors
x16<-quants[which(rownames(quants)=="proj_H1[1]"):which(rownames(quants)=="proj_H1[501]"),3]
x17<-quants[which(rownames(quants)=="proj_H2[1]"):which(rownames(quants)=="proj_H2[501]"),3]
x18<-quants[which(rownames(quants)=="proj_H3[1]"):which(rownames(quants)=="proj_H3[501]"),3]
x19<-quants[which(rownames(quants)=="proj_H4[1]"):which(rownames(quants)=="proj_H4[501]"),3]
x20<-quants[which(rownames(quants)=="proj_H5[1]"):which(rownames(quants)=="proj_H5[501]"),3]

plot(x16~seq(0,1000,2),type="l",ylim=c(0,0.2),xlab= "Reported Catch (1,000 lb)",ylab="Harvest rate",lty=1)
lines(x17~seq(0,1000,2),lty=2)
lines(x18~seq(0,1000,2),lty=3)
lines(x19~seq(0,1000,2),lty=4)
lines(x20~seq(0,1000,2),lty=5)
legend(100,0.18,c("2021","2022","2023","2024", "2025"),lty=c(1,2,3,4,5))




#make big table

big_tab<-cbind(seq(0,1000,2),x1,x2,x3,x4,x6,x7,x8,x9,x10,x11,x12,x13,x14,x16,x17,x18,x19,x20)
big_tab<-as.data.frame(big_tab)
#Tables 15 and 16
#for each probability of overfishing from 0 to 0.5, find corresponding catch, pOFLB, Harvest rate, and Biomass

#Table 15

probs<-seq(0,0.5,0.05)

new_tab<-matrix(NA,18,(length(probs)))

for (i in 1:length(probs)){
#year 1 - biomass metrics not of interest
  maxless1 <- with(big_tab,   max(x6[x6 <= probs[i]]))
  idx <- which(big_tab$x6 == maxless1)
  catch<-big_tab[idx,1]
  Harvest<-big_tab[idx,15 ]

  new_tab[1,i]<-catch
  new_tab[10,i]<-Harvest
  
#year 2 
  maxless2 <- with(big_tab,   max(x7[x7 <= probs[i]]))
  idx <- which(big_tab$x7 == maxless2)
  catch<-big_tab[idx,1 ]
  prob_B<-big_tab[idx,2 ]
  Harvest<-big_tab[idx,16 ]
  Biomass<-big_tab[idx,11 ]
  
  new_tab[2,i]<-catch
  new_tab[11,i]<-Harvest
  new_tab[6,i]<-prob_B
  new_tab[15,i]<-Biomass
  
#year 3 
  maxless3 <- with(big_tab,   max(x8[x8 <= probs[i]]))
  idx <- which(big_tab$x8 == maxless3)
  catch<-big_tab[idx,1 ]
  prob_B<-big_tab[idx,3 ]
  Harvest<-big_tab[idx,17 ]
  Biomass<-big_tab[idx,12 ] 
  
  new_tab[3,i]<-catch
  new_tab[12,i]<-Harvest
  new_tab[7,i]<-prob_B
  new_tab[16,i]<-Biomass
  
#year 4
  maxless4 <- with(big_tab,   max(x9[x9 <= probs[i]]))
  idx <- which(big_tab$x9 == maxless4)
  catch<-big_tab[idx,1 ]
  prob_B<-big_tab[idx,4 ]
  Harvest<-big_tab[idx,18 ]
  Biomass<-big_tab[idx,13 ]  
 
  new_tab[4,i]<-catch
  new_tab[13,i]<-Harvest
  new_tab[8,i]<-prob_B
  new_tab[17,i]<-Biomass
#year 5 
  maxless5 <- with(big_tab,   max(x10[x10 <= probs[i]]))
  idx <- which(big_tab$x10 == maxless5)
  catch<-big_tab[ idx,1 ]
  prob_B<-big_tab[idx,5 ]
  Harvest<-big_tab[idx,19 ]
  Biomass<-big_tab[idx,14 ]  
   
  new_tab[5,i]<-catch
  new_tab[14,i]<-Harvest
  new_tab[9,i]<-prob_B
  new_tab[18,i]<-Biomass
   }

write.csv(new_tab, "C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\proj_table15.csv")


#Table 16
#just the catches for each prob 0-0.5
probs<-seq(0,0.5,0.01)
new_tab<-matrix(NA,length(probs),5)
for (i in 1:length(probs)){
    #year1
    maxless1 <- with(big_tab,   max(x6[x6 <= probs[i]]))
    idx <- max(which(big_tab$x6 == maxless1)) #two similar values occurred
    catch<-big_tab[idx,1]
    new_tab[i,1]<-catch
    
    #year2
    maxless2 <- with(big_tab,   max(x7[x7 <= probs[i]]))
    idx <- max(which(big_tab$x7 == maxless2))
    catch<-big_tab[idx,1]    
    new_tab[i,2]<-catch
    
    #year3
    maxless3 <- with(big_tab,   max(x8[x8 <= probs[i]]))
    idx <- max(which(big_tab$x8 == maxless3))
    catch<-big_tab[idx,1]   
    new_tab[i,3]<-catch
    
    #year4
    maxless4 <- with(big_tab,   max(x9[x9 <= probs[i]]))
    idx <- max(which(big_tab$x9 == maxless4))
    catch<-big_tab[idx,1]    
    new_tab[i,4]<-catch
    
    #year5
    maxless5 <- with(big_tab,   max(x10[x10 <= probs[i]]))
    idx <- max(which(big_tab$x10 == maxless5))
    catch<-big_tab[idx,1]
    new_tab[i,5]<-catch
}


write.csv(new_tab, "C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\proj_table16.csv")



#look at some values across chains - 2019 projection
coda1[1,2579:3079]

x<-colMeans(coda1[,2579:3079]) #proj_Pofl5
y<-colMeans(coda2[,2579:3079])
z<-colMeans(coda3[,2579:3079])
cbind(x,y,z)
v<-seq(0,1000,2)

par(mfrow=c(2,1))

plot(x~v, xlab= "ACL (thousands of lbs.)", ylab = "P(H>HMSY)", main = "2019")
points(y~v, pch = 0, col = 2)
points(z~v, pch = 4, col = 3)

v<-seq(0,1000,5)

abline(h=0.5)

#look at some values across chains - 2022 projection
x<-colMeans(coda1[,4082:4582]) #proj_Pofl5
y<-colMeans(coda2[,4082:4582])
z<-colMeans(coda3[,4082:4582])
cbind(x,y,z)


v<-seq(0,1000,5)

plot(coda1[,3015:3214])

plot(x~v, xlab= "ACL (thousands of lbs.)", ylab = "P(H>HMSY)", main = "2022")
points(y~v, pch = 0, col=2)
points(z~v, pch = 4, col=3)


abline(h=0.5)

#what is biomass before calculating pofl?
coda1[1,1:70]

x<-colMeans(coda1[,1:70])
y<-colMeans(coda2[,1:70])
z<-colMeans(coda3[,1:70])

v<-seq(1949,2018,1)

plot(x~v, xlab="Year", ylab = "Biomass")
points(y~v, pch = 0, col=2)
points(z~v, pch = 4, col=3)




