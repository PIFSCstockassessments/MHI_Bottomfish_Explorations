

#MOUSS camera lengths
cam_length<-read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH 2016-2019 Camera Lengths.csv",header = T)
head(cam_length)


paka_dat<-cam_length[cam_length$SPECIES_CD=="PRFI",]
hist(paka_dat$LENGTH_CM)

tapply(paka_dat$LENGTH_CM,paka_dat$BFISH,mean)
tapply(paka_dat$LENGTH_CM,paka_dat$BFISH,median)
tapply(paka_dat$LENGTH_CM,paka_dat$BFISH,range)
tapply(paka_dat$LENGTH_CM,paka_dat$BFISH,sd)
tapply(paka_dat$LENGTH_CM,paka_dat$BFISH,length)


#research fishing lengths
res_length<-read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH 2016-2019 Research Fishing Lengths.csv",header = T)
head(res_length)


paka_dat_res<-res_length[res_length$SPECIES_CD=="PRFI",]
hist(paka_dat_res$LENGTH_CM)

tapply(paka_dat_res$LENGTH_CM,paka_dat_res$BFISH,mean,na.rm=T)
tapply(paka_dat_res$LENGTH_CM,paka_dat_res$BFISH,median,na.rm=T)
tapply(paka_dat_res$LENGTH_CM,paka_dat_res$BFISH,range,na.rm=T)
tapply(paka_dat_res$LENGTH_CM,paka_dat_res$BFISH,sd,na.rm=T)
tapply(paka_dat_res$LENGTH_CM,paka_dat_res$BFISH,length)

plot(density(paka_dat$LENGTH_CM))

x<-na.omit(paka_dat_res$LENGTH_CM)
lines(density(x))

#overlay histograms
hist(paka_dat$LENGTH_CM, col=rgb(1,0,0,0.5), main="Length frequency BFISH")
hist(paka_dat_res$LENGTH_CM, col=rgb(0,0,1,0.5), add=T)

#####look at for each survey by year - Ignore spring surveys 
length(unique(paka_dat$BFISH)) #4 years, with 5 surveys (2016/17 - spring/fall)


paka_dat16f<-paka_dat[paka_dat$BFISH=="BFISH_2016_F",]
paka_dat17s<-paka_dat[paka_dat$BFISH=="BFISH_2017_S",]
paka_dat17f<-paka_dat[paka_dat$BFISH=="BFISH_2017_F",]
paka_dat18f<-paka_dat[paka_dat$BFISH=="BFISH_2018_F",]
paka_dat19f<-paka_dat[paka_dat$BFISH=="BFISH_2019_F",]

paka_dat_res16f<-paka_dat_res[paka_dat_res$BFISH=="BFISH_2016_F",]
paka_dat_res16s<-paka_dat_res[paka_dat_res$BFISH=="BFISH_2016_S",]
paka_dat_res17f<-paka_dat_res[paka_dat_res$BFISH=="BFISH_2017_F",]
paka_dat_res18f<-paka_dat_res[paka_dat_res$BFISH=="BFISH_2018_F",]
paka_dat_res19f<-paka_dat_res[paka_dat_res$BFISH=="BFISH_2019_F",]


par(mfrow=c(2,2))

brvec<-seq(0,80,5)

hist(paka_dat16f$LENGTH_CM, col=rgb(1,0,0,0.5), ylim=c(0,30),breaks=brvec,main="Length frequency BFISH")
hist(paka_dat_res16f$LENGTH_CM, col=rgb(0,0,1,0.5), ylim=c(0,30),breaks=brvec, add=T)

hist(paka_dat17f$LENGTH_CM, col=rgb(1,0,0,0.5),ylim=c(0,30),breaks=brvec, main="Length frequency BFISH")
hist(paka_dat_res17f$LENGTH_CM, col=rgb(0,0,1,0.5),ylim=c(0,30),breaks=brvec, add=T)

hist(paka_dat18f$LENGTH_CM, col=rgb(1,0,0,0.5),ylim=c(0,30),breaks=brvec, main="Length frequency BFISH")
hist(paka_dat_res18f$LENGTH_CM, col=rgb(0,0,1,0.5),ylim=c(0,30),breaks=brvec, add=T)

hist(paka_dat19f$LENGTH_CM, col=rgb(1,0,0,0.5),ylim=c(0,30),breaks=brvec, main="Length frequency BFISH")
hist(paka_dat_res19f$LENGTH_CM, col=rgb(0,0,1,0.5),ylim=c(0,30),breaks=brvec, add=T)



#pool all bfish lengths and make single hist
all_bfish_len<-c(paka_dat$LENGTH_CM, paka_dat_res$LENGTH_CM)

hist(all_bfish_len)





#look at UFA lengths
#read in population length-frequency sample
ufa_lengths<-read.csv("C:\\Users\\John.Syslo\\Documents\\Opakapaka age comps\\UFA_length_dist.csv",header = T)
hist(ufa_lengths$FL, breaks=seq(250,850,25), main= NULL,xlab="Fork Length (mm)")

mean(ufa_lengths$FL)
median(ufa_lengths$FL)
range(ufa_lengths$FL)
sd(ufa_lengths$FL)
length(ufa_lengths$FL)

#BFISH camera vs UFA
y<-ufa_lengths$FL/10 #rescale to cm
brvec<-seq(0,100,5)
hist(paka_dat$LENGTH_CM, col=rgb(1,0,0,0.5), main="Length frequency BFISH Camera vs UFA",freq=FALSE,ylim=c(0,0.05),breaks=brvec)
hist(y, col=rgb(0,0,1,0.5),freq=FALSE, breaks=brvec,add=T)  #many more - need rel hist or density

legend(60,0.049,legend=c("BFISH Camera", "UFA"),
       col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), lty=1:2,lwd=5, cex=1.5)


#Bfish camera vs BFISH res fishing vs UFA
#BFISH camera vs UFA
y<-ufa_lengths$FL/10 #rescale to cm
brvec<-seq(0,100,5)
hist(paka_dat$LENGTH_CM, col=rgb(1,0,0,0.5), main="Length frequency BFISH Camera vs UFA",freq=FALSE,ylim=c(0,0.05),breaks=brvec)
hist(paka_dat_res$LENGTH_CM, col=rgb(0,1,0,0.5),freq=FALSE,ylim=c(0,0.05),breaks=brvec,add=T)
hist(y, col=rgb(0,0,1,0.5),freq=FALSE, breaks=brvec,add=T)  #many more - need rel hist or density


legend(60,0.049,legend=c("BFISH Camera", "BFISH fishing","UFA"),
       col=c(rgb(1,0,0,0.5),col=rgb(0,1,0,0.5), rgb(0,0,1,0.5)), lty=1:2,lwd=5, cex=1.5)


#BFISH camera vs bfish res fishing
brvec<-seq(0,100,5)
hist(paka_dat$LENGTH_CM, col=rgb(1,0,0,0.5), main="Length frequency BFISH Camera vs BFISH fishing",freq=FALSE,ylim=c(0,0.05),breaks=brvec)
hist(paka_dat_res$LENGTH_CM, col=rgb(0,1,0,0.5),freq=FALSE,ylim=c(0,0.05),breaks=brvec,add=T)

legend(60,0.049,legend=c("BFISH Camera", "BFISH fishing"),
       col=c(rgb(1,0,0,0.5),col=rgb(0,1,0,0.5)), lty=1:2,lwd=5, cex=1.5)

#UFA vs BFISH res fishing
y<-ufa_lengths$FL/10 #rescale to cm
brvec<-seq(0,100,5)
hist(paka_dat_res$LENGTH_CM, col=rgb(0,1,0,0.5),main="Length frequency BFISH Fishing vs UFA",freq=FALSE,ylim=c(0,0.05),breaks=brvec)
hist(y, col=rgb(0,0,1,0.5),freq=FALSE, breaks=brvec,add=T)  #many more - need rel hist or density


legend(60,0.049,legend=c("BFISH fishing","UFA"),
       col=c(rgb(0,1,0,0.5), rgb(0,0,1,0.5)), lty=1:2,lwd=5, cex=1.5)


##all three density plot
plot(density(paka_dat$LENGTH_CM),ylim=c(0,0.05),col=rgb(1,0,0,0.5),lwd=2)
x<-na.omit(paka_dat_res$LENGTH_CM) #BFISH fishing - No NAs
lines(density(x),col=rgb(0,1,0,0.5),lwd=2)
lines(density(y), col=rgb(0,0,1,0.5),lwd=2) #UFA rescaled

legend(-10,0.049,legend=c("BFISH Camera", "BFISH fishing","UFA"),
       col=c(rgb(1,0,0,0.5),col=rgb(0,1,0,0.5), rgb(0,0,1,0.5)), lty=1:2,lwd=5, cex=1.5)

