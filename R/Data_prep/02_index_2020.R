#Adapted from "02_index_brian.Rmd" by JOhn Syslo on 5-7-2020
#Not an Rmd file because of errors with my R markdown

  
  
#From STM:An index of relative abundance was calculated for each distribution. The index (Iy) and its variance were calculated 
#as the mean and variance, respectively, of the predicted values on the scale of the response in year y using the "predict" function
#in R. The variance of the delta-lognormal distribution was calculated as the Taylor series expansion of the variance of the product 
#of two independent random variables (Brodziak and Walsh, 2013; Eq. 7). Bias-correction was applied when back-transforming the 
#positive process of the DLN model from ln(CPUE) to CPUE. Coefficients of variation (CV) were also calculated for each year for all 
#distributions. Pearson correlations between each index pair were calculated to determine similarities among the indices.

  
rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)
library(lme4)
require(Rmisc)

#Follows: 01_data reformating.R
#         dataprepsource.R
#         best_mod_source.R

#Load the correct models

source('C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE standardization\\best_mod_source.R')

#NEW for 2020
#Need to fix TP2 areas as wind data before running the predict functions
for (i in 1:nrow(TP2)){
  TP2$speed[i]<-ifelse(TP2$speed[i]==TP2$area[i],TP2$speed[i]<-NA,TP2$speed[i]<-TP2$speed[i])
  TP2$xdir[i]<-ifelse(TP2$xdir[i]==TP2$area[i],TP2$xdir[i]<-NA,TP2$xdir[i]<-TP2$xdir[i])
  TP2$ydir[i]<-ifelse(TP2$ydir[i]==TP2$area[i],TP2$ydir[i]<-NA,TP2$ydir[i]<-TP2$ydir[i])
  TP2$winddeg[i]<-ifelse(TP2$speed[i]==TP2$area[i],TP2$winddeg[i]<-NA,TP2$winddeg[i]<-TP2$winddeg[i])
  TP2$winddir[i]<-ifelse(TP2$speed[i]==TP2$area[i],TP2$winddir[i]<-NA,TP2$winddir[i]<-TP2$winddir[i])
}

TP2$winddir = cut(TP2$winddeg, breaks = seq(0,360,45), labels = c("W","NW","N","NE","E","SE","S","SW"))
TP2<-TP2[,-10] #effort is line hours so can delete dayeffort
TP2 = TP2[complete.cases(TP2),] #was deleting everything because dayeffort is NA
dim(TP2) #49992 records

range(TP2$speed)


#Prior filtering did not actually account for 2016-2019 FI survey trips before running get_wind. Do that process again here
#-Remove fishing independent survey trips 
#Read in database of fisher, day, areas from the FI surveys
FI=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE\\supporting files\\FI to HDAR grid_ALLFY.csv",header=T)
FI$Date=as.Date(FI$Date,format="%m/%d/%Y")
FI$last=tolower(FI$last)
#Find FI survey trips in filtered database
FI$combo=paste(FI$Date,FI$HDAR.Grid,tolower(FI$Captain)) #do this to combine all unique conditions> If you cond1 & cond2.... then dont get unique conditions
a=TP2[paste(TP2$FISHED,TP2$area,TP2$fisher)%in%FI$combo,]
table(a$fisher) #insure that use of last name was appropriate and only got Captains from survey (it was)
#Check to ensure that a trip doesn't cover multiple areas (it doesn't)
dim(a[!duplicated(a[,c("FISHED","area","fisher")]),])
length(unique(a$trip3))
####----#####
TP2=TP2[!TP2$trip3%in%unique(a$trip3),] #13 records removed - lack of wind data during this time of year in recent data
dim(TP2)

#convert "null null" fisheres to NA
TP1[which(TP1$fisher == "null null"),5]<-NA 



#Code for comparing performance (deviance) among null models and models #with only year effects for all processes and time periods
TP1.B.0null=glm(bin.catch~1,data=TP1,family=binomial)
TP1.B.0year=glm(bin.catch~FYEAR,data=TP1,family=binomial)
(TP1.B.0null$deviance-TP1.B.best$deviance)/TP1.B.0null$deviance*100
(TP1.B.0year$deviance-TP1.B.best$deviance)/TP1.B.0year$deviance*100

TP2.B.0null=glm(bin.catch~1,data=TP2,family=binomial)
TP2.B.0year=glm(bin.catch~FYEAR,data=TP2,family=binomial)
(TP2.B.0null$deviance-TP2.B.best$deviance)/TP2.B.0null$deviance*100
(TP2.B.0year$deviance-TP2.B.best$deviance)/TP2.B.0year$deviance*100

TP1.RLN.test=lmer(log(cpue)~FYEAR+(1|fisher)+area+qtr+area:qtr+log_cum_exp+sqrt_uku_lbs,REML=F,data=TP1[TP1$cpue>0,])
testTP1=TP1[!is.na(TP1$fisher),] #remove NA fishers from TP1 data
TP1.RLN.0null=glm(log(cpue)~1,data=testTP1[testTP1$cpue>0,],na.action=na.omit) #need to have NA fishers removed
TP1.RLN.0year=glm(log(cpue)~FYEAR,data=testTP1[testTP1$cpue>0,]) #need to have NA fishers removed
TP1.RLN.0fish=lmer(log(cpue)~(1|fisher),REML=F,data=TP1[TP1$cpue>0,])


#Doing deviance is the same here when the model structure is the same, but need to do -2 * logLik when it is not
(deviance(TP1.RLN.0null)-deviance(TP1.RLN.test))/deviance(TP1.RLN.0null)

TP2.B.0null=glm(bin.catch~1,data=TP2,family=binomial)
TP2.B.0year=glm(bin.catch~FYEAR,data=TP2,family=binomial)
(TP2.B.0null$deviance-TP2.B.best$deviance)/TP2.B.0null$deviance*100
(TP2.B.0year$deviance-TP2.B.best$deviance)/TP2.B.0year$deviance*100


TP2.RLN.test=lmer(log(cpue)~FYEAR+(1|fisher)+area+qtr+area:FYEAR+log_cum_exp+sqrt_uku_lbs+speed,REML=F,data=TP2[TP2$cpue>0,])
TP2.RLN.0null=glm(log(cpue)~1,data=TP2[TP2$cpue>0,])
(-2*logLik(TP2.RLN.0null)--2*logLik(TP2.RLN.test))/(-2*logLik(TP2.RLN.0null))
TP2.RLN.0fishyear=lmer(log(cpue)~FYEAR+(1|fisher),REML=F,data=TP2[TP2$cpue>0,])
(-2*logLik(TP2.RLN.0fishyear)--2*logLik(TP2.RLN.test))/(-2*logLik(TP2.RLN.0fishyear))
TP2.RLN.0year=glm(log(cpue)~FYEAR,data=TP2[TP2$cpue>0,])
(-2*logLik(TP2.RLN.0year)--2*logLik(TP2.RLN.test))/(-2*logLik(TP2.RLN.0year))
TP2.RLN.0fish=lmer(log(cpue)~(1|fisher),REML=F,data=TP2[TP2$cpue>0,])
(-2*logLik(TP2.RLN.0fish)--2*logLik(TP2.RLN.test))/(-2*logLik(TP2.RLN.0fish))

##
#Determining marginal and conditional R2 for random effects models
#See website https://jonlefcheck.net/2013/03/13/r2-for-linear-mixed-effects-models/
#Marginal is proprotion of variance explained by fixed effects
#Conditional is proprotion of variance explained by fixed and random effects
###
library(piecewiseSEM)
sem.model.fits(TP1.B.best)
sem.model.fits(TP2.B.best)
sem.model.fits(TP1.RLN.best)
sem.model.fits(TP2.RLN.best)

#NEW for 2020
#Need to fix TP2 areas as wind data before running the predict functions




# Positive Process
## Extract predicted values for Positive process ("p") and bind to YEAR (for aggregating)

#TP1 is for first time period
#TP10 is for the first time period but with 1958-1960 data removed
#TP1W is for the first time period but covers only when wind data were available (1987-2003)
#TP2 is for second time period
TP1p = data.frame('LOGCPUE' = predict(TP1.RLN.best), 'FYEAR' = TP1[TP1$cpue>0,][,'FYEAR'])
#TP10p = data.frame('LOGCPUE' = predict(TP1.0.RLN.best), 'FYEAR' = TP1.0[TP1.0$cpue>0,][,'FYEAR'])
#TP1Wp = data.frame('LOGCPUE' = predict(TP1.W.RLN.best),'FYEAR' =TP1.W[TP1.W$cpue>0,][,'FYEAR'])
TP2p = data.frame('LOGCPUE' = predict(TP2.RLN.best),'FYEAR' =TP2[TP2$cpue>0,][,'FYEAR'])

## Backtransform positive process using dispersion from each model following Brian's STM standardization and that from Brodziak and Walsh (2013)

TP1p$trans=exp(TP1p$LOGCPUE + ((summary(TP1.RLN.best)$sigma^2)/2))
#TP10p$trans=exp(TP10p$LOGCPUE + ((summary(TP1.0.RLN.best)$sigma^2)/2))
#TP1Wp$trans=exp(TP1Wp$LOGCPUE + ((summary(TP1.W.RLN.best)$sigma^2)/2))
TP2p$trans=exp(TP2p$LOGCPUE + ((summary(TP2.RLN.best)$sigma^2)/2))

# Bernoulli process
## Extract predicted values for Bernoulli process ("b")- be sure to use type = 'response' which provides the probability of having a non zero tow (Stefansson 1996)

#TP1bnofisher = data.frame('BIN.CATCH' = predict(TP1.B.best_noNAfisher, type = 'response'),'FYEAR' = TP1.adj[,'FYEAR'])
TP1b = data.frame('BIN.CATCH' = predict(TP1.B.best, type = 'response'),'FYEAR' = TP1[,'FYEAR'])
#TP10b = data.frame('BIN.CATCH' = predict(TP1.0.B.best, type = 'response'),'FYEAR' = TP1.0[,'FYEAR'])
#TP1Wb = data.frame('BIN.CATCH' = predict(TP1.W.B.best, type = 'response'),'FYEAR' = TP1.W[,'FYEAR'])
TP2b = data.frame('BIN.CATCH' = predict(TP2.B.best, type = 'response'),'FYEAR' = TP2[,'FYEAR'])

## Use aggregate ('a') to get means and sd for each yearfor the positive process (remember sd2 is var)

TP1pa = aggregate(trans ~ FYEAR, TP1p, function(x) c(mean = mean(x), sd = sd(x), var = var(x)))
#TP10pa = aggregate(trans ~ FYEAR, TP10p, function(x) c(mean = mean(x), sd = sd(x), var = var(x)))
#TP1Wpa = aggregate(trans ~ FYEAR, TP1Wp, function(x) c(mean = mean(x), sd = sd(x), var = var(x)))
TP2pa = aggregate(trans ~ FYEAR, TP2p, function(x) c(mean = mean(x), sd = sd(x), var = var(x)))

## Use aggregate ('a') to get means, sd and variance for each year for the bernoulli process. Var for bernoulli is not standard. Don't transform after this.

bernvar = function(x) mean(x)*(1-mean(x))
#TP1banofisher = aggregate(BIN.CATCH ~ FYEAR, data = TP1bnofisher,  FUN = function(x) c(mean = mean(x), sd = sd(x), var = var(x),bernvar=bernvar(x)))
TP1ba = aggregate(BIN.CATCH ~ FYEAR, data = TP1b,  FUN = function(x) c(mean = mean(x), sd = sd(x), var = var(x),bernvar=bernvar(x)))
#TP10ba = aggregate(BIN.CATCH ~ FYEAR, data = TP10b,  FUN = function(x) c(mean = mean(x), sd = sd(x), var = var(x),bernvar=bernvar(x)))
#TP1Wba = aggregate(BIN.CATCH ~ FYEAR, data = TP1Wb,  FUN = function(x) c(mean = mean(x), sd = sd(x), var = var(x),bernvar=bernvar(x)))
TP2ba = aggregate(BIN.CATCH ~ FYEAR, data = TP2b,  FUN = function(x) c(mean = mean(x), sd = sd(x), var = var(x),bernvar=bernvar(x)))

# Index generation
## Multiply each estimate together and calculate the variance according to Brodziak and Walsh (2013) but ultimately following Goodman (1960)- no bother with the covariance as it is set to 0.

varI=function(pmean,pvar,bmean,bvar){
  index_totvar = bvar*pvar + bvar*(pmean^2) + pvar*(bmean^2)  }

#TP1Inofisher = data.frame('FYEAR' = TP1banofisher[,'FYEAR'], 'INDEX.EST' = TP1banofisher$BIN.CATCH[,'mean'] * TP1pa$trans[,'mean'],'VARIANCE.FORM' = varI(TP1pa$trans[,'mean'],TP1pa$trans[,'var'],TP1banofisher$BIN.CATCH[,'mean'],TP1banofisher$BIN.CATCH[,'var']),'VARIANCE.ADDITIVE'=TP1pa$trans[,'var']+TP1banofisher$BIN.CATCH[,'var'],'BERNVARIANCE.FORM' = varI(TP1pa$trans[,'mean'],TP1pa$trans[,'var'],TP1banofisher$BIN.CATCH[,'mean'],TP1banofisher$BIN.CATCH[,'bernvar']))
#TP1Inofisher$MODEL = 'TP1nofisher'

TP1I = data.frame('FYEAR' = TP1ba[,'FYEAR'], 'INDEX.EST' = TP1ba$BIN.CATCH[,'mean'] * TP1pa$trans[,'mean'],'VARIANCE.FORM' = varI(TP1pa$trans[,'mean'],TP1pa$trans[,'var'],TP1ba$BIN.CATCH[,'mean'],TP1ba$BIN.CATCH[,'var']),'VARIANCE.ADDITIVE'=TP1pa$trans[,'var']+TP1ba$BIN.CATCH[,'var'],'BERNVARIANCE.FORM' = varI(TP1pa$trans[,'mean'],TP1pa$trans[,'var'],TP1ba$BIN.CATCH[,'mean'],TP1ba$BIN.CATCH[,'bernvar']))
TP1I$MODEL = 'TP1'

#TP10I = data.frame('FYEAR' = TP10ba[,'FYEAR'], 'INDEX.EST' = TP10ba$BIN.CATCH[,'mean'] * TP10pa$trans[,'mean'], 'VARIANCE.FORM' = varI(TP10pa$trans[,'mean'],TP10pa$trans[,'var'],TP10ba$BIN.CATCH[,'mean'],TP10ba$BIN.CATCH[,'var']),'VARIANCE.ADDITIVE'=TP10pa$trans[,'var']+TP10ba$BIN.CATCH[,'var'],'BERNVARIANCE.FORM' = varI(TP10pa$trans[,'mean'],TP10pa$trans[,'var'],TP10ba$BIN.CATCH[,'mean'],TP10ba$BIN.CATCH[,'bernvar']))
#TP10I$MODEL = 'TP10'

#TP1WI = data.frame('FYEAR' = TP1Wba[,'FYEAR'], 'INDEX.EST' = TP1Wba$BIN.CATCH[,'mean'] * TP1Wpa$trans[,'mean'], 'VARIANCE.FORM' = varI(TP1Wpa$trans[,'mean'],TP1Wpa$trans[,'var'],TP1Wba$BIN.CATCH[,'mean'],TP1Wba$BIN.CATCH[,'var']),'VARIANCE.ADDITIVE'=TP1Wpa$trans[,'var']+TP1Wba$BIN.CATCH[,'var'], 'BERNVARIANCE.FORM' = varI(TP1Wpa$trans[,'mean'],TP1Wpa$trans[,'var'],TP1Wba$BIN.CATCH[,'mean'],TP1Wba$BIN.CATCH[,'bernvar']))
#TP1WI$MODEL = 'TP1W'

TP2I = data.frame('FYEAR' = TP2ba[,'FYEAR'], 'INDEX.EST' = TP2ba$BIN.CATCH[,'mean'] * TP2pa$trans[,'mean'], 'VARIANCE.FORM' = varI(TP2pa$trans[,'mean'],TP2pa$trans[,'var'],TP2ba$BIN.CATCH[,'mean'],TP2ba$BIN.CATCH[,'var']),'VARIANCE.ADDITIVE'=TP2pa$trans[,'var']+TP2ba$BIN.CATCH[,'var'],'BERNVARIANCE.FORM' = varI(TP2pa$trans[,'mean'],TP2pa$trans[,'var'],TP2ba$BIN.CATCH[,'mean'],TP2ba$BIN.CATCH[,'bernvar']))
TP2I$MODEL = 'TP2'


#set up full df

full.df = rbind(TP1I,TP2I)  #removed TP10I,TP1WI,
full.df$SD = sqrt(full.df$VARIANCE.FORM)
full.df$CV = full.df$SD/full.df$INDEX.EST
full.df$N=c(table(TP1$FYEAR)[1:56],table(TP2$FYEAR)[56:71])
full.df$SE=full.df$SD/sqrt(full.df$N)
full.df$CV_mean=full.df$SE/full.df$INDEX.EST
head(full.df)
write.csv(full.df,'C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE standardization\\Finalized_stdindex_0715_REMLF_NA.csv', row.names = F)

nom_cpue1 = aggregate(cpue ~ FYEAR, data = TP1, FUN = mean)
nom_cpue2 = aggregate(cpue ~ FYEAR, data = TP2, FUN = mean)

full.nom<-cbind(full.df,c(nom_cpue1[,2],nom_cpue2[,2]))
#names(full.nom)[9]<-"nominal"

#re-ran for nominal for Minling on 8/11/21
write.csv(full.nom,'C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE standardization\\Finalized_stdindex_081121_REMLF_w_nom.csv', row.names = F)


# Plotting
#Nominal CPUE for comparison

tempdf = df
tempdf$FYEAR = as.numeric(as.character(df$FYEAR))
tempdf2 = aggregate(cpue ~ FYEAR, data = tempdf, FUN = mean)
nom = ggplot(data = tempdf2, aes(x = FYEAR, y = cpue)) +
  theme_bw() +
  scale_x_continuous(breaks = seq(1948,2018,5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('nominal cpue, yearly averages') +
  geom_line()
nom

#Standardized

full.df$FYEAR = as.numeric(as.character(full.df$FYEAR))
TP1I$FYEAR = as.numeric(as.character(TP1I$FYEAR))
#TP10I$FYEAR = as.numeric(as.character(TP10I$FYEAR))
#TP1WI$FYEAR = as.numeric(as.character(TP1WI$FYEAR))
TP2I$FYEAR = as.numeric(as.character(TP2I$FYEAR))
#TP2Ifixed$FYEAR

ggplot(data = NULL, aes(x = FYEAR, y = INDEX.EST)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  
  scale_x_continuous(breaks = seq(1948,2018,5)) +
  scale_y_continuous(breaks = seq(0,max(full.df$INDEX.EST),25)) +
  # scale_color_brewer(palette = 'Accent') +
  
  ggtitle('Model Estimated CPUE, lbs') +
  ylab('Estimated CPUE Index, lbs') +
  
  ## add indices
  geom_line(data = TP1I, aes(col = 'TP1'), size = 1) +
#  geom_line(data = TP10I[TP10I$FYEAR %in% 1948:1957,], aes(col = 'TP10'), size = 1) +
 # geom_line(data = TP10I[TP10I$FYEAR %in% 1961:2003,], aes(col = 'TP10'), size = 1) +
#  geom_line(data = TP1WI, aes(col = 'TP1W'), size = 1) +
  geom_line(data = TP2I, aes(col = 'TP2'), size = 1) +
  
  #add data
  geom_point(data=tempdf2, aes(x=FYEAR,y=cpue)) +
  
  ## add SD
  geom_errorbar(data = full.df, aes(ymin = INDEX.EST - 2*SD,
                                    ymax = INDEX.EST + 2*SD))


#Zoom in figures

ggplot(data = NULL, aes(x = FYEAR, y = INDEX.EST)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank()) +
  
  scale_x_continuous(breaks = seq(1948,2003,5)) +
  scale_y_continuous(breaks = seq(0,max(full.df$INDEX.EST),10)) +
  # scale_color_brewer(palette = 'Accent') +
  
  ggtitle('Model Estimated CPUE 1948-2003, lbs') +
  ylab('Estimated CPUE Index, lbs') +
  
  ## add data
  geom_line(data = TP1I, aes(col = 'TP1'), size = 1) +
#  geom_line(data = TP1WI, aes(col = 'TP1W'), size = 1) +
#  geom_line(data = TP10I[TP10I$FYEAR %in% 1948:1957,], aes(col = 'TP10'), size = 1) +
#  geom_line(data = TP10I[TP10I$FYEAR %in% 1961:2003,], aes(col = 'TP10'), size = 1) +
   geom_line(data = TP2I, aes(col = 'TP2'), size = 1) +
  
  ## add SD
  geom_errorbar(data = full.df[full.df$FYEAR %in% 1980:2003,], aes(ymin = INDEX.EST - 2*SD,
                                                                   ymax = INDEX.EST + 2*SD))


multiplot(nom,est,zoom,cols = 3)


