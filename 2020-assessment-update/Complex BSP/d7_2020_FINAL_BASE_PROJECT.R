##############################################################
# d7_2020_baseFINAL
# Jon Brodziak, PIFSC, December 2010, updated by Annie Yau, May 2014
# to two-CPUE time series. Updated further by Brian Langseth, April 2017, for 2018 benchmark 
#John Syslo updated June 2020 for 2021 update assessment

# Catch is in million pounds
# CPUE is in lbs/trip up before 10/1/2002, and lbs/hr thereafter

# Time period for two-CPUE indices, 1948-2002 and 2002-2018 (calendar year)
# and so use revised data entry structure.
# The CVs for years where CPUE is not used must still be entered, so
# that the code runs proporly. 

# Single catchability value per index
# Include fitting to survey biomass, with sd of survey on scale of log of data
# Includes capacity to set a weight to the survey std dev to fit the survey
# better (downeight survey sd)
# Exclude 1948 and use actual 2019 catch to set catch for that year
# Use natural mortality of 0.156

# Updated the survey to reflect a prior around the survey catchability
# based on min and max effective radiuas, corresponding to min and max
# scalar of the relative estimate of 46-1415, centered at 195

# Updated November 14, 2017
#############################################################

rm(list=ls())
#DATA = read.csv("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\model code\\SPM code\\d7_2017_data_REMLF.csv",header=T)
DATA = read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\data_NA.csv",header=T)
#DATA = read.csv("C:\\Users\\John.Syslo\\Documents\\2020 Deep 7\\Complex BSP\\data_early_old.csv",header=T)
head(DATA)
DATA=DATA[-1,]

Survey_data = read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Survey_data.csv",header=T)

addname <- 'd7_2020_base_PROJ_NA_B_H'  ##<--------name of model---------- # change accordingly
#src.dir <- paste('D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\model code\\SPM code\\',addname,"\\",sep="") # Change accordingly
src.dir <- paste('C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\',addname,"\\",sep="")

dir.create(src.dir)
dest.dir <- src.dir # where you want files copied to
setwd(src.dir)

library(R2WinBUGS)		# Load the R2WinBUGS library
library(coda)


##########################################################################
#These represent the previous measures:
#TAC_2016 <- 0.326
#prop_caught_2016 <- 0.802 # actual was 261,329
#
RC_2019 = 0.177051     #reported catch from 2016 (millions of lbs.) .  Unreported ratio is 1.07, so similar to 1.08 used in previous assessment's projections
#9/7/17 - run with 1.06

#TAC_2017 <- 0.318  
#prop_caught_2017 <- 0.7 #0.688  # equals 310,000 lbs as of 6/30/17

RC_2020 = 0.218388    #2017 catch is mean of last 3 years empirical data - should this be mean of ESTIMATED total catch??
##########################################################################






nt <- 20     # Thinning rate
ni <- 555000 # Number of total iterations per chain, including burn-in = multiply # of iterations desired by nt. 130,000 runs per chain thinned by 5 starting at 10,000
nb <- 255000 #round(ni*(1/10)) # Number of draws to discard as burn in; discard 10,000 


######################################################################
# DATA  
# model variable set-up
######################################################################
###obs_CPUE_1 = na.rm(DATA$CPUE_1_1)
# In this case, there is one CPUE set split at 1994 into two
# Vector Catch() is total catch weight in thousand metric tons 1949-2015
# Vector S1() is the Main Hawaiian Islands CPUE index 1949-1993
# Vector S1() is the Main Hawaiian Islands CPUE index 1993-2015

# sigma2 is process error
# tau2 is observation error by survey

NTIME <- length(DATA$Catch) 
Reported_Catch <- DATA$Catch  

UnrepCatch <- DATA$UnrepCatch

#CPUE and relCV of CPUE
CPUE_S1 <- DATA$CPUE_1
CPUE_S2 <- DATA$CPUE_2
CPUE_S1_REL_CV <- DATA$CPUE_1_rel_CV[!is.na(DATA$CPUE_1_rel_CV)] #exclude NAs
CPUE_S2_REL_CV <- DATA$CPUE_2_rel_CV[!is.na(DATA$CPUE_2_rel_CV)] #exclude NAs

#Biomass and SE of survey  #will probably need to remove last 2 years from this data frame
BIOSV<- Survey_data$Biomass_kg
SE_SV<- Survey_data$SE_Biomass_kg


#Accounting of time series length and dealing with NAs
NCPUE_S1_1=0
NCPUE_S1_MISS=0
NCPUE_S1_2 <- max(which(!is.na(DATA$CPUE_1))) #end year of first time series
if (match(NA, CPUE_S1)>0 & match(NA,CPUE_S1)!=(NCPUE_S1_2+1)){ #if there is an NA in first time period, prior to when the first time period ends
  NCPUE_S1_1 <- match(NA, CPUE_S1)-1 #last year prior to first NA
  NCPUE_S1_MISS <- length(DATA$CPUE_1[is.na(DATA$CPUE_1)]) + max(which(!is.na(DATA$CPUE_1)))-length(CPUE_S1) # Total missing values within time series (last positive + total NAs - total length)
}
NCPUE_S1_2 <- max(which(!is.na(DATA$CPUE_1))) #end year of first time series
NCPUE_S1_3 <- length(DATA$CPUE_1) #end year of all time series

NSV_1<- length(DATA$CPUE_1)-1 #account for 4 years of survey data: 2016 - 2019 (Fishing years 2017-2020)
NSV_2<- length(DATA$CPUE_1)   #17 and 18 will be within a loop 
  
#set up other survey index bounds to avoid errors within BUGS
NSVa<-NSV_1-1
NSVb<-NSV_2+1
NSVc<-((NSV_2+1)-(NSV_1-1))
NSVd<-NSV_2+2
NSVe<-((NSV_2+2)-(NSV_1-1))

#Survey biomass and SE estimate for 2016-19 calendar years (millions of pounds)
#From B. Richards 27.6 m radius Convert to million lbs
                                                  #27.6^2*pi*104.4653 gives psu area
BioYR <- BIOSV/1000000*2.20462/(25892*104.4653)   #194.89 was a scaler to account for proprotion of PSU sampled based on 20.2-m radius
s_eta2 <- (SE_SV/1000000*2.20462/(25892*104.4653))^2  #194.89 changed to 104.4653
s_CV <- sqrt(s_eta2)/BioYR
s_eta2log <- log(s_CV*s_CV+1)
s_lambda <- 1 #initial weighting on sd of survey estimate


#########################
# model parameters
#########################

Target_K_Prior_avg <- 29
CV_K <- 0.5

Target_r_Prior_avg <- 0.10
CV_r <- 0.25

Target_P1_Prior_avg <- 0.53
CV_P1 <- 0.2

M_shape <- 0.5
M_scale <- 0.5

process_shape <- 0.2
process_scale <- 0.1

observation_shape <- 0.2
observation_scale <- 1.0

q_lo <- 0.00001    
q_hi <- 100000

Target_rad_Prior_avg <-27.6
CV_rad <- 0.5

LB <- 0.6
UB <- 1.4

proj_LB <- 0.6  #should probably match what is in retro. part of model
proj_UB <- 1.4


pLIM_B <- 0.844

UC_ratio <-1.11  #recalculated this - average of 2014-2018

start_TAC <- 0.000
mesh_TAC <- 0.002
#mesh_TAC <- 0.005
NTAC <- 501
#NTAC <- 201


###############################################################
# Bundle Data
###############################################################

win.data <- list(
  
  NTIME = NTIME,
  Reported_Catch = Reported_Catch,
  UnrepCatch = UnrepCatch,
  
  CPUE_S1 = CPUE_S1,
  CPUE_S2 = CPUE_S2,
  CPUE_S1_REL_CV = CPUE_S1_REL_CV,
  CPUE_S2_REL_CV = CPUE_S2_REL_CV,
  NCPUE_S1_1 = NCPUE_S1_1,
  NCPUE_S1_MISS = NCPUE_S1_MISS,
  NCPUE_S1_2 = NCPUE_S1_2,
  NCPUE_S1_3 = NCPUE_S1_3,
  NSV_1 = NSV_1,
  NSV_2 = NSV_2,
  NSVa = NSVa,
  NSVb = NSVb,
  NSVc = NSVc,
  NSVd = NSVd,
  NSVe = NSVe,
  

  
  Target_K_Prior_avg = Target_K_Prior_avg,
  CV_K = CV_K,
  
  Target_r_Prior_avg = Target_r_Prior_avg, 
  CV_r = CV_r,
  
  Target_P1_Prior_avg = Target_P1_Prior_avg,
  CV_P1 = CV_P1,
  
  M_shape = M_shape,
  M_scale = M_scale,
  
  process_shape = process_shape,
  process_scale = process_scale,
  
  observation_shape = observation_shape,
  observation_scale = observation_scale,
  
  q_lo = q_lo,
  q_hi = q_hi,
  
  Target_rad_Prior_avg = Target_rad_Prior_avg,
  CV_rad = CV_rad,
  
  LB = LB,
  UB = UB,
  
  pLIM_B = pLIM_B,
  
  BioYR = BioYR,
  s_eta2log = s_eta2log,
  s_lambda = s_lambda,
  
  proj_LB = proj_LB,
  proj_UB = proj_UB,
  
  RC_2019 = RC_2019,
  
  RC_2020 = RC_2020,
  
  UC_ratio = UC_ratio,
  
  
  start_TAC = start_TAC,
  mesh_TAC = mesh_TAC,
  NTAC = NTAC
  
) # end data list

## END DATA
###################################################################  

# Analysis using WinBUGS - not used at the momment. Instead read separate BUGS file
# Define model written in WinBUGS code ------
model_code=paste0("model ",addname,".txt")
sink(model_code)  # sink diverts R output to a connection. 
cat("
    
    model
    {
    
    ##############################################################
    # PRIOR DISTRIBUTIONS
    ##############################################################
    
    # Lognormal prior for carrying capacity parameter, K
    #(P1)##########################################################
    K_Prior_Precision <- 1.0/log(1.0+CV_K*CV_K)
    K_Prior_avg <- log(Target_K_Prior_avg) - (0.5/K_Prior_Precision)
    K ~ dlnorm(K_Prior_avg,K_Prior_Precision)I(0.001,200.0)
    
    # Lognormal prior for intrinsic growth rate parameter, r
    #(P2)##########################################################
    r_Prior_Precision <- 1.0/log(1.0+CV_r*CV_r)
    r_Prior_avg <- log(Target_r_Prior_avg) - (0.5/r_Prior_Precision)
    r ~ dlnorm(r_Prior_avg,r_Prior_Precision)I(0.01,1.00)
    
    # Gamma prior for production shape parameter, M
    #(P3)##########################################################
    M ~ dgamma(M_shape, M_scale)
    
    # Uniform prior for CPUE catchability coefficients
    # in the interval (0.0001,10000), q1 and q2
    #(P4)##########################################################
    q1 ~ dunif(q_lo, q_hi)
    q2 ~ dunif(q_lo, q_hi)

    # Lognormal prior for the survey radius
    #(P4.b)##########################################################
    rad_Prior_Precision <- 1.0/log(1.0+CV_rad*CV_rad)
    rad_Prior_avg <- log(Target_rad_Prior_avg) - (0.5/rad_Prior_Precision)
    rad ~ dlnorm(rad_Prior_avg,rad_Prior_Precision)I(7.5,60.6)     #rather than 41.6, I think the upper bound would now be 60.6 based on publication
    q3 <- 250000/(rad*rad*3.14159)                                 #rlnorm for 27.6 gives (10.9 and 69.5 for 95% CI)
                                                                   #kept the 7.5 for now - not sure if justified in draft MS
    # Gamma prior for process error variance, sigma2
    #(P5)##########################################################
    isigma2 ~ dgamma(process_shape,process_scale)I(0.000001,1000000)
    sigma2 <- 1/isigma2
    
    # Gamma prior for observation error variance, tau2
    #(P6)##########################################################
    itau2_1   ~ dgamma(observation_shape,observation_scale)I(0.000001,1000000)
    tau2_1  <- 1/itau2_1
    
    itau2_2   ~ dgamma(observation_shape,observation_scale)I(0.000001,1000000)
    tau2_2   <- 1/itau2_2
    
    # Lognormal priors for unobserved states, the time series of proportions of K, P[]
    # MHI time catch series starts in FY1949 and ends in FY2013, n=65
    #(P7)##########################################################
    P1_Prior_Precision <- 1.0/log(1.0+CV_P1*CV_P1)
    P1_Prior_avg <-log(Target_P1_Prior_avg) - (0.5/P1_Prior_Precision)
    P[1] ~ dlnorm(P1_Prior_avg,P1_Prior_Precision) I(0.0001,10000)
    
    # Catch is uninformly distributed on the interval [lower, upper]
    #(P8)##########################################################
    lower[1] <- LB*UnrepCatch[1] + Reported_Catch[1]
    upper[1] <- UB*UnrepCatch[1] + Reported_Catch[1]
    Catch[1] ~ dunif(lower[1],upper[1])
    
    ##############################################################
    # PROCESS DYNAMICS
    ##############################################################
    for (i in 2:NTIME) {
    Pmean[i] <- log(max(P[i-1] + r*P[i-1]*(1-pow(P[i-1],M)) - Catch[i-1]/K,0.0001))
    P[i]  ~ dlnorm(Pmean[i],isigma2)I(0.0001,10000)
    lower[i] <- LB*UnrepCatch[i] + Reported_Catch[i]
    upper[i] <- UB*UnrepCatch[i] + Reported_Catch[i]
    Catch[i] ~ dunif(lower[i],upper[i])
    }

    Pmean2019 <- log(max(P[NTIME] + r*P[NTIME]*(1-pow(P[NTIME],M)) - Catch[NTIME]/K,0.0001))
    P2019 ~ dlnorm(Pmean2019,isigma2)I(0.0001,10000)
    C2019lo <- LB*0.15521783 + 0.177051  #last quantity is reported
    C2019hi <- UB*0.15521783 + 0.177051  #last quantity is reported
    Catch2019 ~ dunif(C2019lo,C2019hi)
    Pmean2020 <- log(max(P2019 + r*P2019*(1-pow(P2019,M)) - Catch2019/K,0.0001))
    P2020 ~ dlnorm(Pmean2020,isigma2)I(0.0001,10000)


    ##############################################################
    # LIKELIHOOD OF OBSERVED CPUE
    ##############################################################
    
    # deep7 bottomfish CPUE lLIKELIHOOD, NOT USED P[1:NCPUE_S1_1]
    #(L1)##########################################################
    # for (i in 1:NCPUE_S1_1) {
    # CPUE_mean[i] <- log(q1*K*P[i])
    # Precision_CPUE[i] <- itau2_1/(CPUE_S1_REL_CV[i]*CPUE_S1_REL_CV[i])
    # CPUE_S1[i] ~ dlnorm(CPUE_mean[i],Precision_CPUE[i])
    # LOG_RESID1[i] <- log(CPUE_S1[i]) - log(q1*K*P[i])
    # }
    
    # deep7 bottomfish CPUE lLIKELIHOOD, 1948-2003 P[(NCPUE_S1_1+NCPUE_S1_MISS+1):NCPUE_S1_2]
    #(L2)##########################################################
    for (i in (NCPUE_S1_1+NCPUE_S1_MISS+1):NCPUE_S1_2) {
    CPUE_mean[i] <- log(q1*K*P[i])
    Precision_CPUE[i] <- itau2_1/(CPUE_S1_REL_CV[i]*CPUE_S1_REL_CV[i])
    CPUE_S1[i] ~ dlnorm(CPUE_mean[i],Precision_CPUE[i])
    LOG_RESID1[i] <- log(CPUE_S1[i]) - log(q1*K*P[i])
    }
    
    # deep7 bottomfish CPUE lLIKELIHOOD, 2003-2019 P[(NCPUE_S1_2+1):NCPUE_S1_3]
    #(L3)##########################################################
    for (i in (NCPUE_S1_2):NCPUE_S1_3) {
    CPUE_mean2[i] <- log(q2*K*P[i])
    Precision_CPUE2[i] <- itau2_2/(CPUE_S2_REL_CV[i]*CPUE_S2_REL_CV[i])
    CPUE_S2[i] ~ dlnorm(CPUE_mean2[i],Precision_CPUE2[i])
    LOG_RESID2[i] <- log(CPUE_S2[i]) - log(q2*K*P[i])
    }

    # survey likelihood, for 2017-18 estimates  ##HOW IS THE PROJECTION HANDLED?? 5-11-20  P2019 and P2020 aren't in this
    #(L4)##########################################################
    for (i in (NSV_1):NSV_2){
    survey_mean[i] <- log(P[i]*K/(q3*25892))
    Precision_survey[i] <- (s_lambda*s_lambda)/s_eta2log[i-NSVa]
    BioYR[i-NSVa] ~ dlnorm(survey_mean[i],Precision_survey[i])
    LOG_RESID3[i] <- log(BioYR[i-NSVa]) - log(P[i]*K/(q3*25892))
    }

   #need to add P2020
  
    # survey likelihood, for 2019 and 2020 estimates
    #(L4)##########################################################
    survey_mean[NSVb] <- log(P2019*K/(q3*25892))
    Precision_survey19 <- (s_lambda*s_lambda)/s_eta2log[NSVc]
    BioYR[NSVc] ~ dlnorm(survey_mean[NSVb],Precision_survey19)
    LOG_RESID3[NSVb] <- log(BioYR[NSVc]) - log(P2019*K/(q3*25892))
    


    survey_mean[NSVd] <- log(P2020*K/(q3*25892))
    Precision_survey20 <- (s_lambda*s_lambda)/s_eta2log[NSVe]
    BioYR[NSVe] ~ dlnorm(survey_mean[NSVd] ,Precision_survey20)
    LOG_RESID3[NSVd] <- log(BioYR[NSVe]) - log(P2020*K/(q3*25892))
    

    
    # Compute LOG_RSS and LOG_RMSE
    ##############################################################
    # LOG_RSS1 <- inprod(LOG_RESID1[1:NCPUE_S1_1], LOG_RESID1[1:NCPUE_S1_1]) +
    LOG_RSS1 <- inprod(LOG_RESID1[(NCPUE_S1_1+NCPUE_S1_MISS+1):NCPUE_S1_2], LOG_RESID1[(NCPUE_S1_1+NCPUE_S1_MISS+1):NCPUE_S1_2]) 

    LOG_RSS2 <- inprod(LOG_RESID2[(NCPUE_S1_2):NCPUE_S1_3], LOG_RESID2[(NCPUE_S1_2):NCPUE_S1_3])
    
    LOG_RSS3 <- inprod(LOG_RESID3[NSV_1:NSVd], LOG_RESID3[NSV_1:NSVd])

    LOG_RMSE1 <- sqrt(LOG_RSS1/(NCPUE_S1_2-NCPUE_S1_MISS))
    
    LOG_RMSE2 <- sqrt(LOG_RSS2/(NCPUE_S1_3-(NCPUE_S1_2-1)))

    LOG_RMSE3 <- sqrt(LOG_RSS3)
    
    
    # Compute standardized log-scale residuals, predicted CPUE, and unscaled residuals
    ##############################################################
    # for (i in 1:NCPUE_S1_1) {
    # STD_LOG_RESID1[i] <- LOG_RESID1[i]/LOG_RMSE1
    # PRED_CPUE[i] <- exp(CPUE_mean[i]) ## PRED_CPUE[i] <- exp(log(CPUE_mean[i]))
    # RESID1[i] <- CPUE_S1[i] - PRED_CPUE[i]
    # }
    
    for (i in (NCPUE_S1_1+NCPUE_S1_MISS+1):NCPUE_S1_2) {
    STD_LOG_RESID1[i] <- LOG_RESID1[i]/LOG_RMSE1
    PRED_CPUE[i] <- exp(CPUE_mean[i]) ## PRED_CPUE[i] <- exp(log(CPUE_mean[i]))
    RESID1[i] <- CPUE_S1[i] - PRED_CPUE[i]
    }
    
    for (i in (NCPUE_S1_2):NCPUE_S1_3) {
    STD_LOG_RESID2[i] <- LOG_RESID2[i]/LOG_RMSE2
    PRED_CPUE2[i] <- exp(CPUE_mean2[i])
    RESID2[i] <- CPUE_S2[i] - PRED_CPUE2[i]
    }


    for (i in (NSV_1):NSVd){
    STD_LOG_RESID3[i] <- LOG_RESID3[i]/LOG_RMSE3
    PRED_Bio[i] <- exp(survey_mean[i])
    RESID3[i] <- BioYR[i-(NSV_1-1)] - PRED_Bio[i]
    }


    
    # Compute RSS and RMSE for MHI CPUE
    ##############################################################
    #RSS1 <- inprod(RESID1[1:NCPUE_S1_1], RESID1[1:NCPUE_S1_1]) +
    RSS1 <- inprod(RESID1[(NCPUE_S1_1+NCPUE_S1_MISS+1):NCPUE_S1_2], RESID1[(NCPUE_S1_1+NCPUE_S1_MISS+1):NCPUE_S1_2])
    
    RSS2 <- inprod(RESID2[(NCPUE_S1_2):NCPUE_S1_3], RESID2[(NCPUE_S1_2):NCPUE_S1_3])
    
    RSS3 <- inprod(RESID3[NSV_1:NSVd],RESID3[NSV_1:NSVd])

    RMSE1 <- sqrt(RSS1/(NCPUE_S1_2-NCPUE_S1_MISS))
    
    RMSE2 <- sqrt(RSS2/(NCPUE_S1_3-(NCPUE_S1_2-1)))

    RMSE3 <- sqrt(RSS3/(NSV_2+2-NSV_1))
    
    ###############################################################
    # STOCK ASSESSMENT QUANTITIES OF INTEREST
    ###############################################################
    
    # Compute exploitation rate and biomass time series
    #(QOI1)#########################################################
    # MHI 1948-2015 P[1:NTIME]
    for (i in 1:NTIME) {
    B[i] <- P[i]*K
    H[i] <- min(Catch[i]/B[i],0.999)
    F[i] <- -log(1-H[i])
    }  
    
    # Compute MSY reference points
    #(QOI2)#########################################################
    BMSY <- K*pow(M+1.0,(-1.0/M))
    MSY <- r*BMSY*(1.0-(1.0/(M+1.0)))
    HMSY <- min(r*(1.0-(1.0/(M+1.0))),0.999)
    PMSY <- BMSY/K
    FMSY <- -log(1-HMSY)
    CPUE_MSY <- q2*BMSY
    
    # Compute relative biomass and harvest, BSTATUS and HSTATUS
    #(QOI3)#########################################################
    for (i in 1:NTIME) {
    BSTATUS[i] <- B[i]/BMSY
    HSTATUS[i] <- H[i]/HMSY
    production[i] <- r*B[i]*(1-pow(P[i],M))
    }
    
    # Compute probabilities of H[i] > HMSY, B[i] < BMSY, 
    # and B[i] < pLIM_B*BMSY, a minimum biomass limit
    #(QOI4)##########################################################
    for (i in 1:NTIME) {
    pOFL_H[i] <- step(HSTATUS[i] - 1.0)
    pBMSY_B[i] <- step(1.0 - BSTATUS[i])
    pOFL_B[i] <- step(pLIM_B - BSTATUS[i])
    }


#################################################################
  #####   PROJECTIONS
    #################################################################
    # Fishing Year 2019 Projection
    proj_Pmean2019<- (max(P[NTIME] + r*P[NTIME]*(1-pow(P[NTIME],M)) - Catch[NTIME]/K,0.0001))
    #    proj_P2019  ~ dlnorm(proj_Pmean2019,isigma2)I(0.0001,10000)
    B[NTIME+1] <- proj_Pmean2019*K
    #    B[NTIME+1] <- proj_P2019*K
    
    
    #    RC[1] <- TAC_2019*prop_caught_2019
    UC[1] <- UC_ratio*RC_2019
    
    lower[NTIME+1] <- proj_LB*UC[1] + RC_2019
    upper[NTIME+1] <- proj_UB*UC[1] + RC_2019
    
    
    proj_C2019 ~ dunif(lower[NTIME+1],upper[NTIME+1])
    
    H[NTIME+1] <- min(proj_C2019/B[NTIME+1],0.999)
    
    BSTATUS[NTIME+1] <- B[NTIME+1]/BMSY 
    HSTATUS[NTIME+1] <- H[NTIME+1]/HMSY
    production[NTIME+1] <- r*B[NTIME+1]*(1-pow(proj_Pmean2019,M))
    
    pOFL_H[NTIME+1] <- step(HSTATUS[NTIME+1] - 1.0)
    pBMSY_B[NTIME+1] <- step(1.0 - BSTATUS[NTIME+1])
    pOFL_B[NTIME+1] <- step(pLIM_B - BSTATUS[NTIME+1])
    
    
    
    # Fishing Year 2020 Projection
    #################################################################
    proj_Pmean2020 <- (max(proj_Pmean2019 + r*proj_Pmean2019*(1-pow(proj_Pmean2019,M)) - proj_C2019/K,0.0001))
    #   proj_P2020  ~ dlnorm(proj_Pmean2020,isigma2)I(0.0001,10000)
    #   B[NTIME+2] <- proj_P2020*K
    B[NTIME+2] <- proj_Pmean2020*K
    
    #   RC[2] <- TAC_2020*prop_caught_2020
    UC[2] <- UC_ratio*RC_2020
    
    lower[NTIME+2] <- proj_LB*UC[2] + RC_2020
    upper[NTIME+2] <- proj_UB*UC[2] + RC_2020
    
    proj_C2020 ~ dunif(lower[NTIME+2],upper[NTIME+2])
    
    H[NTIME+2] <- min(proj_C2020/B[NTIME+2],0.999)
    
    BSTATUS[NTIME+2] <- B[NTIME+2]/BMSY 
    HSTATUS[NTIME+2] <- H[NTIME+2]/HMSY
    production[NTIME+2] <- r*B[NTIME+2]*(1-pow(proj_Pmean2020,M))
    
    pOFL_H[NTIME+2] <- step(HSTATUS[NTIME+2] - 1.0)
    pBMSY_B[NTIME+2] <- step(1.0 - BSTATUS[NTIME+2])
    pOFL_B[NTIME+2] <- step(pLIM_B - BSTATUS[NTIME+2])
    
    # Fishing Year 2021-2022 Projection #################################################################
    proj_lower <- proj_LB*UC_ratio
    proj_upper <- proj_UB*UC_ratio
    
    proj_Pmean <- (max(proj_Pmean2020 + r*proj_Pmean2020*(1-pow(proj_Pmean2020,M)) - proj_C2020/K,0.0001))
    #  proj_P ~ dlnorm(proj_Pmean,isigma2)I(0.0001,10000)
    B[NTIME+3] <- proj_Pmean*K  #2021 biomass
    
    BSTATUS[NTIME+3] <- B[NTIME+3]/BMSY 
    production[NTIME+3] <- r*B[NTIME+3]*(1-pow(proj_Pmean,M)) 
    pBMSY_B[NTIME+3] <- step(1.0 - BSTATUS[NTIME+3])
    pOFL_B[NTIME+3] <- step(pLIM_B - BSTATUS[NTIME+3])
    
    for (j in 1:NTAC)
    {
    #2021-2022
    proj_TAC[j] <- start_TAC+mesh_TAC*(j-1)
    
    proj_UC_ratio1[j] ~ dunif(proj_lower,proj_upper)
    proj_UC1[j] <- proj_UC_ratio1[j]*proj_TAC[j]
    proj_C1[j] <- proj_TAC[j] + proj_UC1[j]
    proj_H1[j] <- min(proj_C1[j]/B[NTIME+3],0.999)
    proj_HSTATUS1[j] <- proj_H1[j]/HMSY
    proj_pOFL_H1[j] <- step(proj_HSTATUS1[j] - 1.0)  
    
    
    #    proj_Pmean_2022[j] <- log(max(proj_P + r*proj_P*(1-pow(proj_P,M)) - proj_C1[j]/K,0.0001))  #JS added process error here and henceforth
    #    proj_P_2022[j] ~ dlnorm(proj_Pmean_2022[j],isigma2)I(0.0001,10000)
    #    proj_B2022[j] <- K*proj_P_2022[j]
    
    proj_P2022[j] <- max(proj_Pmean + r*proj_Pmean*(1-pow(proj_Pmean,M)) - proj_C1[j]/K,0.0001)  
    proj_B2022[j] <- proj_P2022[j]*K
    
    #    proj_BSTATUS[j] <- proj_B[j]/BMSY
    
    proj_BSTATUS[j] <- proj_B2022[j]/BMSY   ##############DO WE WANT THE STATUS UNTIL 2022??#######################
    proj_pOFL_B[j] <- step(pLIM_B - proj_BSTATUS[j])
    proj_UC_ratio2[j] ~ dunif(proj_lower,proj_upper)
    proj_UC2[j] <- proj_UC_ratio2[j]*proj_TAC[j]
    proj_C2[j] <- proj_TAC[j] + proj_UC2[j]
    proj_H2[j] <- min(proj_C2[j]/proj_B2022[j],0.999)
    #    proj_H2[j] <- min(proj_C2[j]/proj_B[j],0.999)
    proj_HSTATUS2[j] <- proj_H2[j]/HMSY
    proj_pOFL_H2[j] <- step(proj_HSTATUS2[j] - 1.0)
    
    
    
    #2023-2025   
    
    
    #    proj_Pmean_2023[j] <- log(max(proj_P_2022[j] + r*proj_P_2022[j]*(1-pow(proj_P_2022[j],M)) - proj_C2[j]/K,0.0001))
    #   proj_P_2023[j] ~ dlnorm(proj_Pmean_2023[j],isigma2)I(0.0001,10000)
    #    proj_B2023[j] <- K*proj_P_2023[j]
    
    
    
    proj_P2023[j] <- max(proj_P2022[j] + r*proj_P2022[j]*(1-pow(proj_P2022[j],M)) - proj_C2[j]/K,0.0001) 
    proj_B2023[j]<- proj_P2023[j]*K
    
    proj_UC_ratio3[j] ~ dunif(proj_lower,proj_upper)
    proj_UC3[j] <- proj_UC_ratio3[j]*proj_TAC[j]
    proj_C3[j] <- proj_TAC[j] + proj_UC3[j]
    proj_H3[j] <- min(proj_C3[j]/proj_B2023[j],0.999)
    proj_HSTATUS3[j] <- proj_H3[j]/HMSY
    proj_pOFL_H3[j] <- step(proj_HSTATUS3[j] - 1.0)
    proj_BSTATUS3[j]<- proj_B2023[j]/BMSY
    proj_pOFL_B3[j] <- step(pLIM_B - proj_BSTATUS3[j])
    
    #    proj_Pmean_2021[j] <- log(max(proj_P_2020[j] + r*proj_P_2020[j]*(1-pow(proj_P_2020[j],M)) - proj_C3[j]/K,0.0001))
    #    proj_P_2021[j] ~ dlnorm(proj_Pmean_2021[j],isigma2)I(0.0001,10000)
    #    proj_B2021[j] <- K*proj_P_2021[j]
    
    proj_P2024[j] <- max(proj_P2023[j] + r*proj_P2023[j]*(1-pow(proj_P2023[j],M)) - proj_C3[j]/K,0.0001) 
    proj_B2024[j] <-proj_P2024[j]*K
    
    proj_UC_ratio4[j] ~ dunif(proj_lower,proj_upper)
    proj_UC4[j] <- proj_UC_ratio4[j]*proj_TAC[j]
    proj_C4[j] <- proj_TAC[j] + proj_UC4[j]
    proj_H4[j] <- min(proj_C4[j]/proj_B2024[j],0.999)
    proj_HSTATUS4[j] <- proj_H4[j]/HMSY
    proj_pOFL_H4[j] <- step(proj_HSTATUS4[j] - 1.0)
    proj_BSTATUS4[j]<- proj_B2024[j]/BMSY
    proj_pOFL_B4[j] <- step(pLIM_B - proj_BSTATUS4[j])
    
    proj_P2025[j] <- max(proj_P2024[j] + r*proj_P2024[j]*(1-pow(proj_P2024[j],M)) - proj_C4[j]/K,0.0001)
    proj_B2025[j] <- proj_P2025[j]*K
    
    #    proj_Pmean_2022[j] <- log(max(proj_P_2021[j] + r*proj_P_2021[j]*(1-pow(proj_P_2021[j],M)) - proj_C4[j]/K,0.0001))
    #    proj_P_2022[j] ~ dlnorm(proj_Pmean_2022[j],isigma2)I(0.0001,10000)
    #    proj_B2022[j] <- K*proj_P_2022[j]
    
    
    proj_UC_ratio5[j] ~ dunif(proj_lower,proj_upper)
    proj_UC5[j] <- proj_UC_ratio5[j]*proj_TAC[j]
    proj_C5[j] <- proj_TAC[j] + proj_UC5[j]
    proj_H5[j] <- min(proj_C5[j]/proj_B2025[j], 0.999)
    proj_HSTATUS5[j] <- proj_H5[j]/HMSY
    proj_pOFL_H5[j] <- step(proj_HSTATUS5[j] - 1.0)
    
    proj_BSTATUS5[j] <- proj_B2025[j]/BMSY   ##############DO WE WANT THE STATUS UNTIL 2022??#######################
    proj_pOFL_B5[j] <- step(pLIM_B - proj_BSTATUS5[j])
    
    }

    
    ########################
    } ## END OF WinBUGS MODEL
    
    ",fill=TRUE)
sink()      # ends the last diversion



###################################################################
# END OF CODE/MODEL
###################################################################


#######################################################################  ------
############ Create list of inits for WinBUGS use #####################                                
#######################################################################

inits <- list(    # create inits list of functions
  
  ## Initial Condition 1 
  
  list(
    
    
    Catch=c(0.890568408,0.780141322,0.827136735,0.76223799,0.661001365,
            0.676712544,0.55389128,0.710439819,0.856515787,0.563433499,
            0.517070196,0.42887519,0.339205311,0.4469463,0.560498174,
            0.542606506,0.585471579,0.468168603,0.660834647,0.524553112,
            0.493078085,0.429513477,0.365754118,0.652472663,0.49518941,
            0.68527923,0.639130203,0.685372642,0.718483749,0.872808368,
            0.788074103,0.739263202,0.948284067,0.93863966,1.192176718,
            0.933642904,1.249188416,1.244535433,1.534016004,1.60281839,
            1.604606184,1.220957021,0.840975704,0.971630742,0.716315491,
            0.893114185,0.976257722,0.765849379,0.850127531,0.781433882,
            0.528538985,0.750220277,0.569860369,0.43489043,0.514820142,
            0.396826878,0.470045616,0.355597079,0.421877928,0.40172532,
            0.538372806,0.438659038,0.586252204,0.445221546,0.445352322,
            0.652537384,0.640379572,0.557556757,0.523622913,0.493875921),
            
    
    Catch2019 = 0.332268904,
    
    r=0.05,
    
    P=c(rep(0.5,32), rep(0.5, NTIME-32)), 
    
    P2019=0.5,
    P2020=0.5,
    
    K=45.0,
    
    M=1.0,
    
    q1=10.0,
    q2=10.0,
    rad=27.6,
    
    isigma2=100,
    
    itau2_1=100,
    itau2_2=100,
    
    proj_C2019= 0.332268904,   #0.3493428 if using 1.11 UC ratio
    
    proj_C2020= 0.449923, 
    
    proj_UC_ratio1=rep(1.11, 501),
    
    proj_UC_ratio2=rep(1.11, 501),
    
    proj_UC_ratio3=rep(1.11, 501),
    
    proj_UC_ratio4=rep(1.11, 501),
    
    proj_UC_ratio5=rep(1.11, 501)
    
    
    
    
  )##END init 1
  
  ## Initial Condition 2 
  
  ,list(
    
    Catch=c(0.890568408,0.780141322,0.827136735,0.76223799,0.661001365,
            0.676712544,0.55389128,0.710439819,0.856515787,0.563433499,
            0.517070196,0.42887519,0.339205311,0.4469463,0.560498174,
            0.542606506,0.585471579,0.468168603,0.660834647,0.524553112,
            0.493078085,0.429513477,0.365754118,0.652472663,0.49518941,
            0.68527923,0.639130203,0.685372642,0.718483749,0.872808368,
            0.788074103,0.739263202,0.948284067,0.93863966,1.192176718,
            0.933642904,1.249188416,1.244535433,1.534016004,1.60281839,
            1.604606184,1.220957021,0.840975704,0.971630742,0.716315491,
            0.893114185,0.976257722,0.765849379,0.850127531,0.781433882,
            0.528538985,0.750220277,0.569860369,0.43489043,0.514820142,
            0.396826878,0.470045616,0.355597079,0.421877928,0.40172532,
            0.538372806,0.438659038,0.586252204,0.445221546,0.445352322,
            0.652537384,0.640379572,0.557556757,0.523622913,0.493875921),
    
    Catch2019 = 0.332268904,
    
    r=0.15,
    
    P=c(rep(0.5,32), rep(0.5, NTIME-32)), 
    
    P2019=0.5,
    P2020=0.5,
    
    K=15.0,
    
    M=1.0,
    
    q1=10.0,
    q2=10.0,
    rad=27.6,
    
    isigma2=100,
    
    itau2_1=100,
    itau2_2=100,
    
    proj_C2019= 0.332268904,   #0.3493428 if using 1.11 UC ratio
    
    proj_C2020= 0.449923, 
    
    proj_UC_ratio1=rep(1.11, 501),
    
    proj_UC_ratio2=rep(1.11, 501),
    
    proj_UC_ratio3=rep(1.11, 501),
    
    proj_UC_ratio4=rep(1.11, 501),
    
    proj_UC_ratio5=rep(1.11, 501)
    
    
  )##END init 2
  
  ## Initial Condition 3 
  
  ,list(
    
    Catch=c(0.890568408,0.780141322,0.827136735,0.76223799,0.661001365,
            0.676712544,0.55389128,0.710439819,0.856515787,0.563433499,
            0.517070196,0.42887519,0.339205311,0.4469463,0.560498174,
            0.542606506,0.585471579,0.468168603,0.660834647,0.524553112,
            0.493078085,0.429513477,0.365754118,0.652472663,0.49518941,
            0.68527923,0.639130203,0.685372642,0.718483749,0.872808368,
            0.788074103,0.739263202,0.948284067,0.93863966,1.192176718,
            0.933642904,1.249188416,1.244535433,1.534016004,1.60281839,
            1.604606184,1.220957021,0.840975704,0.971630742,0.716315491,
            0.893114185,0.976257722,0.765849379,0.850127531,0.781433882,
            0.528538985,0.750220277,0.569860369,0.43489043,0.514820142,
            0.396826878,0.470045616,0.355597079,0.421877928,0.40172532,
            0.538372806,0.438659038,0.586252204,0.445221546,0.445352322,
            0.652537384,0.640379572,0.557556757,0.523622913,0.493875921),
    
    Catch2019 = 0.332268904,
    
    r=0.10,
    
    P=c(rep(0.5,32), rep(0.5, NTIME-32)), 
    
    P2019=0.5,
    P2020=0.5,
    
    K=30.0,
    
    M=1.0,
    
    q1=10.0,
    q2=10.0,
    rad=27.6,
    
    isigma2=100,
    
    itau2_1=100,
    itau2_2=100,
    
    proj_C2019= 0.332268904,   #0.3493428 if using 1.11 UC ratio
    
    proj_C2020= 0.449923, 
    
    proj_UC_ratio1=rep(1.11, 501),
    
    proj_UC_ratio2=rep(1.11, 501),
    
    proj_UC_ratio3=rep(1.11, 501),
    
    proj_UC_ratio4=rep(1.11, 501),
    
    proj_UC_ratio5=rep(1.11, 501)
    
  )##END init 3
)  ## close list of functions

##### end initials function ############################################
########################################################################


## Parameters to estimate
########################################################################

params <- c(
  
  ## model parameters ##
#  "K","r","M", "q1","q2","sigma2","tau2_1","tau2_2","q3","rad",
  
  ## time-series derived variables ##
#  "P","B","H","PRED_CPUE","PRED_CPUE2","PRED_Bio",
  
  ## management metrics ##
#  "MSY","PMSY","BMSY","HMSY","BSTATUS","HSTATUS","FMSY", 
#  "pOFL_H","pOFL_B","pBMSY_B",                                                                 
  
  ## statistics and diagnoses ##
#  "STD_LOG_RESID1", "STD_LOG_RESID2", "STD_LOG_RESID3",                             #"STD_LOG_RESID3",
#  "LOG_RESID1", "LOG_RESID2", "LOG_RESID3","RESID1", "RESID2", "RESID3",              #"LOG_RESID3", "RESID3",
 # "LOG_RSS1", "LOG_RSS2", "LOG_RSS3", "LOG_RMSE1", "LOG_RMSE2", "LOG_RMSE3",             #"LOG_RSS3",
#   "RSS1", "RSS2","RSS3",  "RMSE1", "RMSE2" ,"RMSE3"                          # "LOG_RMSE3",  "RSS3",  ,"RMSE3"
  

  #Need to divide out due to number of nodes being monitored, do
  #B and H for one, and rerun with poflB and poflH as the other
  "B", #always monitor B as a check
  
  "proj_B2022",
  
  "proj_B2023",
  
  "proj_B2024",
  
  "proj_B2025",
  
  "proj_H1",
  
  "proj_H2",
  
  "proj_H3",
  
  "proj_H4",
  
  "proj_H5"
  
#  "proj_pOFL_H1",
  
#  "proj_pOFL_H2",
  
#  "proj_pOFL_H3",
  
#  "proj_pOFL_H4",
  
#  "proj_pOFL_H5",
  
#  "pOFL_B",
  
#  "proj_pOFL_B",
  
#  "proj_pOFL_B3",
  
#  "proj_pOFL_B4",
  
#  "proj_pOFL_B5"
  
  
  #   "proj_C1", "proj_H1", "proj_HSTATUS1", "proj_P", "proj_P2017",
  #   "proj_BSTATUS", "proj_C2", "proj_pOFL_B", 
  #   "proj_H2", "proj_HSTATUS2", "proj_pOFL_H2"#,
  
)


begin_time = proc.time()[3]  
nc <- length(inits)    # Number of Markov chains, default is 3 
################################################################### 
# Start Gibbs sampling, cycle through the initials

bugs(win.data,inits,params,model_code,n.chains=nc,n.iter=ni,n.burnin=nb,n.thin=nt,
     debug=TRUE,codaPkg=FALSE,bugs.directory="c:/WinBUGS/",   #c:/Program Files/WinBUGS14/  #c:/WinBUGS/
     working.directory=src.dir)

###################################################################

end_time = proc.time()[3]
print(paste("RUN_COST = ",(end_time-begin_time)/60," mins",sep=""))

#######################################################################
## Post-model processing!

## Rename (add model name) and copy files from to SAP shared network drive;
##   Make sure you're connected to the network drives

# file.names <- dir(src.dir)
# rem=which(file.names%in%c("codaIndex.txt",paste0("model ",addname,".txt")))
# file.names=file.names[-rem]
# sapply(file.names, function(x) {
#   file.rename(from=paste(src.dir, x, sep=''),
#               to=paste(src.dir,addname,"_",x, sep=''))}) # add model name

     