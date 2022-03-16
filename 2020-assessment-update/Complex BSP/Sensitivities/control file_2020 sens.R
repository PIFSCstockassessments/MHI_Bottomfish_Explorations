##############################################################
# Control file used to automate running scenarios
# for the sensitivity version of the base BSP model "d7_2017_base"
# which is "d7_2017_sensitivity"

# Brian Langseth, PIFSC, June 2017

# Updated by John Syslo, PIFSC, June 2020
#############################################################


#####For sensitivities where values change#####
# Unreported catch
# Unreported catch error
# Natural mortality
#
# Priors:
#   Carrying capacity (K)
#   Intrinsic rate of growth (r)
#   Production shape (M)
#   Process error (tau2)
#   Observation error (sigma2)
#   Initial proportion of carrying capacity (P1)
###############################################

rm(list=ls())



#Unreported catch
####  RUN as of 6-17-20 7-18??
#file that runs the sensitivities
source("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sensitivity_unrepCatch.R")
#need to specify 2019 unreported catch separately from the data file for each scenario:
#unrep_2019<-c(0.442628,0.040786604,0.195672,0)

basedata = read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\data_NA.csv",header=T)
surveydata=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Survey_data.csv",header=T)


#Scenario 1
#data1=read.csv("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\model code\\SPM code\\sensitivities\\d7_2017_data_nrepC1.csv",header=T)
data1=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\data_nrep_C1.csv",header=T)
init_c1=data1[-1,"Catch"]+data1[-1,"UnrepCatch"]
d7_2020_sens(data=data1,survey_data=surveydata,filename='d7_2020_sens_unrepC1',
             K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=0.53,Mscale=0.5,
             tau2_scale=1,sigma2_scale=0.1,nat_mort=0.156,
             unrep_c_err=0.4,dir_err=0,sweight=1,radius_Prior_avg=27.6,
             radius_CV=0.5,initC=init_c1,unrep2019=0.442628)


#Scenario 3
#data3=read.csv("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\model code\\SPM code\\sensitivities\\d7_2017_data_nrepC3.csv",header=T)
data3=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\data_nrep_C3.csv",header=T)
init_c3=data3[-1,"Catch"]+data3[-1,"UnrepCatch"]
d7_2020_sens(data=data3,survey_data=surveydata,filename='d7_2020_sens_unrepC3',
             K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=0.53,Mscale=0.5,
             tau2_scale=1,sigma2_scale=0.1,nat_mort=0.156,
             unrep_c_err=0.4,dir_err=0,sweight=1,radius_Prior_avg=27.6,
             radius_CV=0.5,initC=init_c3,unrep2019=0.040786604)
#Scenario 4
#data4=read.csv("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\model code\\SPM code\\sensitivities\\d7_2017_data_nrepC4.csv",header=T)
data4=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\data_nrep_C4.csv",header=T)
init_c4=data4[-1,"Catch"]+data4[-1,"UnrepCatch"]
d7_2020_sens(data=data4,survey_data=surveydata,filename='d7_2020_sens_unrepC4',
             K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=0.53,Mscale=0.5,
             tau2_scale=1,sigma2_scale=0.1,nat_mort=0.156,
             unrep_c_err=0.4,dir_err=0,sweight=1,radius_Prior_avg=27.6,
             radius_CV=0.5,initC=init_c4,unrep2019=0.195672)
#Scenario 5
#data5=read.csv("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\model code\\SPM code\\sensitivities\\d7_2017_data_nrepC5.csv",header=T)
data5=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\data_nrep_C5.csv",header=T)
init_c5=data5[-1,"Catch"]+data5[-1,"UnrepCatch"]
d7_2020_sens(data=data5,survey_data=surveydata,filename='d7_2020_sens_unrepC5',
             K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=0.53,Mscale=0.5,
             tau2_scale=1,sigma2_scale=0.1,nat_mort=0.156,
             unrep_c_err=0.4,dir_err=0,sweight=1,radius_Prior_avg=27.6,
             radius_CV=0.5,initC=init_c5,unrep2019=0)


#source new model code - will remove data [rm()]

source("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Sensitivity\\d7_2020_sensitivity.R")

basedata = read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\data_NA.csv",header=T)
surveydata=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Survey_data.csv",header=T)



###
#Unreported catch error - uce scenarios run on 6-16-20
uce=c(0.01,0.2,0.6)
for(i in 1:length(uce)) {
  d7_2020_sens(data=basedata,survey_data=surveydata,filename=paste0('sens_unrepCe',uce[i]),
             K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=0.53,Mscale=0.5,
             tau2_scale=1,sigma2_scale=0.1,nat_mort=0.156,
             unrep_c_err=uce[i],dir_err=0,sweight=1,radius_Prior_avg=27.6,
             radius_CV=0.5)
}
#Directional unreported catch errors - 2020 - run 6-16-2020 - needed to shorten file name to avoid errors
d7_2020_sens(data=basedata,survey_data=surveydata,filename='unrepCe_neg25',
             K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=0.53,Mscale=0.5,
             tau2_scale=1,sigma2_scale=0.1,nat_mort=0.156,
             unrep_c_err=0.4,dir_err=-0.25,sweight=1,radius_Prior_avg=27.6,
             radius_CV=0.5)
d7_2020_sens(data=basedata,survey_data=surveydata,filename='unrepCe_pos25',
             K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=0.53,Mscale=0.5,
             tau2_scale=1,sigma2_scale=0.1,nat_mort=0.156,
             unrep_c_err=0.4,dir_err=0.25,sweight=1,radius_Prior_avg=27.6,
             radius_CV=0.5)

###
#Natural mortality
d7_2020_sens(data=basedata,survey_data=surveydata,filename='2020_sens_nat25',
             K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=0.53,Mscale=0.5,
             tau2_scale=1,sigma2_scale=0.1,nat_mort=0.25,
             unrep_c_err=0.4,dir_err=0,sweight=1,radius_Prior_avg=27.6,
             radius_CV=0.5)



###  RAN these on 6-12-20 and into the weekend
#PRIORS
#k=c(14.5,20,25,30,35,40,43.5)

k=c(29-29*0.5,29-29*0.25,29+29*0.25,29+29*0.5) #adjust by +/- 25, 50%

r=c(0.05,0.15,0.25) #adjust by +/- 50, +150%

M=c(0.33,0.4,0.67,1.0)      #0.5/c(1.5,1.25,0.75,0.5)

tau2=c(0.01,0.1,10,100)

sigma2=c(0.001,0.01,1,10) # could do 5 instead of 10..."undefined real"

P1=c(0.265, 0.4, 0.66, 0.795) 



for(i in 1:length(k)) { #Carrying capacity
  d7_2020_sens(data=basedata,survey_data=surveydata,filename=paste0('origq_sens_k',k[i]),
               K_Prior_avg=k[i],r_Prior_avg=0.1,P1_Prior_avg=0.53,Mscale=0.5,
               tau2_scale=1,sigma2_scale=0.1,nat_mort=0.156,
               unrep_c_err=0.4,dir_err=0,sweight=1,radius_Prior_avg=27.6,
               radius_CV=0.5) 
}
for(i in 1:length(r)) { #Intrinsic rate of growth
  d7_2020_sens(data=basedata,survey_data=surveydata,filename=paste0('d7_2020_sens_r',r[i]),
               K_Prior_avg=29.0,r_Prior_avg=r[i],P1_Prior_avg=0.53,Mscale=0.5,
               tau2_scale=1,sigma2_scale=0.1,nat_mort=0.156,
               unrep_c_err=0.4,dir_err=0,sweight=1,radius_Prior_avg=27.6,
               radius_CV=0.5) 
}
for(i in 1:length(M)) { #Production shape parameter
  d7_2020_sens(data=basedata,survey_data=surveydata,filename=paste0('d7_2020_sens_M',M[i]),
               K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=0.53,Mscale=M[i],
               tau2_scale=1,sigma2_scale=0.1,nat_mort=0.156,
               unrep_c_err=0.4,dir_err=0,sweight=1,radius_Prior_avg=27.6,
               radius_CV=0.5) 
}
for(i in 1:length(tau2)) { #Observation error
  d7_2020_sens(data=basedata,survey_data=surveydata,filename=paste0('d7_2020_sens_tau',tau2[i]),
               K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=0.53,Mscale=0.5,
               tau2_scale=tau2[i],sigma2_scale=0.1,nat_mort=0.156,
               unrep_c_err=0.4,dir_err=0,sweight=1,radius_Prior_avg=27.6,
               radius_CV=0.5) 
}
for(i in 1:length(sigma2)) { #Process error 
  d7_2020_sens(data=basedata,survey_data=surveydata,filename=paste0('d7_2020_sens_sig',sigma2[i]),
               K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=0.53,Mscale=0.5,
               tau2_scale=1,sigma2_scale=sigma2[i],nat_mort=0.156,
               unrep_c_err=0.4,dir_err=0,sweight=1,radius_Prior_avg=27.6,
               radius_CV=0.5) 
}

for(i in 1:length(P1)) { #P1 mean 
  d7_2020_sens(data=basedata,survey_data=surveydata,filename=paste0('d7_2020_sens_P1_',P1[i]),
               K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=P1[i],Mscale=0.5,
               tau2_scale=1,sigma2_scale=0.1,nat_mort=0.156,
               unrep_c_err=0.4,dir_err=0,sweight=1,radius_Prior_avg=27.6,
               radius_CV=0.5) 
}      
  
  #Belief in estimate from survey
  d7_2020_sens(data=basedata,survey_data=surveydata,filename='sens_survey',
               K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=0.53,Mscale=0.5,
               tau2_scale=1,sigma2_scale=0.1,nat_mort=0.156,
               unrep_c_err=0.4,dir_err=0,sweight=1,radius_Prior_avg=27.6,
               radius_CV=0.01) 
 

  
  
  
      
#####For sensitivities where structure changes#####
# Random walk q
# Uniform(0,100) on process and observation error STDEV
###################################################

  
  NEED TO UPDATE THE SOURCES:      
###
#Random walk q - base conditions with the exception of treatment of q
source("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\model code\\SPM code\\d7_2017_varq.R")

####
#Uniform prior for error STDEVs - base conditions with exception of distributions of tau2 and sigma2
source("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\model code\\SPM code\\d7_2017_uniform_error.R")


