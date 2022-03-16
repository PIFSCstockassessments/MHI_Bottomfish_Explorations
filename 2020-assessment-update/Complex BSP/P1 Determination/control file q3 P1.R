##############################################################
# Control file used to automate running scenarios
# for the sensitivity version of the base BSP model "d7_2017_base"
# which is "d7_2017_48_q3_sensitivity"

# Brian Langseth, PIFSC, June 2017

# Updated November 15, 2017 to re-estimate p1 for the model with q3
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
library(R2WinBUGS)
#source("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\model code\\SPM code\\d7_2017_48sensitivity.R")

source("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\P1 determination\\d7_2020_48_FINAL_REMLF_q3_sensitivity.R")
##For determining initial proportion of carrying capacity
basedata = read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\data_NA.csv",header=T)
surveydata=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\Complex BSP\\Survey_data.csv",header=T)

###
#PRIOR on Initial Proportion of Carrying Capacity
P1=seq(0.1,1.0,0.1)

begin_time = proc.time()[3]  

for(i in 1:length(P1)) { #Initial proportion of carrying capacity
  d7_2017_sens(data=basedata,survey_data=surveydata,filename=paste0('proper_rad_NA_P1',P1[i]),
               K_Prior_avg=29.0,r_Prior_avg=0.1,P1_Prior_avg=P1[i],Mscale=0.5,
               tau2_scale=1,sigma2_scale=0.1,nat_mort=0.156,
               unrep_c_err=0.4,dir_err=0,sweight=1) 
}


end_time = proc.time()[3]
print(paste("RUN_COST = ",(end_time-begin_time)/60," mins",sep=""))




