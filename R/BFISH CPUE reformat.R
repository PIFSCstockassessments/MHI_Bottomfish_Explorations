#
rm(list=ls())

##RESEARCH FISHING ################################################################ 
#Use drift as base file - 2 drifts per sample ID, but can comibne later if needed
#go through catch and sample files to combine needed variables

BFISH_D<-read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\CRF_DRIFT.csv") #drift-specific information
BFISH_S<-read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\CRF_SAMPLE.csv") #sample-specific (i.e., PSU) information
BFISH_C<-read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\CRF_CATCH.csv") #catches (one row per individual)

#"Island and other PSU variables can be found here:  
PSU_table <- read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\BFISH PSU lookup table.csv")


#start with sample ID file.  This is at the PSU level
#Grab catch from catch file, grab PSU information, including lat and long, from the PSU file
#DRift file has start end lat/long/depth for each drift, depth is already aggregated as mean of means by sample ID in BFISH_S
     #lat long more easily gotten from the PSU file

valid<-function(x) {length(x[!is.na(x)])}
temp1<-tapply(BFISH_C$BFISH,list(BFISH_C$SAMPLE_ID,BFISH_C$SPECIES_CD),valid)  #1543 fish without spp code (i.e., "NONE").  Leave in dataframe for now 
temp2<-cbind(temp1,row.names(temp1))
colnames(temp2)[71]<-"SAMPLE_ID"
rownames(temp2) <- NULL

temp3<-merge(BFISH_S,temp2,by="SAMPLE_ID") 

#write.csv(temp3,"C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\Res_Fish_PSU.csv")


#presence not an issue as it is in camera (i.e., catch is presence.  can probably remove comment field, as well as "none" for spp code
#match up with appropriate PSU data

#RF<-read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\Res_Fish_PSU.csv")

#merge to bring in all PSU data
#RF_PSU<-merge(RF,PSU_table,by="PSU",all.x=TRUE)
#write.csv(RF_PSU,"C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\Res_Fish_PSU.csv")


RF_PSU<-read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\Res_Fish_PSU.csv")

#########################################################################################################################
#now match up with PacIOOS data from Netcdf files
#what are the earliest and latest BFISH dates by year for NEtCDF download?

tapply(RF_PSU$SAMPLE_DATE,RF_PSU$BFISH,range)
tapply(RF_PSU$lat_deg,RF_PSU$BFISH,range,na.rm=T)
tapply(RF_PSU$lon_deg,RF_PSU$BFISH,range,na.rm=T)

#$BFISH_2016_F
#[1] "2016-09-04" "2016-11-12"

#$BFISH_2016_S
#[1] "2016-02-06" "2016-03-23"

#$BFISH_2017_F
#[1] "2017-08-27" "2017-10-25"

#$BFISH_2018_F
#[1] "2018-08-15" "2018-11-28"

#$BFISH_2019_F
#[1] "2019-09-14" "2019-11-15"

#$BFISH_2020_F
#[1] "2020-08-15" "2020-11-15"



#####CAMERA SAMPLING REFORMATTING##############################################################################################################################
#Drop time is concatenated two digit hour, two digit minute, two digit second, in UTC time. will need to reformat
#there are 2 drops to sample a PSU
#will need to combine vaiables for the drops within a PSU 
#drop code is key to link the sample file with the abundance file (Max_N ), so need to bring in abundance and then sum by PSU

#this script keeps one record for PSUs sampled on 2 different days (dplyr takes the first)
#It may be best to align with covariate data at Drop_code level to get values at each day, then run this script to combine within PSU and reduce to unique psus.  Tnhis would  


rm(list=ls())
#BFISH_CAM_S<-read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\CAM_SAMPLE.csv")

#Version below has additional time field "DROP_TIME_NEW" to account for inconsistencies
BFISH_CAM_S<-read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\CAM_SAMPLE_TIME.csv")

BFISH_CAM_N<-read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\CAM_MAXN.csv")

head(BFISH_CAM_S)
dim((BFISH_CAM_S))

length(unique(BFISH_CAM_S$DROP_CD))
length(unique(BFISH_CAM_N$DROP_CD))


head(BFISH_CAM_N)
dim((BFISH_CAM_N))

#strive for 2 drops per PSU, so PSU is the experimental unit

#Deep 7 species codes
#sp<-c("PRFI","PRSI","PRZO","ETCA","ETCO","HYQU","APRU") #P. filamentosus, P. sieboldii, P. zonatus,  E. carbunculus, E. coruscans,  H. querdus, A. rutilans


temp1<-tapply(BFISH_CAM_N$MAXN,list(BFISH_CAM_N$DROP_CD,BFISH_CAM_N$SPECIES_CD),sum, na.rm=T) #this works for the abundance variable

#convert NA to 0 and apply length function
BFISH_CAM_N_temp<-BFISH_CAM_N
BFISH_CAM_N_temp$Ntemp<-1
#1's for presence, NA otherwise 
temp_presence<-tapply(BFISH_CAM_N_temp$Ntemp,list(BFISH_CAM_N_temp$DROP_CD,BFISH_CAM_N_temp$SPECIES_CD),sum, na.rm=T) #this works for the abundance variable
colnames(temp_presence)<-paste0(colnames(temp_presence),"_pres")
#convert the NAs to 0s
temp_presence[is.na(temp_presence)] <- 0


#Deep 7 species codes to potentially sum later
#sp<-c("PRFI","PRSI","PRZO","ETCA","ETCO","HYQU","APRU") #P. filamentosus, P. sieboldii, P. zonatus,  E. carbunculus, E. coruscans,  H. querdus, A. rutilans

#merge the catch with the camera sample file
temp1_drop<-cbind(temp1,row.names(temp1))
colnames(temp1_drop)[124]<-"DROP_CD"
temp2<-merge(BFISH_CAM_S,temp1_drop,by="DROP_CD")


#merge the presence/absence with the file
temp_presence_drop<-cbind(temp_presence,row.names(temp_presence))
colnames(temp_presence_drop)[124]<-"DROP_CD"
temp3<-merge(temp2,temp_presence_drop,by="DROP_CD")

#now combine drops with a PSU -
#create unique identifier
temp3$Cam_ID<-paste0(temp3$BFISH,"_",temp3$PSU)
length(unique(temp3$Cam_ID))
#724 PSUs



#sum to get catch of spp by psu
temp3[, c(21:266)] <- sapply(temp3[, c(21:266)], as.numeric) #convert catch and presence to numeric

temp4<-cbind(temp3[,21:143],temp3$Cam_ID)
names(temp4)[124]<-"Cam_ID"

psu_sum<-aggregate(temp4[1:123],by=list(temp4$Cam_ID),FUN=sum, na.rm=T) #Group.1 is Cam_ID
dim(psu_sum)
head(psu_sum)#Group.1 is Cam_ID

#mean to get depth, temp, lat, long, and presence (any catch>0),Cam_ID
temp5<-cbind(temp3$OBS_LAT,temp3$OBS_LON,temp3$OFFICIAL_DEPTH_M,temp3$OFFICIAL_TEMP_C,temp3$DROP_TIME_NEW ,temp3[,c(148:267)])
temp5[, c(6:124)] <- sapply(temp5[, c(6:124)], as.numeric) #convert catch and presence to numeric

names(temp5)[1:5]<-c("MEAN_OBS_LAT","MEAN_OBS_LON","MEAN_OFFICIAL_DEPTH_M", "MEAN_OFFICIAL_TEMP_C","MEAN_DROP_TIME_NEW")

#DROP_TIME_NEW reformat
library(chron)


temp5$MEAN_DROP_TIME_NEW<-strptime(temp5$MEAN_DROP_TIME_NEW,  "%H:%M:%S")
temp5$MEAN_DROP_TIME_NEW<-times(strftime(temp5$MEAN_DROP_TIME_NEW,"%H:%M:%S"))

psu_mean<-aggregate(temp5[,1:124],by=list(temp5$Cam_ID),FUN=mean,NA.rm=T)
psu_mean[,c(7:125)]<-ifelse(psu_mean[,c(7:125)]>0,1,0) #convert decimal presences to 1.0

temp6<-merge(psu_mean,psu_sum,by="Group.1") #merge means and sums

#merge back with original dataframe for complete sample data
names(BFISH_CAM_S)
BFISH_CAM_S$Cam_ID<-paste0(BFISH_CAM_S$BFISH,"_",BFISH_CAM_S$PSU)

Camera_PSU_all<-merge(BFISH_CAM_S,temp6,by.x="Cam_ID",by.y="Group.1")

#remove non-unique variables t o reduce dataframe
remove<-c("DROP_CD","DROP_TIME","GEAR_ID","OBS_LAT","OBS_LON","OFFICIAL_DEPTH_M", 
                                  "SEQUENCE_NR", "OFFICIAL_TEMP_C", "DROP_TIME_NEW")  #, "VESSEL", "DROP_DATE"

#25 drop dates for a PSU are different, 3 vessels for PSU are different


Camera_PSU_red<-Camera_PSU_all[,!names(Camera_PSU_all)%in%remove]


Camera_PSU<-unique(Camera_PSU_red)
dim(Camera_PSU) #should be 724

library(dplyr)

Camera_PSU_single<-Camera_PSU %>% distinct(Cam_ID, .keep_all = TRUE) #this keeps one Cam_ID (the first one) for PSUs sampled on 2 different days
                                                                 #not a big deal for catches, and presences, which are already combined by PSU.  Might have small effect on other variables 
                                                                 #in this case, it may be best to align with covariate data at Drop_code level, then run this script to combine within PSU and reduce to unique

#write.csv(Camera_PSU_single,"C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\Camera_PSU.csv")
#Camera_PSU<-read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\Camera_PSU.csv")

#merge with PSU_table file here....
Camera_PSU<-merge(Camera_PSU,PSU_table,by="PSU",all.x=TRUE)
write.csv(Camera_PSU,"C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\BFISH variables\\Camera_PSU.csv")

#now match up with PacIOOS data from Netcdf files
#what are the earliest and latest BFISH dates by year (for NEtCDF download)
tapply(Camera_PSU$DROP_DATE,Camera_PSU$BFISH,range)
tapply(Camera_PSU$MEAN_OBS_LAT,Camera_PSU$BFISH,range)
tapply(Camera_PSU$MEAN_OBS_LON,Camera_PSU$BFISH,range)
#$BFISH_2016_F
#[1] 20161014 20161104

#$BFISH_2017_F
#[1] 20171021 20171118

#$BFISH_2017_S
#[1] 20170309 20170322

#$BFISH_2018_F
#[1] 20180920 20181128

#$BFISH_2019_F
#[1] 20190912 20191126

#$BFISH_2020_F
#[1] 20201008 20201013

