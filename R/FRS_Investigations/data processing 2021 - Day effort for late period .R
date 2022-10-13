rm(list=ls())

fulldata=read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\Raw FRS depredation\\Full_FRS_2021.csv",header=T)
dim(fulldata)

#new data has fields in lower case letters.  Change to upprcase to avoid changing the code.

new_names<-c("LICENSE","TRIP_END","FISHED","BUOY","AREA","SUBAREA","GEAR","SETS","HOURS","NO_GEARS_OR_NET_LENGTH",
             "HOOKS","NET_HEIGHT","BANK","BANK_QUAD","DEPTH_BEG","DEPTH_END","WIND_SPEED","WIND_DIR","WAVE_HT",
             "CUR_SPEED","CUR_DIR","NOCATCH","PORT_LAND","LOGBOOK","NUMBER_LOST_TO_PREDATOR","PREDATOR_CODE",
             "PREDATOR","CHARTERED","R_MONTH","R_YEAR","FORMTYPE","TIMELINK","TRIPID","F_AREA","F_BUOY","F_GEAR",
             "F_AR_GR","F_BY_PK","IMP_USER","LNAME","FNAME","VESSEL","USCG_HA_NO","PORT_DPRTR","TRIP_BEG","SYSDATE",
             "SPECIES","CAUGHT","LBS","NUM_SOLD","LBS_SOLD","VALUE","BUYER","LOST","RELEASE","REL_ALIVE","REL_DEAD",
             "DAMAGE","F_SPEC_GR","F_SPEC_AR","F_SEASON","F_PIECES","F_LBS","F_CTCHLIM","ORIG_USER","APPROVED","REC_NUM",
             "BF_WIND_SP","BF_DEALER","CYEAR","FYEAR")

names(fulldata)<-new_names


#NEED TO READ IN DATA FROM 2018 because it has HA_NO
early=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Raw data\\picdr_112849_fy48_15.csv")
late=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Raw data\\picdr_112849_fy16_18.csv")

fulldata=rbind(early,late) #JS 4/21/20 still getting warning of "invalid factor level, NA generated"
dim(fulldata) 



#######################################################################################################################
#code from "01 workshop data processing
######################################################################################################################

# --Use only records from the MHI------------------------------
#Read in key table from last assessment for MHI areas and keep
#only these grids
mhi_area=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE\\supporting files\\mhi_area.csv",header=T)

#2020 - use of reggie areas does not change the number of obs in nwhi:
#mhi_area=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE\\supporting files\\BF_Area_Grid_Reggie.csv",header=T)
#NOTE: This is not the final MHI definition, but was used in early
#processing of the CPUE data. The proper MHI definition was used
#in "CPUE_data_preparation_WPSAR_update.R)

####----#####
mhidata=fulldata[fulldata$AREA%in%mhi_area$area,]
####----#####

dim(fulldata)-dim(mhidata) #removes records from NWHI
#removes 224109 for 2020 update

# --Separate out deepsea handline gear as well as -------------------------
#other gears used on trips using deepsea handline gear
check_gear=table(mhidata$GEAR,useNA="always")

#Gear 3 records...
temp=mhidata[mhidata$GEAR==3,]
dim(temp)
#...add in records from other gears on "trips" that used gear 3
#How to define a "trip" temporally
#Base "trips" on licence number and fished date

#Determine trip number for trips that used gear 3
addrec=unique(temp[,c("LICENSE","FISHED")])
#assign trip number and merge it to dataset
addrec$trip=c(1:nrow(addrec)) #assign trip number
temp2=merge(addrec,mhidata,by=c("LICENSE","FISHED")) #merge to dataset
dim(temp2)
check_gear2=table(temp2$GEAR,useNA="always") 

####----#####
#Records from gear 3 and other gears on trips with gear 3
#although ultimately we use only gear 3 records
mhidata2=temp2 #1011510 records
#1035335 records for 2020 update
#1057708 in 2021
####----#####

#and remove fishing year 2016 for which data are not complete - JS 2020 this would be 2019 so no need
#leave alone for 2021 exploration
####----#####
mhidata3=mhidata2#[mhidata2$FYEAR<=2018,]
dim(mhidata3) 


####----#####



# --Output final datasets -------------------------
#For group step 3 comparison
#Keep some personal information
#temp_out=subset(mhidata3,select=-c(USCG,HA_NO,IMP_USER,O_FNAME,O_LNAME,O_VESSEL,O_USCG,O_HA_NO,ORIG_USER))
#JS changed fomat of these names for 2020

#2021 USCG_HA_NO is a single field, not USCG and HA_NO separately

#temp_out=subset(mhidata3,select=-c(USCG_HA_NO,IMP_USER,VESSEL,USCG_HA_NO,ORIG_USER))  #we will need FNAME and LNAME


temp_out=mhidata3[mhidata3$GEAR==3,] #871002   records of only gear 3 records
write.csv(temp_out,"C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\Raw FRS depredation\\ForGroupCpueAnalysis_gear3Step3_2020.csv",row.names=F)

rm(list=ls())

###################################################################################################################
#code from "CPUE_data_preparation_2020" 
#####################################################################################################################


#modified by JS for update on 4/21/2020


#######################################

#Purpose: Setup dataset for CPUE analysis
#         for bottomfish assessment. This is 
#         a working document based on
#         decisions made in the data workshops.
#Use: Archives for 2018 Deep 7 MHI bottomfish 
#     assessment
#By: Brian Langseth, PIFSC
#Created: 8/8/2016
#Last updated: 3/21/2018
#Rversion 3.2.4

#Follows: 01_WORKSHOP_data_processing_2020.R (need to have run first)

#Supporting files:
# BF_Area_Grid_Reggie.csv (valid MHI grids)
# port_comport_table.csv (assigns ports to common ports)
# comport_area_distance_table.csv (dist between common ports and grids)
# FI to HDAR grid_priortoFY16.csv (location of fishers for survey)
# MHI_Fishchart2007110807_gridpts.dbf (midpoints of management grids)
# getWind.R (code to assign downloaded wind data to records)
# oceanwindsdly_X.NC where X is for separate files 1 to 28 (wind data) #have 31 files for 2020 update

#Input: Dataset output from 01_WORKSHOP_data_processing.R
#       with MHI areas, gear 3 records, fisher names 
#       and fishing years 1948-2015.

#Output: Best available dataset for CPUE time series for
#        deep7 bottomfish and single species paka models.
#        Datasets are both those with all final records
#        used, and combined into fishing events with a 
#        single value for each factor in the standardization.
#        I added fields: trip, lbs, deep7_lbs, deep7_perc,
#        n_ports, n_areas, uku_lbs, cum_exp,

#######################################

rm(list=ls())

# --Read in basic reduced dataset from 01_WORKSHOP_data processing--------------------------
datafull=read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\Raw FRS depredation\\ForGroupCpueAnalysis_gear3Step3_2020.csv")

dim(datafull)


# --FILTER data: keep duplicates, no invalid areas, no license number 0----------
a=duplicated(datafull)
sum(a) #663  #759 in 2020 update     #668?? (5/19/20)

#2792 in 2021 data...

#-Read in the mhi_areas excel file from Reggie that lists which are
#valid entries and which are not, and which are ponds
#Only keep valid entries
areasReggie=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE\\supporting files\\BF_Area_Grid_Reggie.csv",header=T)
#areasReggie=read.csv("C:\\Users\\Brian.Langseth\\Desktop\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\Data exploration\\Workshop_files_Aug16\\Tuesday analyses\\BF_Area_Grid_Reggie.csv",header=T)
valid=areasReggie[which(areasReggie$Valid.==""),]$area
#Confirm subareas are valid for 16123 (which is a valid area)
table(datafull[!is.na(datafull$SUBAREA),]$AREA) #all subareas are for 16123
table(datafull$SUBAREA) #filter out A and B subareas as well as those without a subarea
####----#####
data=datafull[datafull$AREA%in%valid,] #valid areas
data=data[!data$SUBAREA%in%c("A","B"),] #remove known invalid subareas (16123A and 16123B)
data=data[!(data$AREA==16123 & is.na(data$SUBAREA)),] #remove the 16123 records without a subarea distinction
####----#####

#-Remove all instances where license number is 0 because
#these are unlikely to be from the same trip
head(table(data$LICENSE))
table(data[data$LICENSE==0,]$FYEAR) #all in 1997 or before
####----#####
data=data[data$LICENSE>0,]
####----#####
dim(data) #800130 records
#826450 records in update
#867968 in 2021 data

# --STEP01: select deep7 bottomfish records----------------------------

##-##
#-Determine the percentage of deep7 caught for each trip and add to dataset
#Base this on the larger of lbs_kept and lbs_sold (see catchdata for details)
id=c(15,17,19,21,22,36,58,97) #include both ehu codes
data$lbs=pmax(data$LBS,data$LBS_SOLD,na.rm=T)
alllbs=aggregate(data$lbs,by=list(data$trip),sum)
deep7lbs=aggregate(data[data$SPECIES%in%id,]$lbs,by=list(data[data$SPECIES%in%id,]$trip),sum) #lbs 
#merge records since many trips dont catch deep7
lbs_trip=merge(alllbs,deep7lbs,by="Group.1",all=T)
colnames(lbs_trip)=c("trip","all_lbs","deep7_lbs")
lbs_trip$deep7_perc = round(lbs_trip$deep7_lbs/lbs_trip$all_lbs,3)
#replace NAs with 0
lbs_trip[is.na(lbs_trip$deep7_perc),"deep7_perc"]=0 
lbs_trip[is.na(lbs_trip$deep7_lbs),"deep7_lbs"]=0
dim(lbs_trip) #number of unique trips
####----#####
#Add percentage of deep7 fish caught per trip to database
data1=merge(data,lbs_trip[,c("trip","deep7_perc","all_lbs","deep7_lbs")],by=c("trip"))
####----#####


##-##
#-Remove some trips where no deep7 bottomfish were caught
#Of those trips that did not catch deep7, remove all records from 
#trips that caught PMUS species, Uku, or Species=0 from Sept. 2002 back
#All records w/o deep7
data1old=data1[substr(data1$FISHED,1,7)<="2002-09",]
all0=data1old[data1old$deep7_perc==0,]
table(all0$SPECIES%in%id) #12 records where deep7 were recorded, just not as weight 
#14 in 2021 data

dim(all0)
#Records from trips that caught PMUS
pmus.id=c(1:14,106,107,108,320,321,323,324,39,118) #39 is snake mackeral (Gempylidae) and 118 is monchong (Bramidae)
pmus_all0=all0[all0$SPECIES%in%pmus.id,] #records of PMUS from trips w/o deep7
pmus_all0_trips=unique(pmus_all0$trip) #trips with records of PMUS from trips w/o deep7
pmus_all0_records=all0[all0$trip%in%pmus_all0_trips, ]#all records from trips with PMUS w/o deep7
dim(pmus_all0_records) 
#Records from trips that caught Uku
uku_all0=all0[all0$SPECIES==20,] #records of Uku from trips w/o deep7 
uku_all0_trips=unique(uku_all0$trip) #trips with records of Uku from trips w/o deep7
uku_all0_records=all0[all0$trip%in%uku_all0_trips,] #all records from trips with uku w/o deep7
dim(uku_all0_records)
#Records from trips that had Species=0
spec0_all0=all0[all0$SPECIES==0,] #records of Species=0 from trips w/o deep7 
spec0_all0_trips=unique(spec0_all0$trip) #trips with records of Species=0 from trips w/o deep7
spec0_all0_records=all0[all0$trip%in%spec0_all0_trips,] #all records from trips with Species=0 w/o deep7
dim(spec0_all0_records)
#Total number of records removed from records w/o deep7
dim(all0[all0$trip%in%unique(c(pmus_all0_trips,uku_all0_trips,spec0_all0_trips)),])
####----#####
#Remove trips w/o deep7, and with records of PMUS species, uku, and species=0 <= Sept. 2002
data2=(data1[!(data1$deep7_perc==0 & data1$trip%in%pmus_all0_trips) & 
               !(data1$deep7_perc==0 & data1$trip%in%uku_all0_trips) & 
               !(data1$deep7_perc==0 & data1$trip%in%spec0_all0_trips),])
dim(data2)
####----#####

##-##
#-Kona/South Point: Remove trips that caught any PMUS species while catching <50lbs of deep7 in FY1985 and before
BIarea = c(108,100,101,102,128,120,121,122 ) # DAR grids off Kona and South Point
BIdata = data2[data2$AREA%in%BIarea,] 
#Which of these trips had deep7 catch less than 50 lbs and caught a PMUS in FY1985 and before
BIdata_pmus=BIdata[BIdata$SPECIES%in%pmus.id,]
BIdata_pmus50=BIdata_pmus[BIdata_pmus$deep7_lbs<50,]
BIdata_pmus50_1985=BIdata_pmus50[BIdata_pmus50$FYEAR<=1985,]
####----#####
data3=data2[!data2$trip%in%unique(BIdata_pmus50_1985$trip),]
dim(data3)
####----#####
#2020:
#determine number of single reporting days removed:
length(unique(data1$trip))-length(unique(data3$trip))

##-##
#-Remove records (trips) that occurred when the fishery was closed
#Date fishery closed: Apr 16,2008; Jul 6,2009; Apr 20,2010; Mar 12,2011
#Theres a lot, originally many of these were removed when I included
#trips after Sept. 2002 to be removed that didnt catch any deep7 (but
#did catch PMUS, Uku, and species=0) - most (by weight) of these are uku
data3$FISHED=as.Date(data3$FISHED)
####----#####
data4=data3[!(data3$FISHED>"2008-04-16" & data3$FISHED<="2008-08-31") &
              !(data3$FISHED>"2009-07-06" & data3$FISHED<="2009-08-31") &
              !(data3$FISHED>"2010-04-20" & data3$FISHED<="2010-08-31") &
              !(data3$FISHED>"2011-03-12" & data3$FISHED<="2011-08-31"),]
dim(data4)
####----#####
#determine number of single reporting days removed:
length(unique(data3$trip))-length(unique(data4$trip))

# --STEP02: select records that comply to a unit of effort---------------------
#Single reproting day (CML-FISHED) combo up through Sept. 2002   NEED TO DO for full data set now, so change dates to 2018-06
#Hours fished after Sept. 2002

##-##
#-Multiday trips: distance cutoff
#Read in distance values - These are files provide by Jerry Ault and Steve Smith
#from U of Miami for an analysis there were doing on FRS data. 
#distances=read.csv("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\Data exploration\\Grid distances and FI survey fishers\\Grid distance\\comport_area_distance_table.csv",header=T,col.names=c("com_port","AREA","SUBAREA","distance"))
#ports=read.csv("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\Data exploration\\Grid distances and FI survey fishers\\Grid distance\\port_comport_table.csv",header=T,col.names=c("PORT_LAND","com_port"))

distances=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE\\supporting files\\comport_area_distance_table.csv",header=T,col.names=c("com_port","AREA","SUBAREA","distance"))
ports=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE\\supporting files\\port_comport_table.csv",header=T,col.names=c("PORT_LAND","com_port"))
#There are 15 ports in the data that aren't in Steve and 
#Jerrys distance matrix (entries 0,3,6,10,23,51,51,101,201,330,382,403,427,432)
#Jerry and Steve made changes to three of these entries (I describe these in "step 5 from Ault_Smith\\port_comport_work.xls") 
port_changes=matrix(c(101,182,201,291,427,411),3,2,byrow=T)
colnames(port_changes)=c("PORT_LAND","com_port")
ports2=rbind(port_changes,ports)
####----#####
#Add common port to dataset.
data4b=merge(data4[substr(data4$FISHED,1,7)<="2018-06",],ports2,by="PORT_LAND",all.x=T) 
table(data4b$com_port,useNA="always") #5809 records that dont have a corresponding com_port
####----#####
#Remove distance combination for subareas of area 16123 that aren't in the data
#which are all subareas other than subarea "C"
head(data4b[data4b$AREA==16123,])
distances=distances[!(distances$AREA==16123 & distances$SUBAREA!="C"),]
#This section not needed anymore because reduced 16123 records to only subareas A,B,C
#so now only have 3 records of 16123, all in 682-16123 C
# #Remove area-port combos from distance table that are 512/611/682-16123 area 
# #which are the common ports for the PORT_LAND records that fished area 16123.
# #Do this except for 682-16123 C, 611-16123 B, and 512-16123 Z.
# #Do 611-B because 512-B and 512-Z (no subarea) are similar 
# #and so assume B is best among the options for 611.
# distances[distances$AREA==16123,]
# distances2=distances[-c(790,791,792,793,1295,1297,1435),]
# distances2[distances2$AREA==16123,]

##Assign number of known common ports visited and number of areas visited on a trip
#How many areas per trip are there?
multiarea=rowSums(table(data4b$trip,data4b$AREA)>0) 
table(multiarea) #at most 5 areas per trip ##  Looks like some 6 and a 7 in 2020 update
#How many common ports per trip are there?
multiport=rowSums(table(data4b$trip,data4b$com_port)>0)
table(multiport,useNA="always") #161 trips vist 2 common ports (0 are because a common port isn't available for PORT_LAND)  #161 in update
ports_v=data.frame("trip"=names(multiport),"n_com_ports"=multiport)
areas_v=data.frame("trip"=names(multiarea),"n_areas"=multiarea)
ports_areas=merge(ports_v,areas_v,by="trip",sort=FALSE)
####----#####
#number of known common ports and areas per trip
data4c=merge(data4b,ports_areas,by="trip")
####----#####
#There are 11 records with an unknown com_port but these are on trips
#with a known com_port (because their port_land is different than
#the other records within the trip) hence there are 5798 (not 5809) 
#0 com_port records.
table(data4c$n_com_ports,data4c$n_areas)
table(unique(data4c[,c("trip","n_com_ports","n_areas")])$n_com_ports,unique(data4c[,c("trip","n_com_ports","n_areas")])$n_areas)
#There are also 10 records in 1 area trips with 2 ports with 
#a known com_port but one that is not in the database 
#from Jerry and Steve and therefore will have an NA distance
#See code run after data5 is set up

##Assign distances based on combination of ports-areas
#FOR TRIPS WITH 1 AREA
data4c_1area=merge(data4c[data4c$n_areas==1,],distances[,-3],by=c("com_port","AREA"),all.x=T)
#Those with 1 area and 1 common port. Distances are valid.
data4c_1area_1=data4c_1area[data4c_1area$n_com_ports==1,]
#
#Those with 1 area and 2 common ports. 
#Need to separate trips when applicable and use a single distance
data4c_1area_2=data4c_1area[data4c_1area$n_com_ports>1,]
#Separate trips when applicable:
#For ports on same island or with general port code on same island, keep together as 1 trip.
#For ports on distant islands, separate. For ports with very distant 
#island, keep together because likely a database error (define as greater than 3 islands away,
#so Kauai-Maui, or Oahu-Big Island).
u_areas=tapply(data4c_1area_2$com_port,data4c_1area_2$trip,unique)
u_areas2=matrix(as.numeric(unlist(u_areas,use.names=F)),length(u_areas),2,byrow=T,dimnames=list(rownames(u_areas)))
sametrips=rownames(u_areas2)[which(substr(u_areas2[,1],1,1)==substr(u_areas2[,2],1,1))]
difftrips=rownames(u_areas2)[which(substr(u_areas2[,1],1,1)!=substr(u_areas2[,2],1,1))]
#Of the difftrips, which are very far apart (Big Island - Oahu or Maui - Kauai) and thus assumed to be errors. Want to keep together
far_difftrips=rownames(u_areas2)[which(abs(u_areas2[,1]-u_areas2[,2])>300)]
difftrips2=difftrips[!difftrips%in%far_difftrips]
#Take the trips that are to be separate, reclassify trip and add back to database
tosep=data4c_1area_2[data4c_1area_2$trip%in%difftrips2,]
tosep$trip2=paste0(tosep$trip,".",floor(tosep$com_port/100)) #to add a decimal that indicates the first digit of the com_port (needed to distinguish ports 499 and 501, which weren't when rounding)
data4c_1area_2sep=merge(data4c_1area_2,unique(tosep[,c("trip","com_port","trip2")]),by=c("trip","com_port"),all.x=T)
data4c_1area_2sep[is.na(data4c_1area_2sep$trip2),"trip2"]=data4c_1area_2sep[is.na(data4c_1area_2sep$trip2),"trip"]
#Reclassify distance. 
#For each trip2 the most frequent distance is used, which can include NAs.
#Create own function to extract the modal distance. For instances where the number of distances 
#is the same (e.g. two record trip that goes to two com_ports, this function defaults to the shorter distance)
modal=function(x){
  aa=as.numeric(names(which.max(table(x,useNA="always"))))
  return(aa)}
temp=aggregate(data4c_1area_2sep$distance,by=list(data4c_1area_2sep$trip2),modal)
colnames(temp)=c("trip2","distance2")
data4c_1area_2FINAL=merge(data4c_1area_2sep,temp,by="trip2")

#FOR TRIPS WITH 2+ AREA
data4c_manyareas=merge(data4c[data4c$n_areas>1,],distances[,-3],by=c("com_port","AREA"),all.x=T)
#Those with many areas and 1 common port
data4c_manyareas_1=data4c_manyareas[data4c_manyareas$n_com_ports==1,]
head(data4c_manyareas_1[order(data4c_manyareas_1$trip),],20)
#Reclassify distance.
#Take greatest distance among areas as the trip distance
temp=aggregate(data4c_manyareas_1$distance,by=list(data4c_manyareas_1$trip),max)
colnames(temp)=c("trip","distance2")
data4c_manyareas_1FINAL=merge(data4c_manyareas_1,temp,by="trip")
#
#Those with many areas and 2 common ports
#Need to separate trips when applicable and use a single distance
data4c_manyareas_2=data4c_manyareas[data4c_manyareas$n_com_ports>1,]
head(data4c_manyareas_2[order(data4c_manyareas_2$trip),],20)
#Separate by trip. Nearly all have distinct areas for distinct ports
#and vice versa. Hence use trip-com_port as unique trip identifier for all
#Take the trips that are to be separate, reclassify trip and add back to database
data4c_manyareas_2$trip2=paste0(data4c_manyareas_2$trip,".",data4c_manyareas_2$com_port) #use entire com_port to separate out trips
#Reclassify distance. 
#Take max distance among areas as the trip distance
temp=aggregate(data4c_manyareas_2$distance,by=list(data4c_manyareas_2$trip2),max)
colnames(temp)=c("trip2","distance2")
data4c_manyareas_2FINAL=merge(data4c_manyareas_2,temp,by="trip2")

#FOR TRIPS WITH 0 COM_PORT
#No need to add distances since com_port is unavailable
data4c_0port=data4c[data4c$n_com_ports==0,]

##Combine all trips together into a single database
#This includes the records without com_ports 
#Ensure all sub datasets have all names included
data4c_1area_1$trip2=data4c_1area_1$trip #Add trip2 to 1port-1area dataset
data4c_1area_1$distance2=data4c_1area_1$distance #Add distance2 to 1port-1area dataset
data4c_manyareas_1FINAL$trip2=data4c_manyareas_1FINAL$trip #Add trip2 to 1port-manyarea dataset
data4c_0port$distance=data4c_0port$distance2=NA #Add NA distance and distance2 to 0port dataset
data4c_0port$trip2=data4c_0port$trip #Add trip2 to 0port dataset
##Single database
data4c_corrected=rbind(data4c_1area_1,data4c_1area_2FINAL,data4c_manyareas_1FINAL,data4c_manyareas_2FINAL,data4c_0port)
#
#The 10 records in 1 area-2 ports trips with com_port
#not in database
table(data4c_corrected[is.na(data4c_corrected$distance2),]$n_com_ports)
#So NA's in distance are either due to a PORT_LAND that
#doesnt have a corresponding com_port, or a com_port/Area
#combination not in the database

##Now determine distance cutoff
#Corrected distances BY TRIP2 in 10yr increments
data4c_corrected_trip2=unique(data4c_corrected[,c("trip2","FYEAR","distance2")])

hist(data4c_corrected_trip2[data4c_corrected_trip2$FYEAR<1960,]$distance2,breaks=50,main="<1960",xlim=c(0,300),xlab="Distances by trip (nm)")
hist(data4c_corrected_trip2[data4c_corrected_trip2$FYEAR>=1960 & data4c_corrected_trip2$FYEAR<1970,]$distance2,breaks=40,main="1960-1969",xlim=c(0,300),xlab="Distances by trip (nm)")
hist(data4c_corrected_trip2[data4c_corrected_trip2$FYEAR>=1970 & data4c_corrected_trip2$FYEAR<1980,]$distance2,breaks=50,main="1970-1979",xlim=c(0,300),xlab="Distances by trip (nm)")
hist(data4c_corrected_trip2[data4c_corrected_trip2$FYEAR>=1980 & data4c_corrected_trip2$FYEAR<1990,]$distance2,breaks=50,main="1980-1989",xlim=c(0,300),xlab="Distances by trip (nm)")
hist(data4c_corrected_trip2[data4c_corrected_trip2$FYEAR>=1990 & data4c_corrected_trip2$FYEAR<2000,]$distance2,breaks=50,main="1990-1999",xlim=c(0,300),xlab="Distances by trip (nm)")
hist(data4c_corrected_trip2[data4c_corrected_trip2$FYEAR>=2000 & data4c_corrected_trip2$FYEAR<2003,]$distance2,breaks=50,main="2000-2002",xlim=c(0,300),xlab="Distances by trip (nm)")

hist(data4c_corrected_trip2[data4c_corrected_trip2$FYEAR<1955,]$distance2,breaks=50,main="<1955",xlim=c(0,300),xlab="Distances by trip (nm)")
hist(data4c_corrected_trip2[data4c_corrected_trip2$FYEAR>=1955 & data4c_corrected_trip2$FYEAR<1960,]$distance2,breaks=50,main="1955-1959",xlim=c(0,300),xlab="Distances by trip (nm)")
#Looks like 30nm is good cutoff in general. 20nm could be used for <1960
#as 15nm seems too small but 30nm is good too. 
#Use 30nm for cutoff after 1960. Its inclusive of more data, and 
#minimizes the amount of changes we make. There also is a break
#in the histograms at this value.
####----#####
data4c_corrected$days=ceiling(data4c_corrected$distance2/30)
####----#####

##Relink back to original database (which include 
#the records after Sept. 2002)
newdata4=data4[substr(data4$FISHED,1,7)>"2018-06",] #there are no data here, so just rename to data 5 below
newdata4$com_port=newdata4$n_com_ports=newdata4$n_areas=newdata4$distance=newdata4$distance2=newdata4$days=NA
newdata4$trip2=newdata4$trip
####----#####
#data5=rbind(data4c_corrected,newdata4)

data5<-data4c_corrected
####----#####
#Number of records with more than 1 day
table(data5$days,useNA="always")
table(data5$days>1,useNA="always")
#Number of trips with more than 1 day
table(unique(data5[,c("trip2","days")])$days)
table(unique(data5[,c("trip2","days")])$days>1,useNA="always")

##-##
#-Multiday trips: only first-last reports
#Remove trips from fishers who only ever reported on the first or last day of the month
#Use only for records prior to September 2002  CHANGE to 2018
monthday=table(substr(data5[substr(data5$FISHED,1,7)<="2018-06",]$FISHED,6,10))
lastindex=c(31,59,91,121,152,182,213,244,274,305,335,366) #end of each month
monthday[lastindex] #check - exclude Feb. 29
firstindex=c(1,32,61,92,122,153,183,214,245,275,306,336) #first of each month
monthday[firstindex] #check
choices=c(names(monthday[lastindex]),names(monthday[firstindex])) #first and last dates of the month

#Explore the extent to which a fisher only reports on the first or last day of the month
past_data5=unique(data5[substr(data5$FISHED,1,7)<="2018-06",c("trip2","LICENSE","FISHED","deep7_lbs","FYEAR")]) #trips in past
past_data5$firstlast=0 #indicator if trips was first/last day of month
past_data5[which(substr(past_data5$FISHED,6,10)%in%choices),"firstlast"]=1 
tab=table(past_data5$LICENSE,past_data5$FYEAR,past_data5$firstlast) #table of trips by fishing year and license that occurred and did not occur on the first and last days of the month
per_tab=round(tab[,,2]/(tab[,,1]+tab[,,2]),2) #percent of trips by fisher-year that occurred on the first or last day of the month

#Pull out the 646 that were 100%; only occurred on the first or last day of the month 100%
col=floor(which(per_tab==1)/nrow(per_tab)) #column occurring
row=which(per_tab==1)-nrow(per_tab)*col #row occurring
all_firstlast=data.frame("LICENSE"=NA,"FYEAR"=NA,"percent"=NA)
for(i in 1:length(which(per_tab==1))){
  all_firstlast[i,"LICENSE"]=as.numeric(dimnames(per_tab)[[1]][row[i]]) #add license # (row names)
  all_firstlast[i,"FYEAR"]=as.numeric(dimnames(per_tab)[[2]][col[i]+1]) #add fyear (col name). Add one to col due to flooring
  all_firstlast[i,"percent"]=per_tab[row[i],col[i]+1] #add value
}
#merge to obtain trip numbers of license-fished combinations
trips_firstlast=merge(all_firstlast,past_data5,by=c("FYEAR","LICENSE")) #trips on the first and last days of the month from CML-fished combinations that did not record another trip
#merge with original data to obtain all records from license-fished combinations
records_firstlast=merge(data5[substr(data5$FISHED,1,7)<="2018-06",],trips_firstlast,by=c("trip2","LICENSE","FISHED","deep7_lbs","FYEAR"))

####----#####
#exclude records from the first or last day of the month in a year
data6=data5[!data5$trip2%in%trips_firstlast$trip2,] 
dim(data6)
####----#####


# --STEP03: select data representative of CPUE trends---------------------

##-##
#-Representative data: remove records from fishers that never caught weight of Deep7
####----#####
#no punctuation in 2020 data but this works
data6$fisher=gsub("[.]","",tolower(paste(data6$FNAME,data6$LNAME))) #remove punctuation and set all names to lowercase
####----#####
checkFISH=table(data6$fisher,data6$deep7_lbs==0)
per_checkFISH=checkFISH[,1]/rowSums(checkFISH) #percent of recors with deep7
checkFISH_dataframe=data.frame("fisher"=attr(checkFISH,"dimnames")[[1]],"records"=as.numeric(rowSums(checkFISH)),"perc_deep7"=as.numeric(per_checkFISH))
neverDeep7=checkFISH_dataframe[checkFISH_dataframe$perc_deep7==0,]
dim(neverDeep7) #1482 people never caught a Deep7  #1461 in 2021 data
sum(neverDeep7$records) #have 7851 records  #8390 in 2021 data
####----#####
data7=data6[!data6$fisher%in%neverDeep7$fisher,]
####----#####

##-##
#-Remove fishing independent survey trips 
#Read in database of fisher, day, areas from the FI surveys #NOT UPDATE FOR 2021 DATA PULL
#FI=read.csv("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\Data exploration\\Grid distances and FI survey fishers\\FI survey fishers\\FI to HDAR grid_priortoFY16.csv",header=T)
FI=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE\\supporting files\\FI to HDAR grid_ALLFY.csv",header=T)
FI$Date=as.Date(FI$Date,format="%m/%d/%Y")
FI$last=tolower(FI$last)
#Find FI survey trips in filtered database
FI$combo=paste(FI$Date,FI$HDAR.Grid,FI$last) #do this to combine all unique conditions> If you cond1 & cond2.... then dont get unique conditions
a=data7[paste(data7$FISHED,data7$AREA,tolower(data7$LNAME))%in%FI$combo,]
table(a$fisher) #insure that use of last name was appropriate and only got Captains from survey (it was)
#Check to ensure that a trip doesn't cover multiple areas (it doesn't)
dim(a[!duplicated(a[,c("FISHED","AREA","fisher")]),])
length(unique(a$trip2))

####----#####
data8=data7[!data7$trip2%in%unique(a$trip2),] #97 records, 40 trips removed
####----#####


# --STEP04: Factors for CPUE standardization---------------------

##-##
#-Add a measure of experience to each trip; cumulative number of trips
red=data8[,c("trip2","fisher","FISHED")]
ord_uniq=unique(red)[order(unique(red)$fisher,unique(red)$FISHED),] #unique fisher-trip combinations ordered by fisher
library(plyr)
num_exp=count(ord_uniq,vars="fisher") #number of unique fisher-trip combos by fisher
ord_uniq$cum_exp=unlist(lapply(num_exp$freq,FUN=seq)) #add experience value; from 1 to number of unique fisher-trip combinations per fisher 
red2=join(red,ord_uniq,by=c("trip2","FISHED","fisher")) #use join so that the same order is maintained with data6
####----#####
data8$cum_exp=red2$cum_exp
####----#####

##-##
#Determine amount of Uku caught in each trip and add to dataset
ukulbs=aggregate(data8[data8$SPECIES==20,]$lbs,by=list("trip2"=data8[data8$SPECIES==20,]$trip2),sum)
#merge records since many trips dont catch deep7
####----#####
data9=merge(data8,ukulbs,by="trip2",all=TRUE)
data9[is.na(data9$x),"x"]=0
dimnames(data9)[[2]][88]="uku_lbs" #dim was previously 82
####----#####

##-##
#Add quarter to the dataset
#Base it on FYEAR quarter
months=as.numeric(substr(data9$FISHED,6,7))
season=months
season[months%in%c(1:3)]=3
season[months%in%c(4:6)]=4
season[months%in%c(7:9)]=1
season[months%in%c(10:12)]=2
####----#####
data9$qtr=season
####----#####


#SKIP WIND FOR 2021 DATA - skip to line 654

##-##
#Add wind speed and wind direction information to appropriate records. 
#Wind is on a daily basis over 0.25 degree grids starting 1987-07-09
#Use getWind function which is in getWind.R
#-#
#Read in locations of management grids
library(foreign)
#area_loc=read.dbf("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\Data exploration\\Grid distances and FI survey fishers\\Grid distance\\step 5 from Ault_Smith\\ArcGIS_shapefiles\\Fishchart2007110807_gridpts.dbf",as.is=T)
area_loc=read.dbf("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE\\supporting files\\Fishchart2007110807_gridpts.dbf",as.is=T)
area_locMHI=area_loc[area_loc$AREA_ID%in%valid,] #valid areas (from filtering above)
#write.dbf(area_locMHI,"D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\Data exploration\\Grid distances and FI survey fishers\\Grid distance\\step 5 from Ault_Smith\\ArcGIS_shapefiles\\MHI_Fishchart2007110807_gridpts.dbf")
#There are duplicate area values in Jerry and Steve's dbf file for
#area_id's 307,313,323,401,406,and 16123. Since data for only 16123 
#subarea C exists in the dataset, can remove all other subareas for 
#16123. For all others remove the records with small length components
#which are not actually midpoitns of grids.
area_locMHI2=area_locMHI[-c(232,233,243:248, #16123
                            282,290,292,136,140,267),] #307,313(two of them),323,401,406

#Assign wind information to all records within a trip
#Accounts for trips with more than one area.
#Not all locations will have viable wind data
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#This takes 4-5 hours! <<<<<<<<<<<<<<<<<<<<<<< NOTE <<<<<<<<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#source("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\Data exploration\\winddata\\getWind.R")
source("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE\\supporting files\\getWind_2020.R")
temp=data9[data9$FISHED>="1987-07-09",]
temp[,c("speed","xdir","ydir")]=NA
u_trips2=unique(temp$trip2)
for(i in 1:length(u_trips2)){
  temp_data=temp[temp$trip2%in%u_trips2[i],]
  wind=getWind(temp_data,area_locMHI2)
  temp[temp$trip2%in%u_trips2[i],c("speed","xdir","ydir")]=
    wind[match(temp[temp$trip2%in%u_trips2[i],c("AREA")],wind$AREA),c("speed","xdir","ydir")]
  print(paste(i,"out of",length(u_trips2)))
  #  print(wind[match(temp[temp$trip2%in%u_trips2[i],c("AREA")],wind$AREA),c("speed","xdir","ydir")])
  #  print(wind)
}

#2020 - the code above inserts Area into columns "speed","xdir","ydir" when there is a gap in 
#temporal coverage with 2014-2018 wind files.  This appears to be a funtion of the start/end 
#dates of the oceanwindsdly files, which were possibly misspecified.  Corrections were made as 
#possible, but still some cases remain.  The following code replaces these values with "NAs"

for (i in 1:nrow(temp)){
  temp$speed[i]<-ifelse(temp$speed[i]==temp$area[i],temp$speed[i]<-NA,temp$speed[i]<-temp$speed[i])
  temp$xdir[i]<-ifelse(temp$xdir[i]==temp$area[i],temp$xdir[i]<-NA,temp$xdir[i]<-temp$xdir[i])
  temp$ydir[i]<-ifelse(temp$ydir[i]==temp$area[i],temp$ydir[i]<-NA,temp$ydir[i]<-temp$ydir[i])
  # temp$winddeg[i]<-ifelse(temp$speed[i]==temp$area[i],TP2$winddeg[i]<-NA,TP2$winddeg[i]<-TP2$winddeg[i])
  # temp$winddir[i]<-ifelse(temp$speed[i]==temp$area[i],TP2$winddir[i]<-NA,TP2$winddir[i]<-TP2$winddir[i])
}

#temp$winddir = cut(temp$winddeg, breaks = seq(0,360,45), labels = c("W","NW","N","NE","E","SE","S","SW"))



#Relink with old (prior to wind data) records and form updated dataset
data9_old=data9[data9$FISHED<"1987-07-09",]
data9_old[,c("speed","xdir","ydir")]=NA
####----#####
data10=rbind(data9_old,temp)
####----#####

#2021 data pull - skipped wind

data10=data9


# --OUTPUT: finalized record dataset for CPUE standardization---------------------
#write.csv(data10,"D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\UptoDate files\\Finalized_recordCPUE_dataset.csv",row.names=F)
#write.csv(data10,"C:\\Users\\John.Syslo\\Documents\\2020 Deep 7\\2020 Data prep files\\CPUE\\Finalized_recordCPUE_dataset.csv",row.names=F)
write.csv(data10,"C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\Raw FRS depredation\\Finalized_recordCPUE_dataset.csv",row.names=F)

# --POST PRCOESSING: Establish trip-based dataset for CPUE standardization---------------
rm(list=ls())
#Read in data, set up fields I work with later.
#fulldata=read.csv("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\UptoDate files\\Finalized_recordCPUE_dataset.csv",header=T)
#THIS IS 2020 HERE: 
#fulldata=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE\\Finalized_recordCPUE_dataset.csv",header=T)
fulldata=read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\Raw FRS depredation\\Finalized_recordCPUE_dataset.csv",header=T)

dim(fulldata)
#Prior filtering did not actually account for 2016-2019 FI survey trips before running get_wind. Do that process again here
#-Remove fishing independent survey trips 
#Read in database of fisher, day, areas from the FI surveys
#FI=read.csv("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\Data exploration\\Grid distances and FI survey fishers\\FI survey fishers\\FI to HDAR grid_priortoFY16.csv",header=T)
FI=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE\\supporting files\\FI to HDAR grid_ALLFY.csv",header=T)
FI$Date=as.Date(FI$Date,format="%m/%d/%Y")
FI$last=tolower(FI$last)
#Find FI survey trips in filtered database
FI$combo=paste(FI$Date,FI$HDAR.Grid,FI$last) #do this to combine all unique conditions> If you cond1 & cond2.... then dont get unique conditions
a=fulldata[paste(fulldata$FISHED,fulldata$AREA,tolower(fulldata$LNAME))%in%FI$combo,]
table(a$fisher) #insure that use of last name was appropriate and only got Captains from survey (it was)
#Check to ensure that a trip doesn't cover multiple areas (it doesn't)
dim(a[!duplicated(a[,c("FISHED","AREA","fisher")]),])
length(unique(a$trip2))
####----#####
data8=fulldata[!fulldata$trip2%in%unique(a$trip2),] #191 records removed
dim(data8)


#
#THIS is to bring in 2017, finish filteirng, and eventually compare nominal:
#fulldata=read.csv("C:\\Users\\John.Syslo\\Documents\\2017 Deep 7 BSP\\Raw Data\\Finalized_recordCPUE_dataset.csv",header=T)

#not sure where data8 name came from, has same dim as fulldata for old data right now, so keep using fulldata



fulldata$FISHED=as.Date(fulldata$FISHED)
fulldata$fisher=as.character(fulldata$fisher)

#Set fisher names without a first name to just the last name, and without a last name to NA
#Otherwise the fisher name has NA in it, or its name is "NA NA"
fulldata[is.na(fulldata$FNAME),]$fisher=gsub("[.]","",tolower(fulldata[is.na(fulldata$FNAME),]$LNAME)) #only last name for records with no first name
fulldata[is.na(fulldata$LNAME),]$fisher=NA #any records without a last name, have no first name

#Reduce number of columns to those used as factors in standardization

#Dont use SUBAREA since only apply to 16123 (all are 16123-C)  #JS comment out this subsetting
tdata=fulldata#[,c("FYEAR","FISHED","trip2","AREA","n_areas","HOURS","days","SPECIES","lbs","deep7_lbs","fisher","cum_exp","qtr","uku_lbs", "NO_GEARS_OR_NET_LENGTH", "NUMBER_LOST_TO_PREDATOR", "PREDATOR_CODE","PREDATOR"  )]   #,"speed","xdir","ydir")] no wind for 2021

names(tdata)[82]="n_areas_e" #n_areas prior to Oct. 2002 (early)

##For single species analysis
#tdata=fulldata[,c("FYEAR","FISHED","trip2","AREA","n_areas","HOURS","days","SPECIES","lbs","deep7_lbs","fisher","LICENSE","cum_exp","qtr","uku_lbs","speed","xdir","ydir")]
#names(tdata)[5]="n_areas_e" #n_areas prior to Oct. 2002 (early)

##-##
#Correct hours in recent records (>= Oct. 2002) on trips
#that have more than one hour value recorded, or that 
#have 0 hrs reported.
##-##
#Initially I summed together unique hours but this is wrong, and biases CPUE.
#CPUEs should be treated separatly (thus averaged) for unique hours in unique
#areas, which is same as treating as its own fishing event. So, treat as own trip
#Sum together weight of deep7 for records with the same hours (use that hour for 
#effort) within the same area, and assign new trip category (trip3). Calculte
#new total catch of deep7 within these trips, and assign new catch category.
#
#Determine multihour trips
recent_tdata=tdata[tdata$FISHED>="2002-10-01",]
temp=rowSums(table(recent_tdata$trip2,recent_tdata$HOURS)>0)
table(temp) #782 trips with multiple hours reported [799+44+3 in 2020 update]
multihour_trips=names(which(temp>1))
multihour_red=tdata[tdata$trip2%in%multihour_trips,c("trip2","lbs","SPECIES","deep7_lbs","HOURS","AREA"),]
head(multihour_red[order(multihour_red$trip2),])
#Sum lbs of deep7 and non-deep7 with similar hours within area/trip 
id=c(15,17,19,21,22,36,58,97) #include both ehu codes
lbs_hour_area_trip=aggregate(multihour_red$lbs,by=list(multihour_red$trip2,multihour_red$AREA,multihour_red$HOURS,multihour_red$SPECIES%in%id),FUN=sum)
lbs_hour_area_trip$d7=lbs_hour_area_trip$Group.4*lbs_hour_area_trip$x 
head(lbs_hour_area_trip[order(lbs_hour_area_trip$Group.1),],20)
#Sum d7 lbs across similar hours within an area/trip 
d7lbs_hour_area_trip=aggregate(lbs_hour_area_trip$d7,by=list(lbs_hour_area_trip$Group.1,lbs_hour_area_trip$Group.2,lbs_hour_area_trip$Group.3),FUN=sum)
head(d7lbs_hour_area_trip[order(d7lbs_hour_area_trip$Group.1),],20)
names(d7lbs_hour_area_trip)=c("trip2","AREA","HOURS","d7lbs_mhcorrect")
#Assign new trip definition (trip3) to these fishing events
ordered_d7lhat=d7lbs_hour_area_trip[order(d7lbs_hour_area_trip$trip2),]
ordered_d7lhat$tripADD=unlist(apply(table(d7lbs_hour_area_trip$trip2),MARGIN=1,FUN=seq)) #these are already in trip order
ordered_d7lhat$trip3=ordered_d7lhat$trip2+ordered_d7lhat$tripADD/10
####----#####
tdata2=merge(tdata,ordered_d7lhat[,-5],by=c("trip2","AREA","HOURS"),all=T)
tdata2[is.na(tdata2$trip3),"trip3"]=tdata2[is.na(tdata2$trip3),"trip2"] #non multihour trips have the same trip defintion (trip2)
####----#####
#
#Remove trips with 0 hours reported. Thus for recent trips. 
#There are 97 trips (267 records) with 0 hours recorded #2020- 99 and 273
recent_tdata2=tdata2[tdata2$FISHED>="2002-10-01",]
trips_0hr=unique(recent_tdata2[which(recent_tdata2$HOURS==0),]$trip3)
####----#####
tdata3=tdata2[-which(tdata2$trip3%in%trips_0hr),] #doesnt work with 2021 data
####----#####
#assign tdata3 for 2021 data
tdata3=tdata2


#Assign final deep7 catch amount. Since d7lbs_mhcorrect
#for multihour trips is always less than deep7_lbs, 
#apply minimum of these to the determine the final catch 
#amount (d7catch)
####----#####
tdata3$d7catch=pmin(tdata3$deep7_lbs,tdata3$d7lbs_mhcorrect,na.rm=T)
####----#####
#
#Assigned final effort metric for recent data. So that 
#trips with multiple hours prior to 2002-10-01 are 
#appropriately counted as dayseffort
####----#####
tdata3$houreffort=NA
tdata3[tdata3$FISHED>="2002-10-01",]$houreffort=tdata3[tdata3$FISHED>="2002-10-01",]$HOURS

##-##
#Assign final effort metric (DAYeffort) for early data.
#Effort metric is days, although some days have NA. Some are 
#NAs within a trip with days assigned (happens when the 
#port is misreported - 11 trips), but most have all records in a
#trip as NA. For the later, assume they are one day (cant
#treat them as more than one single-reporting-day). 
##-##
old_tdata3=tdata3[tdata3$FISHED<"2002-10-01",]
nadays=old_tdata3[which(is.na(old_tdata3$days)),]
nadays_trips=unique(nadays$trip3)
nadays_records=old_tdata3[old_tdata3$trip3%in%nadays_trips,]
maxdays=aggregate(nadays_records$days,by=list(nadays_records$trip3),FUN=unique)
names(maxdays)[1]="trip3"
maxdays$nvals=unlist(lapply(maxdays$x,FUN=length)) #which records have more than just NA
maxdays$days_nacorrect=maxdays$nvals #trips with only NA are now 1. 
temp=unlist(lapply(maxdays[which(maxdays$nvals>1),"x"],FUN=max,na.rm=T)) #take max of trips with NA and another value
maxdays[maxdays$trip3%in%maxdays[which(maxdays$nvals>1),"trip3"],]$days_nacorrect=temp #reassign DAYS to be max of NA and the other value
table(maxdays$days_nacorrect)
####----#####
tdata4=merge(tdata3,maxdays[,-c(2,3)],by="trip3",all=T)
tdata4$dayeffort=tdata4$days #set dayeffort to days first
tdata4[which(tdata4$days_nacorrect>0),]$dayeffort=tdata4[which(tdata4$days_nacorrect>0),]$days_nacorrect #correct NA records (assign 1 if all records in trip are NA, assign value if some records have nonNA value)
####----#####

##-##
#Assign area to trip for multiarea trips. Use
#area with most amount of deep7 catch, or if not
#distinguishable use lower numbered grid. Can use
#deep7_lbs here because no multihour trips with 
#corrected weight now have multiple areas
##-##
#First add factor indicating number of areas for recent data
#to help calculations later on. Combine recent and early
#number of areas into a single facotr for number of areas
recent_tdata4=tdata4[tdata4$FISHED>="2002-10-01",]
areas_by_trip_recent=rowSums(table(recent_tdata4$trip3,recent_tdata4$AREA)>0)
table(areas_by_trip_recent) #93 trips with multiple areas in recent time period [102 in 2020 update]
nareas_recent=data.frame("trip3"=names(areas_by_trip_recent),"n_areas_r"=areas_by_trip_recent)
####----#####
tdata5=merge(tdata4,nareas_recent,by="trip3",all=T)
tdata5$n_areas_all=pmax(tdata5$n_areas_e,tdata5$n_areas_r,na.rm=TRUE) #combine n_areas_e with n_areas_r
####----#####
#
#Now assign areas to trips for multiarea trips
table(tdata5$n_areas_all,useNA="always") #no NAs here
all_multiareas=tdata5[tdata5$n_areas_all>1,]
length(unique(all_multiareas$trip3))
wbyarea=aggregate(all_multiareas$lbs,by=list(all_multiareas$trip3,all_multiareas$AREA,all_multiareas$SPECIES%in%id),FUN=sum) #weight by area
length(unique(all_multiareas$trip3)) #1877 trips with multiple areas over all years [1893 in 2020 update]
names(wbyarea)=c("trip3","AREA","deep7","weight")
wbyarea$d7weight=wbyarea$deep7*wbyarea$weight  #weight by area of only deep7 catch
#determine maximum value of deep7 weight for each trip
max_per_trip=aggregate(wbyarea$d7weight,by=list(wbyarea$trip3),FUN=max)
#and assign area correspondnig to that maximum value (this has additional rows because some trips have more than one area with maximum deep7 weight
area_of_max=merge(max_per_trip,wbyarea[,-c(3,4)],by.x=c("Group.1","x"),by.y=c("trip3","d7weight")) #area(s) corresponding to maximum value
names(area_of_max)=c("trip3","d7weight","maxarea")
multiareas=merge(wbyarea,area_of_max[,-2],by=c("trip3"),all=T) #check to ensure accuracy
#
#There are 137 trips with the same weight of deep7 in multiple areas. 
#This includes multiarea trips that did not catch any deep7 (weight of 0).
#Looking at number of deep7 records reported for each tied area doesn't 
#separate them further. 
tied_trips=names(which(table(area_of_max$trip3)>1))
length(tied_trips) #137 in 2020 update  159 n 2021 data
length(names(which(table(area_of_max[area_of_max$d7weight==0,]$trip3)>1))) #number of multiarea trips without deep7 catch
#So choose arbitrary area. Go with lower numbered area: nearer to shore.
sorted_area_of_max=area_of_max[order(area_of_max$trip3,area_of_max$maxarea),] #so that lowest numbered areas within a trip come first
final_key=sorted_area_of_max[!duplicated(sorted_area_of_max[,-3]),] #remove the duplicates (higher numbered areas) when weight is the same between areas within a trip
final_key[final_key$trip3%in%tied_trips,]
names(final_key)[3]="AREA_macorrect"
####----#####
tdata6=merge(tdata5,final_key[,-2],by="trip3",all=T) #correct areas for multiarea trips
tdata6$area=NA #assign AREA and AREA_macorrect to a single area variable (area)
tdata6[is.na(tdata6$AREA_macorrect),"area"]=tdata6[is.na(tdata6$AREA_macorrect),"AREA"] 
tdata6[!is.na(tdata6$AREA_macorrect),]$area=tdata6[!is.na(tdata6$AREA_macorrect),"AREA_macorrect"]
####----#####

##-##
#Calculate CPUE. For records < "2002-10-01" effort is number 
#of single-reporting days (dayeffort), for records >= "2002-10-01" 
#effort is reported hours (houreffort). Catch is lbs of deep7 (d7catch).
##-##
####----#####
tdata6$cpue=NA
tdata6[tdata6$FISHED<"2002-10-01","cpue"]=tdata6[tdata6$FISHED<"2002-10-01","d7catch"]/tdata6[tdata6$FISHED<"2002-10-01","dayeffort"]
tdata6[tdata6$FISHED>="2002-10-01","cpue"]=tdata6[tdata6$FISHED>="2002-10-01","d7catch"]/tdata6[tdata6$FISHED>="2002-10-01","houreffort"]
####----#####
write.csv(tdata6,"C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\Raw FRS depredation\\tdata6_2021.csv",row.names=F) 


#SKIP to line 956 for complex-only data set

# --STEPHENS and MACCALL 2004 method: For classifying trips targeting paka-------
######################
#Set up for method
######################
#Assign variable (yj) that is 1 if on trip j, a paka was caught, and 0 otherwise
sm_data = tdata6
sm_data[sm_data$SPECIES==36,"SPECIES"]=21 #Set ehu code 36 to new code 21
sm_data[sm_data$SPECIES==70,"SPECIES"]=24 #Set weke nono code 70 to new code 24
pakatrips=unique(sm_data[sm_data$SPECIES==19,"trip3"]) #114287 trips with paka
sm_data$y = 0
sm_data[sm_data$trip3%in%pakatrips,"y"]=1

#Determine which species to include in the xij variables. 
#Which is 0 or 1 for each trip j, if species i is caught in it.
#Use those species i making up 99% of the cumulative catch
species_catch=aggregate(sm_data$lbs,by=list(sm_data$SPECIES),FUN=sum)
sco=species_catch[order(species_catch$x),]
sco$cum=cumsum(sco$x)
sco$perc=round(sco$cum/sum(sco$x),4)
x_species=sco[which(sco$perc>=0.01),"Group.1"] #Use these species (I checked: these all occur with paka)
x_species=x_species[-which(x_species==19)] #remove paka from the x variables

#Assign variable x<species code> to the dataset indicating which trips
#caught (1) or did not catch (0) the corresponding species
for(i in 1:length(x_species)){
  temp_trips=unique(sm_data[sm_data$SPECIES==x_species[i],"trip3"])
  temp_var=rep(0,nrow(sm_data))
  temp_var[which(sm_data$trip3%in%temp_trips)]=1
  sm_data$temp_var=temp_var
  names(sm_data)[ncol(sm_data)]=paste0("x",x_species[i])
}
#spec_list=read.csv("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\DAR_assessment_data\\raw\\fisher_species_5_6_16.csv",header=T)
spec_list=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\CPUE\\supporting files\\fisher_species_5_6_16.csv",header=T)
spec_names=spec_list[spec_list$SPECIES_PK%in%x_species,c("SPECIES_PK","SPECIES_NAME")]
spec_names=spec_names[match(x_species,spec_names$SPECIES_PK),]

##################################
#Calculate coefficiences and critical value
#################################
#FIRST reduce dataset to trip level
sm_data_trip=unique(sm_data[,c(1,29:67)])

#Generate coefficients from the model.
sm_model=glm(y~.,data=sm_data_trip[,-1],family=binomial(link='logit'))
#Check significance of predictors. 
#step(sm_model,direction="backward") #Removes 65,2,7,16,78 takes a long time
mcoef=data.frame("coef"=coef(sm_model)[-c(1,8,11,13,16,37)],"spec_names"=spec_names$SPECIES_NAME[-c(7,10,12,15,36)],"spec_pk"=spec_names$SPECIES_PK[-c(7,10,12,15,36)])
mcoef=mcoef[order(mcoef$coef),]
#Regression coefficient figure
par(mar=c(5,4,4,2))
a=barplot(mcoef$coef,horiz=T,xlab="Regression coefficients",yaxt="n",xlim=c(-3,2))
axis(2,at=a,labels=mcoef$spec_names,tick=F,las=1,cex.axis=0.9,hadj=0)

#Calculate critical value. 
prediction=predict.glm(sm_model,type="response")
sm_data_trip$y_predict=prediction
hist(sm_data_trip$y_predict) #good range of values
critval = function(data,val){
  obs=data$y
  include=as.numeric(data$y_predict>=val)
  obj=sum(abs(obs-include))
  return(c("obj"=obj,"trips_in"=sum(include)))
}
crits=sapply(seq(0,1,0.01),critval,data=sm_data_trip) #minimum is 0.51 [0.53 in update]
#Plots of critical value
par(mar=c(5,4,4,5))
plot(x=seq(0,1,0.01),y=crits["obj",],type="o",xlab="Probability",ylab="Trips: abs(obs-pred)")
plot(x=seq(0,1,0.01),y=crits["obj",],type="o",xlab="Probability",ylab="# of incorrect predictions")
abline(v=0.53)
par(new=T)
plot(x=seq(0,1,0.01),y=crits["trips_in",]/nrow(sm_data_trip),type="l",lty=2,xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext(side=4,line=3,"Percent of data included")
hist(sm_data_trip$y_predict,main="",xlab="Probability",ylab="Frequency")
abline(v=0.53,lwd=2)

#Output dataset of records from trips targetting paka based on Stephens and MacCall 2004
tdata7=merge(tdata6,sm_data_trip[,c("trip3","y","y_predict")],by="trip3")
pakadata=tdata7[tdata7$y_predict>=0.53,] #was 0.51
#write.csv(pakadata,"D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\UptoDate files\\Finalized_record_pakaCPUE_dataset.csv",row.names=F)
write.csv(pakadata,"C:\\Users\\John.Syslo\\Documents\\2020 Deep 7\\2020 Data prep files\\CPUE\\Finalized_record_pakaCPUE_dataset.csv",row.names=F)

# --OUTPUT: finalized trip dataset for CPUE complex and paka standardization---------------------
##COMPLEX
finalrecords=tdata6

##-##
#Ensure wind data values are only for the management grid used
#as the final area (area). Do this by removing non-final-area 
#records (from the multitrips) in the finalized record dataset.
##-##
for_wind=tdata6[tdata6$AREA==tdata6$area,]

##-##
#Reduce to finalized trip dataset
##-##
red_to_finaltrip=for_wind[,c("trip3","FYEAR","FISHED","d7catch","fisher","cum_exp","qtr",
                             "uku_lbs","houreffort","dayeffort","area","cpue","speed","xdir","ydir")]

#for 2021 data which has no wind variables:
red_to_finaltrip=tdata6[,c("FYEAR","FISHED","trip2","AREA","HOURS","days","lbs","deep7_lbs","fisher","cum_exp","qtr","uku_lbs")] 
                           
#the following does not reduce to single trips when applying unique function - SPECIES
red_to_finaltrip=tdata6[,c("FYEAR","FISHED","trip2","AREA","n_areas_e","HOURS","days","SPECIES","lbs","deep7_lbs","fisher","cum_exp","qtr","uku_lbs")] #, 
                          # "NO_GEARS_OR_NET_LENGTH", "NUMBER_LOST_TO_PREDATOR", "PREDATOR_CODE","PREDATOR")]



finaltrips=unique(red_to_finaltrip)
#write.csv(finaltrips,"D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\UptoDate files\\Finalized_tripCPUE_dataset_forStandardization.csv",row.names=F)
#write.csv(finaltrips,"C:\\Users\\John.Syslo\\Documents\\2020 Deep 7\\2020 Data prep files\\CPUE\\Finalized_tripCPUE_dataset_forStandardization_2017.csv",row.names=F) #bad added 5/20/20

write.csv(finaltrips,"C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\Raw FRS depredation\\Finalized_tripCPUE_dataset_forStandardization_2021_TWO.csv",row.names=F) #2021 data pull

#PAKA ONLY
finalrecords_paka=read.csv("C:\\Users\\John.Syslo\\Documents\\2020 Deep 7\\2020 Data prep files\\CPUE\\Finalized_record_pakaCPUE_dataset.csv",header=T)
finalrecords_paka$FISHED=as.Date(finalrecords_paka$FISHED)

#Calculate paka catch and cpue
#catch
pcatch_trip=aggregate(finalrecords_paka[finalrecords_paka$SPECIES==19,]$lbs,by=list(finalrecords_paka[finalrecords_paka$SPECIES==19,]$trip3),FUN=sum)
finalrecords_paka2=merge(finalrecords_paka,pcatch_trip,by.x="trip3",by.y="Group.1",all.x=T)
finalrecords_paka2[is.na(finalrecords_paka2$x),"x"]=0
names(finalrecords_paka2)[which(names(finalrecords_paka2)=="x")]="paka_catch"
#cpue
names(finalrecords_paka2)[which(names(finalrecords_paka2)=="cpue")]="d7cpue"
finalrecords_paka2$pcpue=NA
finalrecords_paka2[finalrecords_paka2$FISHED<"2002-10-01","pcpue"]=finalrecords_paka2[finalrecords_paka2$FISHED<"2002-10-01","paka_catch"]/finalrecords_paka2[finalrecords_paka2$FISHED<"2002-10-01","dayeffort"]
finalrecords_paka2[finalrecords_paka2$FISHED>="2002-10-01","pcpue"]=finalrecords_paka2[finalrecords_paka2$FISHED>="2002-10-01","paka_catch"]/finalrecords_paka2[finalrecords_paka2$FISHED>="2002-10-01","houreffort"]

#-Adjust cum_exp to reflect experience of paka trip
temp_paka=finalrecords_paka2[,c("trip3","fisher","FISHED")]
ord_uniq=unique(temp_paka)[order(unique(temp_paka)$fisher,unique(temp_paka)$trip3),] #unique fisher-trip combinations ordered by fisher
library(plyr)
num_exp=count(ord_uniq,vars="fisher") #number of unique fisher-trip combos by fisher
ord_uniq$cum_exp=unlist(lapply(num_exp$freq,FUN=seq)) #add experience value; from 1 to number of unique fisher-trip combinations per fisher 
temp_paka2=join(temp_paka,ord_uniq,by=c("trip3","FISHED","fisher")) #use join so that the same order is maintained with finalrecords_paka2
####----#####
finalrecords_paka2$pcum_exp=temp_paka2$cum_exp
####----#####

#Ensure wind data values are only for the management grid used
#as the final area (area). Do this by removing non-final-area 
#records (from the multitrips) in the finalized record dataset.
for_wind_paka=finalrecords_paka2[finalrecords_paka2$AREA==finalrecords_paka2$area,]

#Finalized trip dataset
red_to_finaltrip_paka=for_wind_paka[,c("trip3","FYEAR","FISHED","d7catch","paka_catch","fisher","cum_exp","pcum_exp","qtr",
                                       "uku_lbs","houreffort","dayeffort","area","pcpue","speed","xdir","ydir")]
finaltrips_paka=unique(red_to_finaltrip_paka)
#write.csv(finaltrips_paka,"D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\UptoDate files\\Finalized_trip_pakaCPUE_dataset.csv",row.names=F)
write.csv(finaltrips_paka,"C:\\Users\\John.Syslo\\Documents\\2020 Deep 7\\2020 Data prep files\\CPUE\\Finalized_trip_pakaCPUE_dataset.csv",row.names=F)

