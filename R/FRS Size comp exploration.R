#
rm(list=ls())

# --Combine input datasets into single dataset--------------------------
library(foreign)

#early=read.dbf("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\DAR_assessment_data\\PaulTao_data\\HFY48_93E.dbf",as.is=T)
early=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Raw data\\picdr_112849_fy48_15.csv")

#late=read.dbf("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\DAR_assessment_data\\PaulTao_data\\hf94c15F.dbf",as.is=T)
late=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Raw data\\picdr_112849_fy16_18.csv")

#add 2019 data here
d2019=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Raw data\\picdr_112970_fy19.csv")

#dim(early)
#dim(late)
#dim(d2019)
#names(early)
#names(late)
#names(d2019)

#not needed because fields match for update data:
#early[,c("O_LICENSE","O_LNAME","O_FNAME","O_VESSEL","O_USCG","O_HA_NO","CHNG")]=NA
fulldata=rbind(early,late,d2019)
#dim(fulldata) 

#Add lbs field as maximum of LBS and LBS_SOLD
fulldata$lbs=pmax(fulldata$LBS,fulldata$LBS_SOLD,na.rm=T)

################################################################
#2020- RELATED TO CPUE BUT EASY TO DO HERE:
#For report table 6, proportion of records with some name info. by year
XX<-aggregate(fulldata$FNAME,by=list(fulldata$FYEAR),FUN=length)


#dim(fulldata[which(fulldata$FNAME == "NULL" & fulldata$LNAME == "NULL"),])
no_name<-fulldata[which(fulldata$FNAME == "NULL" & fulldata$LNAME == "NULL"),]
y<-aggregate(no_name$FNAME,by=list(no_name$FYEAR),FUN=length)

perc_no_name<-1-y$x/XX$x[1:55]

#do not need because there is no 2019 data:
#Remove FY2016
####----#####
#reddata=fulldata[fulldata$FYEAR<=2015,]
####----#####

reddata=fulldata #included to reduce number of code changes

# --Use only records of Deep 7 species ----------------------------------
#Use common names of deep7 species since Gindai has two scientific names
com_names=c("HAPUUPUU","KALEKALE","OPAKAPAKA","EHU","ONAGA","ehu","LEHI","GINDAI")
##Species code 36 is also ehu according to Appendix II.4g in 
#Moffitt et al. 2011 "Bottomfish CPUE standardization workshop 
#proceedings August 4-6, 2008", but was discontinued in 1989. 
id=c(15,17,19,21,22,36,58,97)


####----#####
deep7=reddata[reddata$SPECIES%in%id,]
#dim(deep7)
####----#####

#Assign ehu species code 36 to ehu code 21.
deep7[deep7$SPECIES==36,"SPECIES"]=21

#calculate proprotion of Deep7 catch on gear3
a<-aggregate(deep7$lbs,by=list(deep7$GEAR),FUN=sum)
b<-aggregate(deep7$lbs,by=list(deep7$GEAR),FUN=length)

# --Use only records from the MHI------------------------------------
#Read in key table from last assessment for MHI areas (mhi_areas.csv)
#but keep only those grids Reg Kokbun identified as valid from the
#CPUE data workshop
#areasReggie=read.csv("D:\\Stock assessments\\Hawaii\\MHI_bottomfish_2017\\Data exploration\\Workshop_files_Aug16\\Tuesday analyses\\BF_Area_Grid_Reggie.csv",header=T)
areasReggie=read.csv("C:\\Users\\John.Syslo\\Documents\\2020_Deep_7\\2020 Data prep files\\Catch\\BF_Area_Grid_Reggie.csv",header=T)
valid=areasReggie[which(areasReggie$Valid.==""),]$area

#Remove non-valid subareas (A and B) from area 16123 as well as 
#records from 16123 without a subarea specified
####----#####
mhidata=deep7[deep7$AREA%in%valid,] #valid areas
mhidata=mhidata[!mhidata$SUBAREA%in%c("A","B"),] #remove known invalid subareas (16123A and 16123B)
mhidata=mhidata[!(mhidata$AREA==16123 & is.na(mhidata$SUBAREA)),] #remove the 16123 records without a subarea distinction
####----#####



# --Partition invalid areas to known mhi ares-----------------------------
#THE FOLLOWING STATS HAVE NOT BEEN UPDATED FOR 2020, but 1280 records are not valid
#Based on proportion of species/year in mhi area compared to nwhi ares
#Invalid areas include 0's, ponds, and non valid MHI areas
#5 areas (17 records) are ponds, 1249 records are area=0, 317 areas (274 records) aren't valid and aren't area=0
invalid=areasReggie[which(areasReggie$Valid.!=""),]$area
notvalid=deep7[deep7$AREA%in%invalid,]

#Need NWHI areas as well
#NWHI areas are those not in mhi_areas.csv
####----#####
nwhidata=deep7[!deep7$AREA%in%areasReggie$area,]
####----#####

#Calculate proportion of fish by species and year in MHI vs NWHI
agg_mhi=aggregate(mhidata$lbs,by=list("FYEAR"=mhidata$FYEAR,"SPECIES"=mhidata$SPECIES),FUN=sum)
mhi_table=matrix(agg_mhi$x,length(1948:2019),length(id[-6]),dimnames=list(1948:2019,id[-6]))

agg_nwhi=aggregate(nwhidata$lbs,by=list("FYEAR"=nwhidata$FYEAR,"SPECIES"=nwhidata$SPECIES),FUN=sum)
comb_agg=merge(agg_mhi,agg_nwhi,by=c("FYEAR","SPECIES"),all.x=T)
comb_agg[is.na(comb_agg$x.y),"x.y"]=0
comb_agg$percMHI=comb_agg$x.x/(comb_agg$x.x+comb_agg$x.y)


agg_inval=aggregate(notvalid$lbs,by=list("FYEAR"=notvalid$FYEAR,"SPECIES"=notvalid$SPECIES),FUN=sum)
comb_all=merge(comb_agg,agg_inval,by=c("FYEAR","SPECIES"),all=T)
comb_all$add_catch=comb_all$percMHI*comb_all$x



#################################7/6/21 -- New Stuff here####################################################

#column "CAUGHT" appears to contain the number caught per interview
#reduce to interviews with a single fish for each of the deep 7 species


(mhidata$SPECIES)
one_fish<-mhidata[mhidata$CAUGHT==1,]
dim(mhidata)
dim(one_fish) #reduced from 466493 records to 83592

 
valid<-function(x){sum(!is.na(x))} 
by_sp<-tapply(one_fish$LBS,one_fish$SPECIES,valid) #don't appear to be LBS of NA

#com_names=c("HAPUUPUU","KALEKALE","OPAKAPAKA","ONAGA","ehu","LEHI","GINDAI")
#id=c(15,17,19,21,22,58,97)

#number of records by sp for all years
#15824  7169 18102 10892 12965 11591  7049 

#################################
#look at paka #################
#####################################


paka_one<-one_fish[one_fish$SPECIES==19,]
mean(paka_one$LBS) #6.22318

hist(paka_one$LBS) 
range(paka_one$LBS) #0.5 to 278 lbs - need a maximum cutoff - 

#reduce extreme outliers and re-examine
paka_one<-one_fish[one_fish$SPECIES==19 & one_fish$LBS < 100,]
hist(paka_one$LBS) 
range(paka_one$LBS) #still one 98 pounder, histogram frequencies not visible after 40 lbs

#state records according to http://www.hawaiifishingnews.com/records.cfm

#Opakapaka 18.5 lbs
#Onaga 34.1875 lbs
#ehu 11.375 lbs
#Hapu 70 lbs
#Kalekale 3.5625 lbs
#Lehi 32.4375 lbs
#Gindai 4.29 lbs

paka_one<-one_fish[one_fish$SPECIES==19 & one_fish$LBS < 19.0,] #round up nearest lb to state record


#average weight
mean(paka_one$LBS) #5.798457 lbs
hist(paka_one$LBS) 

paka_wt_by_year<-tapply(paka_one$LBS, paka_one$FYEAR,mean)
plot(paka_wt_by_year) # fairly consistent

paka_n_by_year<-tapply(onaga_one$LBS, onaga_one$FYEAR,valid)
plot(onaga_n_by_year) #distinct peak - corresponds with high vatch and effort period


####################
#onaga##################
########################

onaga_one<-one_fish[one_fish$SPECIES==21 & one_fish$LBS < 35.0,] #round up nearest lb to state record

#average weight
mean(onaga_one$LBS) #2.97 lbs
hist(onaga_one$LBS)

onaga_wt_by_year<-tapply(onaga_one$LBS, onaga_one$FYEAR,mean)
plot(onaga_wt_by_year) # decreases through time 

onaga_n_by_year<-tapply(onaga_one$LBS, onaga_one$FYEAR,valid)
plot(onaga_n_by_year) #distinct peak - corresponds with high vatch and effort period

#######################
#ehu###################
######################

ehu_one<-one_fish[one_fish$SPECIES==22 & one_fish$LBS < 12.0,] #round up nearest lb to state record

#average weight
mean(ehu_one$LBS) #4.917945 lbs
hist(ehu_one$LBS)

ehu_wt_by_year<-tapply(ehu_one$LBS, ehu_one$FYEAR,mean)
plot(ehu_wt_by_year) # fairly consistent 

ehu_n_by_year<-tapply(ehu_one$LBS, ehu_one$FYEAR,valid)
plot(ehu_n_by_year) #distinct peak - corresponds with high vatch and effort period


#######################
#hapu###################
######################

hapu_one<-one_fish[one_fish$SPECIES==15 & one_fish$LBS < 70.0,] #round up nearest lb to state record

#average weight
mean(hapu_one$LBS) #9.8 lbs
hist(hapu_one$LBS)

hapu_wt_by_year<-tapply(hapu_one$LBS, hapu_one$FYEAR,mean)
plot(hapu_wt_by_year) # fairly consistent 

hapu_n_by_year<-tapply(hapu_one$LBS, hapu_one$FYEAR,valid)
plot(hapu_n_by_year) #distinct peak - corresponds with high vatch and effort period


#######################
#kalekale###################
######################

kale_one<-one_fish[one_fish$SPECIES==17 & one_fish$LBS < 4.0,] #round up nearest lb to state record

#average weight
mean(kale_one$LBS) #1.74 lbs
hist(kale_one$LBS)

kale_wt_by_year<-tapply(kale_one$LBS, kale_one$FYEAR,mean)
plot(kale_wt_by_year) # fairly consistent 

kale_n_by_year<-tapply(kale_one$LBS, kale_one$FYEAR,valid)
plot(kale_n_by_year) #distinct peak - corresponds with high vatch and effort period


#######################
#lehi###################
######################

lehi_one<-one_fish[one_fish$SPECIES==58 & one_fish$LBS < 33.0,] #round up nearest lb to state record

#average weight
mean(lehi_one$LBS) #8.87 lbs
hist(lehi_one$LBS)

lehi_wt_by_year<-tapply(lehi_one$LBS, lehi_one$FYEAR,mean)
plot(lehi_wt_by_year) # fairly consistent 

lehi_n_by_year<-tapply(lehi_one$LBS, lehi_one$FYEAR,valid)
plot(lehi_n_by_year) #distinct peak - corresponds with high vatch and effort period


#######################
#gindai###################
######################

gindai_one<-one_fish[one_fish$SPECIES==97 & one_fish$LBS < 5.0,] #round up nearest lb to state record

#average weight
mean(gindai_one$LBS) #2.06 lbs
hist(gindai_one$LBS)

gindai_wt_by_year<-tapply(gindai_one$LBS, gindai_one$FYEAR,mean)
plot(gindai_wt_by_year) # fairly consistent 

gindai_n_by_year<-tapply(gindai_one$LBS, gindai_one$FYEAR,valid)
plot(gindai_n_by_year) #distinct peak - corresponds with high vatch and effort period


