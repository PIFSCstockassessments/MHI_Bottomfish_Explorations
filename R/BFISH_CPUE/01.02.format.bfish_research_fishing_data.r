

# Nicholas Ducharme-Barth
# 2023/05/31
# Format BFISH research fishing data
# Copyright (c) 2023 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(suncalc)


#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"

#_____________________________________________________________________________________________________________________________
# define data_flag
	# data_flag = "" # only loads data up through 2020
	data_flag = "2022_" # includes data through 2021

#_____________________________________________________________________________________________________________________________
# define helper function for converting DRIFT_START_TIME & DRIFT_END_TIME to 0-24 decimal
	convert_drift_time = function(x)
	# last two digits give minutes
	# first digit if 3 digits long is hour
	# first two digits if 4 digits long is hour
	{	
		if(class(x)!="integer")
		{
			stop("Bad data type. Must be 3 or 4 digit integer.")
		}
		if(is.na(x)){return(NA)}
		# fishing does not take place over night so assume data entry error and shift to afternoon
		# e.g. 145 should probably be 1345
		if(x<500)
		{
			x = x + 1200
		}

		chr_n = nchar(x)
		if(chr_n==3)
		{
			hour = as.numeric(substr(x, 1, 1))
			min = as.numeric(substr(x, 2, 3))
			if(min>59){stop("Bad data input. Minutes can not be greater than 59.")}
			time = hour + min/60
		} else {
			hour = as.numeric(substr(x, 1, 2))
			min = as.numeric(substr(x, 3, 4))
			if(hour>23){stop("Bad data input. Hours can not be greater than 23.")}
			if(min>59){stop("Bad data input. Minutes can not be greater than 59.")}
			time = hour + min/60
		}
		return(time)
	}

#_____________________________________________________________________________________________________________________________
# bring in research fishing data
	# drift-specific information
	BFISH_D = fread(paste0(proj.dir,"Data/",data_flag,"CRF_DRIFT.csv")) %>%
			  .[,.(BFISH,SAMPLE_ID,DRIFT_START_TIME,START_DEPTH_M,START_LON,START_LAT,DRIFT_END_TIME,END_DEPTH_M,END_LON,END_LAT)] %>%
			  .[,DRIFT_START_TIME:=sapply(DRIFT_START_TIME,convert_drift_time)] %>%
			  .[,DRIFT_END_TIME:=sapply(DRIFT_END_TIME,convert_drift_time)] %>%
			  .[,.(TIME_MIN=min(DRIFT_START_TIME,na.rm=TRUE),START_DEPTH_M=mean(START_DEPTH_M,na.rm=TRUE),START_LON=mean(START_LON,na.rm=TRUE),START_LAT=mean(START_LAT,na.rm=TRUE),TIME_MAX=max(DRIFT_END_TIME,na.rm=TRUE),END_DEPTH_M=mean(END_DEPTH_M,na.rm=TRUE),END_LON=mean(END_LON,na.rm=TRUE),END_LAT=mean(END_LAT,na.rm=TRUE)),by=.(BFISH,SAMPLE_ID)] %>%
			  .[,.(TIME_MIN,TIME_MAX,TIME_MEAN=mean(c(TIME_MIN,TIME_MAX),na.rm=TRUE),DEPTH_M = mean(c(START_DEPTH_M,END_DEPTH_M),na.rm=TRUE),LON=mean(c(START_LON,END_LON),na.rm=TRUE),LAT=mean(c(START_LAT,END_LAT),na.rm=TRUE)),by=.(BFISH,SAMPLE_ID)]

	# sample-specific (i.e., PSU) information
	BFISH_S = fread(paste0(proj.dir,"Data/",data_flag,"CRF_SAMPLE.csv")) %>%
			  .[,.(BFISH,SAMPLE_ID,PSU,SAMPLE_DATE,VESSEL,CAPTAIN_CD,OBSERVER,WIND_SPEED_KT,WAVE_HEIGHT_FT,CURRENT_CD)] %>%
			  .[,SAMPLE_DATE:=as.POSIXct(SAMPLE_DATE,tz="HST",format=c("%Y-%m-%d"))] %>%
			  .[,YEAR:=format(SAMPLE_DATE,format="%Y")] %>%
			  .[,MONTH:=format(SAMPLE_DATE,format="%m")] %>%
			  .[,DAY:=format(SAMPLE_DATE,format="%d")] %>%
			  .[,JD:=as.numeric(format(SAMPLE_DATE,format="%j"))] %>%
			  .[,YEAR_continuous:=as.numeric(YEAR)+(JD-1)/366] %>%
			  .[,LUNAR_PHASE:=getMoonIllumination(format(SAMPLE_DATE, format="%Y-%m-%d"))$fraction]

	# PSU specific information
	PSU_table = fread(paste0(proj.dir,"Data/BFISH PSU lookup table.csv")) %>%
				.[,.(PSU,Island,lon_deg,lat_deg,STRATA,STRATA_2020,Depth_MEDIAN_m,med_slp,med_acr,BS_pct_over_136j,pctHB,pctHS)] %>%
				.[,substrate:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][1])] %>%
				.[,slope:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][2])]

	# recalculate PSU based on actual location of drop
	tmp_dt = merge(BFISH_S[,.(BFISH,SAMPLE_ID,PSU)],BFISH_D[,.(BFISH,SAMPLE_ID,LON,LAT)],by=c("BFISH","SAMPLE_ID"))
	close_psu_vec = rep(NA,nrow(tmp_dt))
	for(i in seq_along(close_psu_vec))
	{
		tmp_lon = tmp_dt$LON[i]
		tmp_lat = tmp_dt$LAT[i]

		tmp_dist = geosphere::distHaversine(c(tmp_lon,tmp_lat),as.matrix(PSU_table[,.(lon_deg,lat_deg)]))
		close_psu_vec[i] = PSU_table$PSU[which(tmp_dist==min(tmp_dist))]
		rm(list=c("tmp_lon","tmp_lat","tmp_dist"))
	}
	mean(tmp_dt$PSU == close_psu_vec)
	tmp_dt$PSU = close_psu_vec
	BFISH_S = BFISH_S %>%
			  .[,PSU:=NULL] %>%
			  merge(.,tmp_dt[,.(BFISH,SAMPLE_ID,PSU)],by=c("BFISH","SAMPLE_ID"))

	# bring in conversion factor data
	species_dt = fread(paste0(proj.dir,"Data/SPECIES_TABLE.csv")) %>%
				 .[SPECIES_CD %in% c("ETCO","ETCA","PRSI","PRFI","PRZO","HYQU","APRU")] %>%
				 unique(.)


	# catches (one row per length measurement)
	BFISH_C = fread(paste0(proj.dir,"Data/",data_flag,"CRF_CATCH.csv")) %>%
			  .[,.N,by=.(BFISH,SAMPLE_ID,BAIT_CD,SPECIES_CD,LENGTH_CM)] %>%
			  .[SPECIES_CD %in% c("ETCO","ETCA","PRSI","PRFI","PRZO","HYQU","APRU")] %>%
			  merge(.,species_dt[,.(SPECIES_CD,A,B,GCF)]) %>%
			  .[,STD_N:=N/GCF] %>% 
			  .[,LENGTH_CM:=round(as.numeric(LENGTH_CM))] %>%
			  .[,Biomass:=A * (LENGTH_CM^B)] %>%
			  .[,KG:=Biomass*N] %>%
			  .[,STD_KG:=Biomass*STD_N] %>%
			  .[,Length_category:=as.character(NA)] %>%
			  .[LENGTH_CM>=29,Length_category:="Exploitable"] %>%
			  .[LENGTH_CM<29,Length_category:="Un-exploitable"] %>%
			  .[is.na(LENGTH_CM),Length_category:="Unknown"] %>%
			  .[,.(N=sum(N),STD_N=sum(STD_N),KG=sum(KG),STD_KG=sum(STD_KG)),by=.(BFISH,SAMPLE_ID,BAIT_CD,SPECIES_CD,Length_category)]
	missing_lengths = unique(BFISH_C[Length_category=="Unknown"]$SAMPLE_ID)
	BFISH_C = BFISH_C %>%
			  # keep only exploitable "sized" biomass
			  .[Length_category=="Exploitable"] %>%
			  .[,.(BFISH,SAMPLE_ID,BAIT_CD,SPECIES_CD,N,KG,STD_N,STD_KG)]
			  # .[SPECIES_CD=="SQSP",SPECIES_CD:="SQMI"] # %>%
			  # .[,SPECIES_GRP:="Deep7"] %>%
			  # .[SPECIES_CD %in% c("APVI","SEDU","SQMI","SQSP"), SPECIES_GRP:="Other"]
	BFISH_C = BFISH_C %>% 
			  .[,.(BFISH,SAMPLE_ID,BAIT_CD,SPECIES_CD,KG)] %>%
			  # keep only biomass measurement
			  dcast(.,BFISH+SAMPLE_ID+BAIT_CD~SPECIES_CD,value.var="KG",fill=0,fun.aggregate=sum) %>%
			  .[,.(BFISH,SAMPLE_ID,BAIT_CD,ETCO,ETCA,PRSI,PRFI,PRZO,HYQU,APRU)] %>%
			  # collapse across bait categories
			  .[,.(ETCO=sum(ETCO),ETCA=sum(ETCA),PRSI=sum(PRSI),PRFI=sum(PRFI),PRZO=sum(PRZO),HYQU=sum(HYQU),APRU=sum(APRU)),by=.(BFISH,SAMPLE_ID)]

	research_fishing_dt = merge(BFISH_S,BFISH_D,by=c("BFISH","SAMPLE_ID")) %>%
						  # this next line drops samples with bad PSUs
						  merge(.,PSU_table[,.(PSU,Island,STRATA,STRATA_2020,Depth_MEDIAN_m,substrate,slope,med_slp,med_acr,BS_pct_over_136j,pctHB,pctHS)],by="PSU") %>%
						  .[,depth_strata:=ifelse(Depth_MEDIAN_m<75,"Z",ifelse(Depth_MEDIAN_m<200,"S",ifelse(Depth_MEDIAN_m<300,"M",ifelse(Depth_MEDIAN_m<400,"D","ZZ"))))] %>%
						  .[,depth_strata_2020:=ifelse(Depth_MEDIAN_m<75,"Z",ifelse(Depth_MEDIAN_m<110,"D1",ifelse(Depth_MEDIAN_m<170,"D2",ifelse(Depth_MEDIAN_m<200,"D3",ifelse(Depth_MEDIAN_m<330,"D4",ifelse(Depth_MEDIAN_m<400,"D5","ZZ"))))))] %>%
						  .[,complexity:=ifelse(med_acr<4,"MA1",ifelse(med_acr<9,"MA2","MA3"))] %>%
						  .[,hardness:=ifelse(BS_pct_over_136j<0.24,"HB1",ifelse(BS_pct_over_136j<0.46,"HB2","HB3"))] %>%
						  # .[,STRATA:=paste0(substrate,"_",slope,"_",depth_strata)] %>%
						  # .[,STRATA_2020:=as.character(NA)] %>%
						  # .[depth_strata_2020=="D1"&complexity%in%c("MA1","MA2")&hardness%in%c("HB1","HB2","HB3"),STRATA_2020:="S01"] %>%
						  # .[depth_strata_2020=="D1"&complexity%in%c("MA3")&hardness%in%c("HB1"),STRATA_2020:="S02"] %>%
						  # .[depth_strata_2020=="D1"&complexity%in%c("MA3")&hardness%in%c("HB2"),STRATA_2020:="S03"] %>%
						  # .[depth_strata_2020=="D1"&complexity%in%c("MA3")&hardness%in%c("HB3"),STRATA_2020:="S04"] %>%
						  # .[depth_strata_2020=="D2"&complexity%in%c("MA1")&hardness%in%c("HB1","HB2"),STRATA_2020:="S05"] %>%
						  # .[depth_strata_2020=="D2"&complexity%in%c("MA1")&hardness%in%c("HB3"),STRATA_2020:="S06"] %>%
						  # .[depth_strata_2020=="D2"&complexity%in%c("MA2")&hardness%in%c("HB1"),STRATA_2020:="S07"] %>%
						  # .[depth_strata_2020=="D2"&complexity%in%c("MA2")&hardness%in%c("HB2"),STRATA_2020:="S08"] %>%
						  # .[depth_strata_2020=="D2"&complexity%in%c("MA2")&hardness%in%c("HB3"),STRATA_2020:="S09"] %>%
						  # .[depth_strata_2020=="D2"&complexity%in%c("MA3")&hardness%in%c("HB1"),STRATA_2020:="S10"] %>%
						  # .[depth_strata_2020=="D2"&complexity%in%c("MA3")&hardness%in%c("HB2"),STRATA_2020:="S11"] %>%
						  # .[depth_strata_2020=="D2"&complexity%in%c("MA3")&hardness%in%c("HB3"),STRATA_2020:="S12"] %>%
						  # .[depth_strata_2020=="D3"&complexity%in%c("MA1","MA2")&hardness%in%c("HB1","HB2","HB3"),STRATA_2020:="S13"] %>%
						  # .[depth_strata_2020=="D3"&complexity%in%c("MA3")&hardness%in%c("HB1"),STRATA_2020:="S14"] %>%
						  # .[depth_strata_2020=="D3"&complexity%in%c("MA3")&hardness%in%c("HB2"),STRATA_2020:="S15"] %>%
						  # .[depth_strata_2020=="D3"&complexity%in%c("MA3")&hardness%in%c("HB3"),STRATA_2020:="S16"] %>%		
						  # .[depth_strata_2020=="D4"&complexity%in%c("MA1","MA2")&hardness%in%c("HB1","HB2"),STRATA_2020:="S17"] %>%
						  # .[depth_strata_2020=="D4"&complexity%in%c("MA1")&hardness%in%c("HB3"),STRATA_2020:="S18"] %>%
						  # .[depth_strata_2020=="D4"&complexity%in%c("MA2")&hardness%in%c("HB3"),STRATA_2020:="S19"] %>%
						  # .[depth_strata_2020=="D4"&complexity%in%c("MA3")&hardness%in%c("HB1"),STRATA_2020:="S20"] %>%
						  # .[depth_strata_2020=="D4"&complexity%in%c("MA3")&hardness%in%c("HB2"),STRATA_2020:="S21"] %>%
						  # .[depth_strata_2020=="D4"&complexity%in%c("MA3")&hardness%in%c("HB3"),STRATA_2020:="S22"] %>%
						  # .[depth_strata_2020=="D5"&complexity%in%c("MA1","MA2")&hardness%in%c("HB1","HB2","HB3"),STRATA_2020:="S23"] %>%
						  # .[depth_strata_2020=="D5"&complexity%in%c("MA3")&hardness%in%c("HB1","HB2","HB3"),STRATA_2020:="S24"] %>%
						  merge(.,BFISH_C,by=c("BFISH","SAMPLE_ID"),all=TRUE) %>%
						  # this next line drops samples with missing PSUs
						  .[!is.na(PSU)] %>%
						  # drops PSUs outside of EFH 75-400m
						  .[!(depth_strata%in%c("Z","ZZ"))] %>%
						  .[is.na(APRU),APRU:=0] %>%
						  .[is.na(ETCA),ETCA:=0] %>%
						  .[is.na(ETCO),ETCO:=0] %>%
						  .[is.na(HYQU),HYQU:=0] %>%
						  .[is.na(PRFI),PRFI:=0] %>%
						  .[is.na(PRSI),PRSI:=0] %>%
						  .[is.na(PRZO),PRZO:=0] %>%
						  # drop samples with missing values
						  na.omit(.) %>%
			   			  # exclude samples where lengths are missing
			  			  .[!(SAMPLE_ID %in% c(missing_lengths))] %>%
						  # define sampling_unit variable
						  .[,SEASON:=ifelse(as.numeric(MONTH)>6,"Fall","Spring")] %>%
						  .[,design_sampling_unit:=SAMPLE_ID] %>%
						  .[,model_sampling_unit:=SAMPLE_ID]

						
	
	# save formatted data
		save(BFISH_D,file=paste0(proj.dir,"Data/",data_flag,"BFISH_D.RData"))
		save(BFISH_S,file=paste0(proj.dir,"Data/",data_flag,"BFISH_S.RData"))
		save(BFISH_C,file=paste0(proj.dir,"Data/",data_flag,"BFISH_C.RData"))
		save(research_fishing_dt,file=paste0(proj.dir,"Data/",data_flag,"research_fishing_dt.RData"))
