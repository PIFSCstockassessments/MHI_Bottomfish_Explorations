

# Nicholas Ducharme-Barth
# 03/09/2022
# Format BFISH research fishing data
# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)

#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"
	setwd(proj.dir)
	source("D:/HOME/SAP/Code/Utilities/turbo.r")

	plot_dir = paste0(proj.dir,"Plot/BFISH_CPUE/")
	dir.create(plot_dir,recursive=TRUE,showWarnings=FALSE)

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
	BFISH_D = fread(paste0(proj.dir,"Data/CRF_DRIFT.csv")) %>%
			  .[,.(BFISH,SAMPLE_ID,DRIFT_START_TIME,START_DEPTH_M,START_LON,START_LAT,DRIFT_END_TIME,END_DEPTH_M,END_LON,END_LAT)] %>%
			  .[,DRIFT_START_TIME:=sapply(DRIFT_START_TIME,convert_drift_time)] %>%
			  .[,DRIFT_END_TIME:=sapply(DRIFT_END_TIME,convert_drift_time)] %>%
			  .[,.(TIME_MIN=min(DRIFT_START_TIME,na.rm=TRUE),START_DEPTH_M=mean(START_DEPTH_M,na.rm=TRUE),START_LON=mean(START_LON,na.rm=TRUE),START_LAT=mean(START_LAT,na.rm=TRUE),TIME_MAX=max(DRIFT_END_TIME,na.rm=TRUE),END_DEPTH_M=mean(END_DEPTH_M,na.rm=TRUE),END_LON=mean(END_LON,na.rm=TRUE),END_LAT=mean(END_LAT,na.rm=TRUE)),by=.(BFISH,SAMPLE_ID)] %>%
			  .[,.(TIME_MIN,TIME_MAX,TIME_MEAN=mean(c(TIME_MIN,TIME_MAX),na.rm=TRUE),DEPTH_M = mean(c(START_DEPTH_M,END_DEPTH_M),na.rm=TRUE),LON=mean(c(START_LON,END_LON),na.rm=TRUE),LAT=mean(c(START_LAT,END_LAT),na.rm=TRUE)),by=.(BFISH,SAMPLE_ID)]

	# sample-specific (i.e., PSU) information
	BFISH_S = fread(paste0(proj.dir,"Data/CRF_SAMPLE.csv")) %>%
			  .[,.(BFISH,SAMPLE_ID,PSU,SAMPLE_DATE,VESSEL,CAPTAIN_CD,OBSERVER,WIND_SPEED_KT,WAVE_HEIGHT_FT,CURRENT_CD)] %>%
			  .[,SAMPLE_DATE:=as.POSIXct(SAMPLE_DATE,tz="HST",format=c("%Y-%m-%d"))] %>%
			  .[,YEAR:=format(SAMPLE_DATE,format="%Y")] %>%
			  .[,MONTH:=format(SAMPLE_DATE,format="%m")] %>%
			  .[,DAY:=format(SAMPLE_DATE,format="%m")] %>%
			  .[,JD:=format(SAMPLE_DATE,format="%j")] %>%
			  .[,YEAR_continuous:=as.numeric(YEAR)+(as.numeric(JD)-1)/366]
	# catches (one row per individual)
	BFISH_C = fread(paste0(proj.dir,"Data/CRF_CATCH.csv")) %>%
			  .[,.N,by=.(BFISH,SAMPLE_ID,BAIT_CD,SPECIES_CD)] %>%
			  .[SPECIES_CD %in% c("APVI","SEDU","SQMI","SQSP","ETCO","ETCA","PRSI","PRFI","PRZO","HYQU","APRU")] %>%
			  .[SPECIES_CD=="SQSP",SPECIES_CD:="SQMI"] # %>%
			  # .[,SPECIES_GRP:="Deep7"] %>%
			  # .[SPECIES_CD %in% c("APVI","SEDU","SQMI","SQSP"), SPECIES_GRP:="Other"]
	samples_missing_bait = unique(BFISH_C[is.na(BAIT_CD)]$SAMPLE_ID)
	# exclude 15 samples where bait was missing for recorded fish
	BFISH_C = BFISH_C %>% .[!(SAMPLE_ID %in% c(samples_missing_bait))] %>%
			  dcast(.,BFISH+SAMPLE_ID+BAIT_CD~SPECIES_CD,value.var="N",fill=0,fun.aggregate=sum) %>%
			  .[,.(BFISH,SAMPLE_ID,BAIT_CD,APVI,SEDU,SQMI,ETCO,ETCA,PRSI,PRFI,PRZO,HYQU,APRU)]
	
	# PSU specific information
	PSU_table = fread(paste0(proj.dir,"Data/BFISH PSU lookup table.csv")) %>%
				.[,.(PSU,Island,lon_deg,lat_deg,STRATA,STRATA_2020,Depth_MEDIAN_m,Depth_MEAN_m,Depth_MAJORITY_m,Depth_MIN_m,Depth_MAX_m,Depth_STD_m,med_slp,med_acr,Depth_px_n,Slp_px_n,BS_px_nj,BS_px_over_136j,BS_pct_over_136j,pctHB,pctHS)]

	research_fishing_dt = merge(BFISH_S,BFISH_D,by=c("BFISH","SAMPLE_ID")) %>%
						  # this next line drops 17 samples with bad PSUs
						  merge(.,PSU_table[,.(PSU,Island,STRATA,STRATA_2020)],by="PSU") %>%
						  .[,BAIT_CD:="S"]
	tmp_dt = copy(research_fishing_dt)
	tmp_dt$BAIT_CD = "F"
	research_fishing_dt = rbind(research_fishing_dt,tmp_dt) %>%
						  merge(.,BFISH_C,by=c("BFISH","SAMPLE_ID","BAIT_CD"),all=TRUE) %>%
						  .[is.na(APRU),APRU:=0] %>%
						  .[is.na(APVI),APVI:=0] %>%
						  .[is.na(ETCA),ETCA:=0] %>%
						  .[is.na(ETCO),ETCO:=0] %>%
						  .[is.na(HYQU),HYQU:=0] %>%
						  .[is.na(PRFI),PRFI:=0] %>%
						  .[is.na(PRSI),PRSI:=0] %>%
						  .[is.na(PRZO),PRZO:=0] %>%
						  .[is.na(SEDU),SEDU:=0] %>%
						  .[is.na(SQMI),SQMI:=0]
	
	# save formatted data
		save(BFISH_D,file=paste0(proj.dir,"Data/BFISH_D.RData"))
		save(BFISH_S,file=paste0(proj.dir,"Data/BFISH_S.RData"))
		save(BFISH_C,file=paste0(proj.dir,"Data/BFISH_C.RData"))
		save(PSU_table,file=paste0(proj.dir,"Data/PSU_table.RData"))
		save(research_fishing_dt,file=paste0(proj.dir,"Data/research_fishing_dt.RData"))

