

# Nicholas Ducharme-Barth
# 04/04/2022
# Format BFISH camera data
# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Notes:
# There are two drops per PSU which is the sampling unit
# However, there are some PSUs that are not sampled on the same day
# Need to link drop information to abundance (Max N)
# Drop time in UTC needs to be recoded to HST

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
	data_flag = "2021_" # includes data through 2021
#_____________________________________________________________________________________________________________________________
# define helper function for converting DROP_TIME in 0-24 hours
	convert_drop_time = function(x)
	{	
		if(class(x)!="integer")
		{
			stop("Bad data type. Must be integer.")
		}
		
		# convert to 0-24 UTC
		if(!is.na(x))
		{
			if(nchar(x)==3)
			{
				hour = 0
				min = as.numeric(substr(x,1,1))/60
				sec = as.numeric(substr(x,2,3))/(60^2)
				time = hour + min + sec
			} else if(nchar(x)==4){
				hour = 0
				min = as.numeric(substr(x,1,2))/60
				sec = as.numeric(substr(x,3,4))/(60^2)
				time = hour + min + sec
			} else if(nchar(x)==5){
				hour = as.numeric(substr(x,1,1))
				min = as.numeric(substr(x,2,3))/60
				sec = as.numeric(substr(x,4,5))/(60^2)
				time = hour + min + sec
			} else if(nchar(x)==6){
				hour = as.numeric(substr(x,1,2))
				min = as.numeric(substr(x,3,4))/60
				sec = as.numeric(substr(x,5,6))/(60^2)
				time = hour + min + sec
			} else {
				time = NA
			}

			# convert to HST (subtract 10 hours)
				time = time - 10
				if(!is.na(time))
				{
					if(time<0)
					{
						time = time + 24
					}

					# assume data-entry error if HST earlier than 6am or later than 6pm
					# label UTC time as HST time in these cases
					if(time<6|time>18)
					{
						time = time + 10
						if(time>24)
						{
							time = time - 24
						}
					}
				}
		} else {
			time = NA
		} 
		return(time)
	}


#_____________________________________________________________________________________________________________________________
# bring in camera data
	# bring in conversion factor data
	species_dt = fread(paste0(proj.dir,"Data/SPECIES_TABLE.csv")) %>%
				 .[SPECIES_CD %in% c("ETCO","ETCA","PRSI","PRFI","PRZO","HYQU","APRU")] %>%
				 unique(.)

	# PSU specific information
	PSU_table = fread(paste0(proj.dir,"Data/BFISH PSU lookup table.csv")) %>%
				.[,.(PSU,Island,lon_deg,lat_deg,STRATA,STRATA_2020,Depth_MEDIAN_m,med_slp,med_acr,BS_pct_over_136j,pctHB,pctHS)] %>%
				.[,substrate:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][1])] %>%
				.[,slope:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][2])]

	# catches (one row per length measurement)
	# remove these drops/PSUs
	dark_drops = fread(paste0(proj.dir,"Data/",data_flag,"CAM_MAXN.csv"))[SPECIES_CD == "DARK"]$DROP_CD
	BFISH_CAM_COUNT = fread(paste0(proj.dir,"Data/",data_flag,"CAM_MAXN.csv")) %>%
			  .[,.(MAXN=sum(MAXN)),by=.(DROP_CD,SPECIES_CD)] %>%
			  .[SPECIES_CD %in% c("ETCO","ETCA","PRSI","PRFI","PRZO","HYQU","APRU")] %>%
			  .[!(DROP_CD %in% dark_drops)]
	BFISH_CAM_LENGTHS = fread(paste0(proj.dir,"Data/",data_flag,"CAM_LENGTHS.csv")) %>%
			      .[,LENGTH_CM:=round(MEAN_MM/10)] %>%
				  .[,.(DROP_CD,SPECIES_CD,LENGTH_CM)] %>%
				  .[,.N,by=.(DROP_CD,SPECIES_CD,LENGTH_CM)] %>%
				  .[,TOTAL_N:=sum(N),by=.(DROP_CD,SPECIES_CD)] %>%
				  .[,PROP_N:=N/TOTAL_N] %>%
				  .[,.(DROP_CD,SPECIES_CD,LENGTH_CM,PROP_N)] %>%
				  .[!(DROP_CD %in% dark_drops)]

	BFISH_CAM_C = merge(BFISH_CAM_LENGTHS,BFISH_CAM_COUNT,by=c("DROP_CD","SPECIES_CD"),all.x=TRUE,all.y=TRUE) %>%
			  	  merge(.,species_dt[,.(SPECIES_CD,A,B,GCF)],by="SPECIES_CD") %>%
			  .[,N:=MAXN*PROP_N] %>%
			  .[,STD_N:=N] %>% 
			  .[,Biomass:=A * (LENGTH_CM^B)] %>%
			  .[,KG:=Biomass*N] %>%
			  .[,STD_KG:=Biomass*STD_N] %>%
			  .[,Length_category:=as.character(NA)] %>%
			  .[LENGTH_CM>=29,Length_category:="Exploitable"] %>%
			  .[LENGTH_CM<29,Length_category:="Un-exploitable"] %>%
			  .[is.na(LENGTH_CM),Length_category:="Unknown"]
	# remove drops where length is missing
	missing_lengths = unique(BFISH_CAM_C[Length_category == "Unknown"]$DROP_CD)

	BFISH_CAM_C = BFISH_CAM_C %>%
			  .[,.(N=sum(N,na.rm=TRUE),STD_N=sum(N,na.rm=TRUE),KG=sum(KG,na.rm=TRUE),STD_KG=sum(STD_KG,na.rm=TRUE)),by=.(DROP_CD,SPECIES_CD,Length_category)] %>%
			  # keep only exploitable "sized" biomass
			  .[Length_category=="Exploitable"] %>%
			  .[,.(DROP_CD,SPECIES_CD,N,KG,STD_N,STD_KG)]

	BFISH_CAM_C_long = copy(BFISH_CAM_C)
	BFISH_CAM_C = BFISH_CAM_C %>% .[!(DROP_CD %in% c(missing_lengths))] %>%
			  .[,.(DROP_CD,SPECIES_CD,KG)] %>%
			  # keep only biomass measurement
			  dcast(.,DROP_CD~SPECIES_CD,value.var="KG",fill=0,fun.aggregate=sum) %>%
			  .[,.(DROP_CD,ETCO,ETCA,PRSI,PRFI,PRZO,HYQU,APRU)]
			

	# drop-specific information
	BFISH_CAM_S = fread(paste0(proj.dir,"Data/",data_flag,"CAM_SAMPLE_TIME.csv")) %>%
			  .[,.(DROP_CD,DROP_DATE,DROP_TIME,VESSEL,GEAR_ID,PSU,OBS_LON,OBS_LAT,OFFICIAL_DEPTH_M,OFFICIAL_TEMP_C)] %>%
			  .[,DROP_TIME_HST:=sapply(DROP_TIME,convert_drop_time)] %>%
			  .[,SAMPLE_DATE:=as.POSIXct(as.character(DROP_DATE),format=c("%Y%m%d"))] %>%
			  .[,YEAR:=format(SAMPLE_DATE,format="%Y")] %>%
			  .[,MONTH:=format(SAMPLE_DATE,format="%m")] %>%
			  .[,DAY:=format(SAMPLE_DATE,format="%d")] %>%
			  .[,JD:=as.numeric(format(SAMPLE_DATE,format="%j"))] %>%
			  .[,YEAR_continuous:=as.numeric(YEAR)+(JD-1)/366] %>%
			  .[,LUNAR_PHASE:=getMoonIllumination(format(SAMPLE_DATE, format="%Y-%m-%d"))$fraction] %>%
			  .[,.(DROP_CD,SAMPLE_DATE,YEAR,MONTH,DAY,JD,YEAR_continuous,LUNAR_PHASE,DROP_TIME_HST,VESSEL,PSU,OBS_LON,OBS_LAT,OFFICIAL_DEPTH_M,OFFICIAL_TEMP_C)]  %>%
			  .[!(DROP_CD %in% dark_drops)] %>%
			  .[!(DROP_CD %in% c(missing_lengths))]

	camera_dt = merge(BFISH_CAM_S,BFISH_CAM_C,by=c("DROP_CD"),all=TRUE) %>%
						  # this next line drops samples with bad PSUs
						  merge(.,PSU_table[,.(PSU,Island,STRATA,STRATA_2020,Depth_MEDIAN_m,substrate,slope,med_slp,med_acr,BS_pct_over_136j,pctHB,pctHS)],by="PSU") %>%
						  .[,depth_strata:=ifelse(Depth_MEDIAN_m<75,"Z",ifelse(Depth_MEDIAN_m<200,"S",ifelse(Depth_MEDIAN_m<300,"M",ifelse(Depth_MEDIAN_m<400,"D","ZZ"))))] %>%
						  .[,depth_strata_2020:=ifelse(Depth_MEDIAN_m<75,"Z",ifelse(Depth_MEDIAN_m<110,"D1",ifelse(Depth_MEDIAN_m<170,"D2",ifelse(Depth_MEDIAN_m<200,"D3",ifelse(Depth_MEDIAN_m<330,"D4",ifelse(Depth_MEDIAN_m<400,"D5","ZZ"))))))] %>%
						  .[,complexity:=ifelse(med_acr<4,"MA1",ifelse(med_acr<9,"MA2","MA3"))] %>%
						  .[,hardness:=ifelse(BS_pct_over_136j<0.24,"HB1",ifelse(BS_pct_over_136j<0.46,"HB2","HB3"))] %>%
						  .[!is.na(PSU)] %>%
						  .[!(depth_strata%in%c("Z","ZZ"))] %>%
						  .[is.na(APRU),APRU:=0] %>%
						  .[is.na(ETCA),ETCA:=0] %>%
						  .[is.na(ETCO),ETCO:=0] %>%
						  .[is.na(HYQU),HYQU:=0] %>%
						  .[is.na(PRFI),PRFI:=0] %>%
						  .[is.na(PRSI),PRSI:=0] %>%
						  .[is.na(PRZO),PRZO:=0] %>%
						  # drop samples with missing values
						  na.omit(.)

	# save formatted data
		save(BFISH_CAM_S,file=paste0(proj.dir,"Data/",data_flag,"BFISH_CAM_S.RData"))
		save(BFISH_CAM_C,file=paste0(proj.dir,"Data/",data_flag,"BFISH_CAM_C.RData"))
		save(camera_dt,file=paste0(proj.dir,"Data/",data_flag,"camera_dt.RData"))
