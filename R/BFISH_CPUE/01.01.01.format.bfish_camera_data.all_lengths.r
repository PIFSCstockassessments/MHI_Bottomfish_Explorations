

# Nicholas Ducharme-Barth
# 2023/05/22
# Format BFISH camera data: Option 1
# The sampling_unit is defined as each individual drop
# Drops are excluded if they are "dark" or if they have missing length measurements
# Make sure to convert camera drop dates and times from UTC to HST
# Copyright (c) 2023 Nicholas Ducharme-Barth
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
	# data_flag = "2022_" # includes data through 2022
#_____________________________________________________________________________________________________________________________
# define helper function for
	convert_date_time = function(input_date,input_time,tz_input="UTC",tz_output="HST")
	{
		hours = trunc(input_time)
		minutes_raw = (input_time-hours)*60
		minutes = trunc(minutes_raw)
		seconds = trunc((minutes_raw-minutes)*60)
		if(nchar(hours)<2)
		{
			hours = paste0("0",hours)
		} else {
			hours = as.character(hours)
		}
		if(nchar(minutes)<2)
		{
			minutes = paste0("0",minutes)
		} else {
			minutes = as.character(minutes)
		}
		if(nchar(seconds)<2)
		{
			seconds = paste0("0",seconds)
		} else {
			seconds = as.character(seconds)
		}
		time_chr = paste0(hours,":",minutes,":",seconds)
		utc_date_time = strptime(paste0(input_date," ",time_chr),format="%Y%m%d %H:%M:%S",tz=tz_input)
		new_time = utc_date_time = as.POSIXct(utc_date_time, tz=tz_input)
		attr(new_time,"tzone")=tz_output
		output_date = as.POSIXct(strptime(gsub("-","",strsplit(as.character(new_time)," ")[[1]][1]),format="%Y%m%d",tz=tz_output), tz=tz_output)
		output_time_raw = as.numeric(strsplit(strsplit(as.character(new_time)," ")[[1]][2],":")[[1]])
		output_time_raw[2] = output_time_raw[2]/60
		output_time_raw[3] = output_time_raw[3]/60^2
		output_time = sum(output_time_raw)
		return(list(date=output_date,time=output_time))
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
			  .[SPECIES_CD %in% c("ETCO","ETCA","PRSI","PRFI","PRZO","HYQU","APRU")]
	BFISH_CAM_LENGTHS = fread(paste0(proj.dir,"Data/",data_flag,"CAM_LENGTHS.csv")) %>%
			      .[,LENGTH_CM:=round(MEAN_MM/10)] %>%
				  .[,.(DROP_CD,SPECIES_CD,LENGTH_CM)] %>%
				  .[,.N,by=.(DROP_CD,SPECIES_CD,LENGTH_CM)] %>%
				  .[,TOTAL_N:=sum(N),by=.(DROP_CD,SPECIES_CD)] %>%
				  .[,PROP_N:=N/TOTAL_N] %>%
				  .[,.(DROP_CD,SPECIES_CD,LENGTH_CM,PROP_N)]

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
			  # keep all sizes for biomass
			  .[Length_category!="Unknown"] %>%
			  .[,.(DROP_CD,SPECIES_CD,N,KG,STD_N,STD_KG)]

	BFISH_CAM_C_long = copy(BFISH_CAM_C)
	BFISH_CAM_C = BFISH_CAM_C %>% 
			  .[,.(DROP_CD,SPECIES_CD,KG)] %>%
			  # keep only biomass measurement
			  dcast(.,DROP_CD~SPECIES_CD,value.var="KG",fill=0,fun.aggregate=sum) %>%
			  .[,.(DROP_CD,ETCO,ETCA,PRSI,PRFI,PRZO,HYQU,APRU)]
			

	# drop-specific information
	BFISH_CAM_S = fread(paste0(proj.dir,"Data/",data_flag,"CAM_SAMPLE_TIME.csv"))
	BFISH_CAM_S = BFISH_CAM_S %>%
			  .[,.(DROP_CD,DROP_DATE,DROP_TIME,VESSEL,GEAR_ID,PSU,OBS_LON,OBS_LAT,OFFICIAL_DEPTH_M,OFFICIAL_TEMP_C)] %>%
			  .[,c("SAMPLE_DATE","DROP_TIME_HST"):= convert_date_time(DROP_DATE,DROP_TIME),by=seq_len(nrow(BFISH_CAM_S))] %>%
			  .[,YEAR:=format(SAMPLE_DATE,format="%Y")] %>%
			  .[,MONTH:=format(SAMPLE_DATE,format="%m")] %>%
			  .[,DAY:=format(SAMPLE_DATE,format="%d")] %>%
			  .[,JD:=as.numeric(format(SAMPLE_DATE,format="%j"))] %>%
			  .[,YEAR_continuous:=as.numeric(YEAR)+(JD-1)/366] %>%
			  .[,LUNAR_PHASE:=getMoonIllumination(format(SAMPLE_DATE, format="%Y-%m-%d"))$fraction] %>%
			  .[,.(DROP_CD,SAMPLE_DATE,YEAR,MONTH,DAY,JD,YEAR_continuous,LUNAR_PHASE,DROP_TIME_HST,VESSEL,PSU,OBS_LON,OBS_LAT,OFFICIAL_DEPTH_M,OFFICIAL_TEMP_C)]
	
	# recalculate PSU based on actual location of drop
	close_psu_vec = rep(NA,nrow(BFISH_CAM_S))
	for(i in seq_along(close_psu_vec))
	{
		tmp_lon = BFISH_CAM_S$OBS_LON[i]
		tmp_lat = BFISH_CAM_S$OBS_LAT[i]

		tmp_dist = geosphere::distHaversine(c(tmp_lon,tmp_lat),as.matrix(PSU_table[,.(lon_deg,lat_deg)]))
		close_psu_vec[i] = PSU_table$PSU[which(tmp_dist==min(tmp_dist))]
		rm(list=c("tmp_lon","tmp_lat","tmp_dist"))
	}
	mean(BFISH_CAM_S$PSU == close_psu_vec)
	BFISH_CAM_S$PSU = close_psu_vec

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
						  na.omit(.) %>%
						  .[!(DROP_CD %in% dark_drops)] %>%
						  .[!(DROP_CD %in% missing_lengths)] %>%
						  # define sampling_unit variable
						  .[,SEASON:=ifelse(as.numeric(MONTH)>6,"Fall","Spring")] %>%
						  .[,design_sampling_unit:=paste0(YEAR,"_",SEASON,"_",PSU)] %>%
						  .[,model_sampling_unit:=DROP_CD]
	fwrite(unique(camera_dt[,.(DROP_CD)]),file=paste0(proj.dir,"Data/",data_flag,"01.camera_dt.all_lengths.drop_cd.csv"))

	# save formatted data
		save(BFISH_CAM_S,file=paste0(proj.dir,"Data/",data_flag,"01.BFISH_CAM_S.all_lengths.RData"))
		save(BFISH_CAM_C,file=paste0(proj.dir,"Data/",data_flag,"01.BFISH_CAM_C.all_lengths.RData"))
		save(dark_drops,file=paste0(proj.dir,"Data/",data_flag,"01.camera.dark_drops.all_lengths.RData"))
		save(missing_lengths,file=paste0(proj.dir,"Data/",data_flag,"01.camera.missing_lengths.all_lengths.RData"))
		save(camera_dt,file=paste0(proj.dir,"Data/",data_flag,"01.camera_dt.all_lengths.RData"))
