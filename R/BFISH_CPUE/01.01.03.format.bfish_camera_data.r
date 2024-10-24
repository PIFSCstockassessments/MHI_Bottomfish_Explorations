

# Nicholas Ducharme-Barth
# 12/20/2022
# Format BFISH camera data: Option 3
# The sampling_unit is defined as two camera drops within a PSU during the same sampling season
# Sampling units are excluded if they contain drops that are "dark" or if they have missing length measurements
# Values are averaged within sampling unit (PSU). 
# Note: The number of camera drops within sampling unit may not be equal to 2...
# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)

#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"

#_____________________________________________________________________________________________________________________________
# define data_flag
	# data_flag = "" # only loads data up through 2020
	data_flag = "2021_" # includes data through 2021

#_____________________________________________________________________________________________________________________________
# bring in camera data
	# PSU specific information
	PSU_table = fread(paste0(proj.dir,"Data/BFISH PSU lookup table.csv")) %>%
				.[,.(PSU,Island,lon_deg,lat_deg,STRATA,STRATA_2020,Depth_MEDIAN_m,med_slp,med_acr,BS_pct_over_136j,pctHB,pctHS)] %>%
				.[,substrate:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][1])] %>%
				.[,slope:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][2])]

	load(file=paste0(proj.dir,"Data/",data_flag,"01.camera.dark_drops.RData"))
	load(file=paste0(proj.dir,"Data/",data_flag,"01.camera.missing_lengths.RData"))
    load(file=paste0(proj.dir,"Data/",data_flag,"01.BFISH_CAM_S.RData"))
	load(file=paste0(proj.dir,"Data/",data_flag,"01.BFISH_CAM_C.RData"))

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
						  # define sampling_unit variable
						  .[,SEASON:=ifelse(as.numeric(MONTH)>6,"Fall","Spring")] %>%
						  .[,design_sampling_unit:=paste0(YEAR,"_",SEASON,"_",PSU)] %>%
						  .[,model_sampling_unit:=paste0(YEAR,"_",SEASON,"_",PSU)]
    
    # remove sampling units that were dark or impacted by missing length measurements
    bad_sampling_units = camera_dt[,.(model_sampling_unit,DROP_CD)] %>%
                         .[(DROP_CD %in% dark_drops) | (DROP_CD %in% missing_lengths)]

    camera_dt = camera_dt %>% .[!(model_sampling_unit %in% bad_sampling_units$model_sampling_unit)]

    	                    

    # combine secondary sampling units (DROP_CD) within each primary sampling unit (model_sampling_unit)
    mode = function(x)
    # grabs the first item in the vector in the case of a tie
    {
        tmp_class = class(x)[1]
        if(!(tmp_class %in% c("integer","numeric")))
        {
            x = as.character(x)
        }
        unique_x = unique(x)
        return(unique_x[which.max(tabulate(match(x, unique_x)))])
    }

    # we will remove the following columns:
    # DROP_CD
    camera_dt = camera_dt %>%
                .[order(YEAR_continuous,DROP_TIME_HST,model_sampling_unit)] %>%
                .[,.(PSU=mode(PSU),SAMPLE_DATE=mode(SAMPLE_DATE),YEAR=mode(YEAR),MONTH=mode(MONTH),DAY=mode(DAY),JD=round(mean(JD)),YEAR_continuous=mean(YEAR_continuous),LUNAR_PHASE=mean(LUNAR_PHASE),DROP_TIME_HST=mean(DROP_TIME_HST),VESSEL=mode(VESSEL),OBS_LON=mean(OBS_LON),OBS_LAT=mean(OBS_LAT),OFFICIAL_DEPTH_M=mean(OFFICIAL_DEPTH_M),OFFICIAL_TEMP_C=mean(OFFICIAL_TEMP_C),ETCO=mean(ETCO),ETCA=mean(ETCA),PRSI=mean(PRSI),PRFI=mean(PRFI),PRZO=mean(PRZO),HYQU=mean(HYQU),APRU=mean(APRU),Island=mode(Island),STRATA=mode(STRATA),STRATA_2020=mode(STRATA_2020),Depth_MEDIAN_m=mean(Depth_MEDIAN_m),substrate=mode(substrate),slope=mode(slope),med_slp=mean(med_slp),med_acr=mean(med_acr),BS_pct_over_136j=mean(BS_pct_over_136j),pctHB=mean(pctHB),pctHS=mean(pctHS),depth_strata=mode(depth_strata),depth_strata_2020=mode(depth_strata_2020),complexity=mode(complexity),hardness=mode(hardness),SEASON=mode(SEASON),design_sampling_unit=mode(design_sampling_unit)),by=model_sampling_unit]


	# save formatted data
		save(camera_dt,file=paste0(proj.dir,"Data/",data_flag,"03.camera_dt.RData"))
