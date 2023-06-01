

# Nicholas Ducharme-Barth
# 2023/06/01
# Combine research_fishing_dt and camera_dt data sets into a combined analysis ready data set
# define common columns
# id: model_sampling_unit, design_sampling_unit
# strata: psu, island, strata, strata_2020, substrate, slope, depth_strata, depth_strata_2020, complexity, hardness
# temporal: date, year, season, month, jd, year_continuous
# spatial: lon, lat
# covar: depth, time, lunar_phase
# q_covar: gear_type, platform, obs_wind, obs_wave, obs_current
# catch: etco, etca, prsi, prfi, przo, hyqu, apru
# Copyright (c) 2023 Nicholas Ducharme-Barth
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
	data_flag = "2022_" # includes data through 2022
#_____________________________________________________________________________________________________________________________
# bring research fishing data
	load(file=paste0(proj.dir,"Data/",data_flag,"research_fishing_dt.RData"))

#_____________________________________________________________________________________________________________________________
# format research_fishing_dt to have common columns
	common_research_fishing_dt = research_fishing_dt %>% 
					   .[,gear_type:="research_fishing"] %>%
					   .[,season:="fall"] %>%
					   .[MONTH %in% c("02","03"),season:="spring"] %>%
					   # anonymize vessel designation
					   .[,VESSEL:=LETTERS[as.numeric(as.factor(VESSEL))]] %>%
					   # rename variables
					   setnames(.,c("SAMPLE_ID",
					   				"PSU","Island","STRATA","STRATA_2020","substrate","slope","depth_strata","depth_strata_2020","complexity","hardness",
					   				"SAMPLE_DATE","YEAR","season","MONTH","JD","YEAR_continuous",
					   				"LON","LAT",
					   				"DEPTH_M","TIME_MIN","LUNAR_PHASE",
					   				"gear_type", "VESSEL","WIND_SPEED_KT", "WAVE_HEIGHT_FT", "CURRENT_CD",
					   				"ETCO","ETCA","PRSI","PRFI","PRZO","HYQU","APRU"),
					   				c("sample_id",
					   				  "psu", "island", "strata", "strata_2020", "substrate", "slope", "depth_strata", "depth_strata_2020", "complexity", "hardness",
					   				  "date", "year","season", "month", "jd", "year_continuous",
					   				  "lon", "lat",
					   				  "depth", "time", "lunar_phase",
					   				  "gear_type", "platform", "obs_wind", "obs_wave", "obs_current",
					   				  "etco", "etca", "prsi", "prfi", "przo", "hyqu", "apru")) %>%
					   # select proper columns
					   .[,.(model_sampling_unit, design_sampling_unit,psu, island, strata, strata_2020, substrate, slope, depth_strata, depth_strata_2020, complexity, hardness,date, year, season, month, jd, year_continuous,lon, lat,depth, time, lunar_phase,gear_type, platform, obs_wind, obs_wave, obs_current, etco, etca, prsi, prfi, przo, hyqu, apru)]
	# save
	save(common_research_fishing_dt,file=paste0(proj.dir,"Data/",data_flag,"common_research_fishing_dt.RData"))


#_____________________________________________________________________________________________________________________________
# loop over 5 different camera data options, format, combine, and save
	combined_dt.list = as.list(rep(NA,5))

	for(i in 1:5)
	{
		# load data
		load(file=paste0(proj.dir,"Data/",data_flag,"0",i,".camera_dt.RData"))
		# format camera_dt to have common columns
		common_camera_dt = camera_dt %>% 
						.[,gear_type:="camera"] %>%
						.[,season:="fall"] %>%
						.[MONTH %in% c("02","03"),season:="spring"] %>%
						.[,obs_wind:=NA] %>%
						.[,obs_wave:=NA] %>%
						.[,obs_current:=NA] %>% 
						# rename variables
						setnames(.,c(	"PSU","Island","STRATA","STRATA_2020","substrate","slope","depth_strata","depth_strata_2020","complexity","hardness",
										"SAMPLE_DATE","YEAR","season","MONTH","JD","YEAR_continuous",
										"OBS_LON","OBS_LAT",
										"OFFICIAL_DEPTH_M","DROP_TIME_HST","LUNAR_PHASE",
										"gear_type", "VESSEL",
										"ETCO","ETCA","PRSI","PRFI","PRZO","HYQU","APRU"),
										c("psu", "island", "strata", "strata_2020", "substrate", "slope", "depth_strata", "depth_strata_2020", "complexity", "hardness",
										"date", "year","season", "month", "jd", "year_continuous",
										"lon", "lat",
										"depth", "time", "lunar_phase",
										"gear_type", "platform",
										"etco", "etca", "prsi", "prfi", "przo", "hyqu", "apru")) %>%
						# select proper columns
					   .[,.(model_sampling_unit, design_sampling_unit,psu, island, strata, strata_2020, substrate, slope, depth_strata, depth_strata_2020, complexity, hardness,date, year, season, month, jd, year_continuous,lon, lat,depth, time, lunar_phase,gear_type, platform, obs_wind, obs_wave, obs_current, etco, etca, prsi, prfi, przo, hyqu, apru)]
		# combine
		if(i>1)
		{
			# match class due to aggregation
			common_research_fishing_dt$date = as.character(common_research_fishing_dt$date)
		}
		combined_dt.list[[i]] = bfish_combined_wide_dt = rbind(common_camera_dt,common_research_fishing_dt) %>%
								.[order(year_continuous,gear_type,psu)]

		bfish_combined_long_dt = bfish_combined_wide_dt %>%
								melt(.,id.vars=c("model_sampling_unit", "design_sampling_unit",
										"psu", "island", "strata", "strata_2020", "substrate", "slope", "depth_strata", "depth_strata_2020", "complexity", "hardness",
										"date", "year", "season", "month", "jd", "year_continuous",
										"lon", "lat",
										"depth", "time", "lunar_phase",
										"gear_type", "platform", "obs_wind", "obs_wave", "obs_current")) %>%
								setnames(.,c("variable","value"),c("species_cd","weight_kg"))
		
		# save formatted data
		save(common_camera_dt,file=paste0(proj.dir,"Data/",data_flag,"0",i,".common_camera_dt.RData"))
		save(bfish_combined_wide_dt,file=paste0(proj.dir,"Data/",data_flag,"0",i,".bfish_combined_wide_dt.RData"))
		save(bfish_combined_long_dt,file=paste0(proj.dir,"Data/",data_flag,"0",i,".bfish_combined_long_dt.RData"))

		# clean-up
		rm(list=c("camera_dt","common_camera_dt","bfish_combined_wide_dt","bfish_combined_long_dt"))
	}
	
	# version 01 of the combined data contains all unique model_sampling_units & design_sampling_units
	# version 02 of the combined data contains all unique design_sampling_units
	# Note that in versions 02 -- 05, model_sampling_unit == design_sampling_unit
	# t = lapply(combined_dt.list,function(x)unique(x$design_sampling_unit))
	# mean(t[[1]] %in% t[[2]])
	# [1] 1
	# mean(t[[3]] %in% t[[2]])
	# [1] 1
	# mean(t[[4]] %in% t[[2]])
	# [1] 1
	# mean(t[[5]] %in% t[[2]])
	# [1] 1
	 
	# sapply(combined_dt.list,nrow)
	# [1] 4686 3872 3770 3680 3677






