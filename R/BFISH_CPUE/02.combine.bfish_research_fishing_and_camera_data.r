

# Nicholas Ducharme-Barth
# 04/14/2022
# Combine research_fishing_dt and camera_dt data sets into a single analysis ready data set
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
# bring in data
	load(file=paste0(proj.dir,"Data/camera_dt.RData"))
	load(file=paste0(proj.dir,"Data/research_fishing_dt.RData"))

#_____________________________________________________________________________________________________________________________
# define common columns
# id: sample_id
# strata: psu, island, strata, strata_2020, substrate, slope, depth_strata, depth_strata_2020, complexity, hardness
# temporal: date, year, season, month, jd, year_continuous
# spatial: lon, lat
# covar: depth, time, lunar_phase
# q_covar: gear_type, platform, bait_type
# catch: etco, etca, prsi, prfi, przo, hyqu, apru

#_____________________________________________________________________________________________________________________________
# format camera_dt to have common columns
	common_camera_dt = camera_dt %>% 
					   .[,gear_type:="camera"] %>%
					   .[,bait_type:="C"] %>%
					   .[,season:="fall"] %>%
					   .[MONTH %in% c("02","03"),season:="spring"] %>% 
					   # rename variables
					   setnames(.,c("DROP_CD",
					   				"PSU","Island","STRATA","STRATA_2020","substrate","slope","depth_strata","depth_strata_2020","complexity","hardness",
					   				"SAMPLE_DATE","YEAR","season","MONTH","JD","YEAR_continuous",
					   				"OBS_LON","OBS_LAT",
					   				"OFFICIAL_DEPTH_M","DROP_TIME_HST","LUNAR_PHASE",
					   				"gear_type", "VESSEL", "bait_type",
					   				"ETCO","ETCA","PRSI","PRFI","PRZO","HYQU","APRU"),
					   				c("sample_id",
					   				  "psu", "island", "strata", "strata_2020", "substrate", "slope", "depth_strata", "depth_strata_2020", "complexity", "hardness",
					   				  "date", "year","season", "month", "jd", "year_continuous",
					   				  "lon", "lat",
					   				  "depth", "time", "lunar_phase",
					   				  "gear_type", "platform", "bait_type",
					   				  "etco", "etca", "prsi", "prfi", "przo", "hyqu", "apru")) %>%
					   # select proper columns
					   .[,.(sample_id,psu, island, strata, strata_2020, substrate, slope, depth_strata, depth_strata_2020, complexity, hardness,date, year, season, month, jd, year_continuous,lon, lat,depth, time, lunar_phase,gear_type, platform, bait_type,etco, etca, prsi, prfi, przo, hyqu, apru)]

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
					   				"gear_type", "VESSEL", "BAIT_CD",
					   				"ETCO","ETCA","PRSI","PRFI","PRZO","HYQU","APRU"),
					   				c("sample_id",
					   				  "psu", "island", "strata", "strata_2020", "substrate", "slope", "depth_strata", "depth_strata_2020", "complexity", "hardness",
					   				  "date", "year","season", "month", "jd", "year_continuous",
					   				  "lon", "lat",
					   				  "depth", "time", "lunar_phase",
					   				  "gear_type", "platform", "bait_type",
					   				  "etco", "etca", "prsi", "prfi", "przo", "hyqu", "apru")) %>%
					   # select proper columns
					   .[,.(sample_id,psu, island, strata, strata_2020, substrate, slope, depth_strata, depth_strata_2020, complexity, hardness,date, year, season, month, jd, year_continuous,lon, lat,depth, time, lunar_phase,gear_type, platform, bait_type,etco, etca, prsi, prfi, przo, hyqu, apru)]

#_____________________________________________________________________________________________________________________________
# combine
	bfish_combined_wide_dt = rbind(common_camera_dt,common_research_fishing_dt) %>%
							 .[order(year_continuous,gear_type,psu)]

	bfish_combined_long_dt = bfish_combined_wide_dt %>%
							 melt(.,id.vars=c("sample_id",
					   				  "psu", "island", "strata", "strata_2020", "substrate", "slope", "depth_strata", "depth_strata_2020", "complexity", "hardness",
					   				  "date", "year", "season", "month", "jd", "year_continuous",
					   				  "lon", "lat",
					   				  "depth", "time", "lunar_phase",
					   				  "gear_type", "platform", "bait_type")) %>%
							 setnames(.,c("variable","value"),c("species_cd","weight_kg"))
#_____________________________________________________________________________________________________________________________
# save formatted data
	save(common_camera_dt,file=paste0(proj.dir,"Data/common_camera_dt.RData"))
	save(common_research_fishing_dt,file=paste0(proj.dir,"Data/common_research_fishing_dt.RData"))
	save(bfish_combined_wide_dt,file=paste0(proj.dir,"Data/bfish_combined_wide_dt.RData"))
	save(bfish_combined_long_dt,file=paste0(proj.dir,"Data/bfish_combined_long_dt.RData"))

