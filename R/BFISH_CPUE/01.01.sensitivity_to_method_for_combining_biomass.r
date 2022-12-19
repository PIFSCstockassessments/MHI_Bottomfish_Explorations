

# Nicholas Ducharme-Barth
# 12/19/2022
# Format BFISH camera data
# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

# This is a sensitivity analysis to see how much of a difference it makes to PSU level biomass if: 
# a) biomass expansion is done at the ssu level and then averaged OR
# b) counts are averaged at the PSU level and then biomass expansion is done with total psu length composition

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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