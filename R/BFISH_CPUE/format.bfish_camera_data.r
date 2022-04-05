

# Nicholas Ducharme-Barth
# 04/04/2022
# Format BFISH camera data
# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
	setwd(proj.dir)
	source("D:/HOME/SAP/Code/Utilities/turbo.r")

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
	# drop-specific information
	BFISH_CAM_S = fread(paste0(proj.dir,"Data/CAM_SAMPLE_TIME.csv")) %>%
			  .[,.(DROP_CD,DROP_DATE,DROP_TIME,VESSEL,GEAR_ID,PSU,OBS_LON,OBS_LAT,OFFICIAL_DEPTH_M,OFFICIAL_TEMP_C)] %>%
			  .[,DROP_TIME_HST:=sapply(DROP_TIME,convert_drop_time)] %>%
			  .[,SAMPLE_DATE:=as.POSIXct(as.character(DROP_DATE),format=c("%Y%m%d"))] %>%
			  .[,YEAR:=format(SAMPLE_DATE,format="%Y")] %>%
			  .[,MONTH:=format(SAMPLE_DATE,format="%m")] %>%
			  .[,DAY:=format(SAMPLE_DATE,format="%m")] %>%
			  .[,JD:=format(SAMPLE_DATE,format="%j")] %>%
			  .[,YEAR_continuous:=as.numeric(YEAR)+(as.numeric(JD)-1)/366] %>%
			  .[,LUNAR_PHASE:=getMoonIllumination(format(SAMPLE_DATE, format="%Y-%m-%d"))$fraction] %>%
			  .[,SAMPLE_ID:=paste0(PSU,"_",DROP_DATE)] %>%
			  .[,.(SAMPLE_ID,DROP_CD,SAMPLE_DATE,YEAR,MONTH,DAY,JD,YEAR_continuous,LUNAR_PHASE,DROP_TIME_HST,VESSEL,GEAR_ID,PSU,OBS_LON,OBS_LAT,OFFICIAL_DEPTH_M,OFFICIAL_TEMP_C)]

			# need to manually edit SAMPLE_ID to combine samples across multiple days
			mult_PSU = BFISH_CAM_S[,.N,by=PSU][N>2]$PSU
			as.data.frame(BFISH_CAM_S[PSU %in% mult_PSU][order(PSU),.(SAMPLE_ID,PSU,SAMPLE_DATE)])
				# 20   8304_20170312  8304  2017-03-12
				# 21   8304_20170313  8304  2017-03-13

				# 28  20261_20181029 20261  2018-10-29
				# 29  20261_20181029 20261  2018-10-29
				# 30  20261_20181030 20261  2018-10-30

				# 54  25772_20201012 25772  2020-10-12
				# 55  25772_20201012 25772  2020-10-12
				# 56  25772_20201013 25772  2020-10-13
				# 57  25772_20201013 25772  2020-10-13

				# 70  31830_20170315 31830  2017-03-15
				# 71  31830_20170315 31830  2017-03-15
				# 72  31830_20170316 31830  2017-03-16

				# 83  33107_20171117 33107  2017-11-17
				# 84  33107_20171118 33107  2017-11-18

				# 143 38422_20161019 38422  2016-10-19
				# 144 38422_20161020 38422  2016-10-20
		  