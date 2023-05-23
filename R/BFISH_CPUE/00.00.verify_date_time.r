

# Nicholas Ducharme-Barth
# 22/05/2023
# confirm that we can just use the date/time from DROP_CD
# DROP_CD records the date/time in UTC
# Copyright (c) 2023 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"

#_____________________________________________________________________________________________________________________________
# bring in data
	BFISH_CAM_S = read.csv(paste0(proj.dir,"Data/2022_CAM_SAMPLE.csv"),stringsAsFactors=FALSE)

	convert_time = function(x)
	{
			if(is.na(x))
			{
				time=NA
			} else if (nchar(x)==3){
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

			return(time)
	}

	convert_time2 = function(x)
	{
			if(is.na(x))
			{
				time=NA
			} else if (nchar(x)==3){
				hour = 0
				min = as.numeric(substr(x,1,1))/60
				time = hour + min 
			} else if(nchar(x)==4){
				hour = 0
				min = as.numeric(substr(x,1,2))/60
				time = hour + min 
			} else if(nchar(x)==5){
				hour = as.numeric(substr(x,1,1))
				min = as.numeric(substr(x,2,3))/60
				time = hour + min 
			} else if(nchar(x)==6){
				hour = as.numeric(substr(x,1,2))
				min = as.numeric(substr(x,3,4))/60
				time = hour + min 
			} else {
				time = NA
			}

			if(!is.na(time))
			{
				time = time - 10
				if(time<0)
				{
					time = 24+time
				}
			}

			return(time)
	}

	convert_time3 = function(x)
	{
			if(is.na(x))
			{
				time=NA
			} else if (nchar(x)==3){
				hour = as.numeric(substr(x,1,1))
				min = as.numeric(substr(x,2,3))/(60)
				time = hour + min
			} else if(nchar(x)==4){
				hour = as.numeric(substr(x,1,2))
				min = as.numeric(substr(x,3,4))/(60)
				time = hour + min
			} else {
				time = NA
			}

			return(time)
	}

	BFISH_CAM_S$test = sapply(as.numeric(sapply(BFISH_CAM_S$DROP_CD,function(x)strsplit(x,"_")[[1]][2])),convert_time)
	BFISH_CAM_S$test2 = sapply(as.numeric(gsub(":","",BFISH_CAM_S$DROP_TIME)),convert_time)

	# pull out cases where the DROP_CD time doesn't match-up with DROP_TIME
	new_data = subset(BFISH_CAM_S,test!=test2)
	new_data$test = sapply(as.numeric(sapply(new_data$DROP_CD,function(x)strsplit(x,"_")[[1]][2])),convert_time2)
	new_data$test2 = sapply(as.numeric(gsub(":","",new_data$DROP_TIME)),convert_time3)
	# only 3 cases where the times are greater than 30 minutes different
	# all from Ao Shibi IV
	# 2 cases in recent years where it looks like GMT was recorded instead of HST for DROP_TIME
	# 1 case where the difference was ~50 minutes and the time from DROP_CD was later
	# all cases would suggest that times derived from DROP_CD are more consistent and reliable
	
