

# Nicholas Ducharme-Barth
# 2023/05/22
# Combine 2016-2022 BFISH data
# also reformat DROP_DATE and DROP_TIME to be based on the DROP_CD column
# now both DROP_DATE and DROP_TIME are in UTC
# Copyright (c) 2023 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"

#_____________________________________________________________________________________________________________________________
# define helper function for reformating drop_time to old format so it doesn't break downstream analysis
	convert_drop_time = function(x)
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

	reformat_drift_time = function(x)
	{	
		if(x=="")
		{
			time = as.integer(NA)
		} else {
			time = as.integer(gsub(":","",x,fixed=TRUE))
			if(time<500)
			{
				time = time + 1200
			}
		}
		return(time)
	}
	
	reformat_drift_date = function(x)
	{
		return(as.character(as.POSIXct(x,format="%Y-%m-%d")))
	}

#_____________________________________________________________________________________________________________________________
# bring in data
	BFISH_CAM_COUNT = read.csv(paste0(proj.dir,"Data/2021_CAM_MAXN.csv"),stringsAsFactors=FALSE)
	BFISH_CAM_LENGTHS = read.csv(paste0(proj.dir,"Data/2021_CAM_LENGTHS.csv"),stringsAsFactors=FALSE)
	BFISH_CAM_S = read.csv(paste0(proj.dir,"Data/2021_CAM_SAMPLE_TIME.csv"),stringsAsFactors=FALSE)

	BFISH_CAM_COUNT_2022 = read.csv(paste0(proj.dir,"Data/CAM_MAXN_2022.csv"),stringsAsFactors=FALSE)
	BFISH_CAM_COUNT_2022 = subset(BFISH_CAM_COUNT_2022,BFISH=="BFISH_2022_F")
	BFISH_CAM_LENGTHS_2022 = read.csv(paste0(proj.dir,"Data/CAM_LENGTHS_2022.csv"),stringsAsFactors=FALSE)
	BFISH_CAM_LENGTHS_2022 = subset(BFISH_CAM_LENGTHS_2022,BFISH=="BFISH_2022_F")
	BFISH_CAM_S_2022 = read.csv(paste0(proj.dir,"Data/CAM_SAMPLE_2022.csv"),stringsAsFactors=FALSE)
	BFISH_CAM_S_2022 = subset(BFISH_CAM_S_2022,BFISH=="BFISH_2022_F")
		# re-format DROP_TIME
		BFISH_CAM_S_2022$DROP_TIME = sapply(as.numeric(sapply(BFISH_CAM_S_2022$DROP_CD,function(x)strsplit(x,"_")[[1]][2])),convert_drop_time)
		# re-format DROP_DATE
		BFISH_CAM_S_2022$DROP_DATE = as.integer(sapply(BFISH_CAM_S_2022$DROP_CD,function(x)strsplit(x,"_")[[1]][1]))
		# add new DROP_TIME_NEW column
		BFISH_CAM_S_2022$DROP_TIME_NEW = NA

	BFISH_D = read.csv(paste0(proj.dir,"Data/2021_CRF_DRIFT.csv"),stringsAsFactors=FALSE)
	BFISH_S = read.csv(paste0(proj.dir,"Data/2021_CRF_SAMPLE.csv"),stringsAsFactors=FALSE)
	BFISH_C = read.csv(paste0(proj.dir,"Data/2021_CRF_CATCH.csv"),stringsAsFactors=FALSE)

	BFISH_D_2022 = read.csv(paste0(proj.dir,"Data/CRF_DRIFT_2022.csv"),stringsAsFactors=FALSE)
	BFISH_D_2022 = subset(BFISH_D_2022,BFISH=="BFISH_2022_F")	
		# re-format drift times	
		BFISH_D_2022$DRIFT_START_TIME = sapply(BFISH_D_2022$DRIFT_START_TIME,reformat_drift_time)
		BFISH_D_2022$DRIFT_END_TIME = sapply(BFISH_D_2022$DRIFT_END_TIME,reformat_drift_time)
	BFISH_S_2022 = read.csv(paste0(proj.dir,"Data/CRF_SAMPLE_2022.csv"),stringsAsFactors=FALSE)
	BFISH_S_2022 = subset(BFISH_S_2022,BFISH=="BFISH_2022_F")	
		# re-format SAMPLE_DATE
		BFISH_S_2022$SAMPLE_DATE = sapply(BFISH_S_2022$SAMPLE_DATE,reformat_drift_date)
	
	BFISH_C_2022 = read.csv(paste0(proj.dir,"Data/CRF_CATCH_2022.csv"),stringsAsFactors=FALSE)
	BFISH_C_2022 = subset(BFISH_C_2022,BFISH=="BFISH_2022_F")	

#_____________________________________________________________________________________________________________________________
# define function to update data
	update = function(old,new)
	{
		original_old_cols = colnames(old)
		colnames(old) = toupper(colnames(old))
		colnames(new) = toupper(colnames(new))

		select_cols = colnames(old)
		missing_cols = setdiff(select_cols,colnames(new))
		# make pad
			pad_df = matrix(NA,nrow=nrow(new),ncol=length(missing_cols))
			colnames(pad_df) = missing_cols
			pad_df = as.data.frame(pad_df)
		# get original class of each missing column
			for(i in 1:length(missing_cols))
			{
				class(pad_df[,missing_cols[i]]) = class(old[,missing_cols[i]])
			}

		new_pad = cbind(new,pad_df)
		updated = rbind(old,new_pad[,select_cols])
		colnames(updated) = original_old_cols

		# report NAs
		na_matrix = matrix(0,nrow=3,ncol=length(select_cols))
		colnames(na_matrix) = original_old_cols
		rownames(na_matrix) = c("old","new","updated")
		for(i in 1:ncol(na_matrix))
		{
			na_matrix[1,i] = sum(ifelse(is.na(old[,i])|old[,i]=="NA",1,0))
			na_matrix[2,i] = sum(ifelse(is.na(new_pad[,select_cols][,i])|new_pad[,select_cols][,i]=="NA",1,0))
			na_matrix[3,i] = sum(ifelse(is.na(updated[,i])|updated[,i]=="NA",1,0))
		}

		print(missing_cols)
		print(na_matrix)
		return(updated)
	}

#_____________________________________________________________________________________________________________________________
# update data
	update.BFISH_CAM_COUNT = update(old=BFISH_CAM_COUNT,new=BFISH_CAM_COUNT_2022)
	update.BFISH_CAM_LENGTHS = update(old=BFISH_CAM_LENGTHS,new=BFISH_CAM_LENGTHS_2022)
	update.BFISH_CAM_S = update(old=BFISH_CAM_S,new=BFISH_CAM_S_2022)
	# re-format VESSEL
	update.BFISH_CAM_S$VESSEL = ifelse(update.BFISH_CAM_S$VESSEL=="Sette","Oscar Elton Sette",update.BFISH_CAM_S$VESSEL)

	update.BFISH_D = update(old=BFISH_D,new=BFISH_D_2022)
	update.BFISH_S = update(old=BFISH_S,new=BFISH_S_2022)
	update.BFISH_C = update(old=BFISH_C,new=BFISH_C_2022)

#_____________________________________________________________________________________________________________________________
# write-out
	write.csv(update.BFISH_CAM_COUNT,file=paste0(proj.dir,"Data/2022_CAM_MAXN.csv"))
	write.csv(update.BFISH_CAM_LENGTHS,file=paste0(proj.dir,"Data/2022_CAM_LENGTHS.csv"))
	write.csv(update.BFISH_CAM_S,file=paste0(proj.dir,"Data/2022_CAM_SAMPLE_TIME.csv"))

	write.csv(update.BFISH_D,file=paste0(proj.dir,"Data/2022_CRF_DRIFT.csv"))
	write.csv(update.BFISH_S,file=paste0(proj.dir,"Data/2022_CRF_SAMPLE.csv"))
	write.csv(update.BFISH_C,file=paste0(proj.dir,"Data/2022_CRF_CATCH.csv"))
	
