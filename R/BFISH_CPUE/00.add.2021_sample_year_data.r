

# Nicholas Ducharme-Barth
# 05/06/2022
# Combine 2016-2021 BFISH data
# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"

#_____________________________________________________________________________________________________________________________
# define helper function for reformating drop_time to old format so it doesn't break downstream analysis
	reformat_new_drop_time = function(x)
	{	
			hour = as.numeric(strsplit(x,":")[[1]][1])
			min = strsplit(x,":")[[1]][2]
		# convert to 0-24 UTC
			utc_hour = hour + 10
			if(utc_hour>=24)
			{
				utc_hour = utc_hour - 24
			}
			utc_hour = as.character(utc_hour)
			time = as.integer(paste0(utc_hour,min,"00"))
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

	reformat_drop_date = function(x)
	{
		return(gsub("-","",as.character(as.POSIXct(x,format="%m/%e/%Y"))))
	}
	
	reformat_sample_date = function(x)
	{
		return(as.character(as.POSIXct(x,format="%m/%e/%Y")))
	}

#_____________________________________________________________________________________________________________________________
# bring in data
	BFISH_CAM_COUNT = read.csv(paste0(proj.dir,"Data/CAM_MAXN.csv"),stringsAsFactors=FALSE)
	BFISH_CAM_LENGTHS = read.csv(paste0(proj.dir,"Data/CAM_LENGTHS.csv"),stringsAsFactors=FALSE)
	BFISH_CAM_S = read.csv(paste0(proj.dir,"Data/CAM_SAMPLE_TIME.csv"),stringsAsFactors=FALSE)

	BFISH_CAM_COUNT_2021 = read.csv(paste0(proj.dir,"Data/CAM_MAXN_2021.csv"),stringsAsFactors=FALSE)
	BFISH_CAM_LENGTHS_2021 = read.csv(paste0(proj.dir,"Data/CAM_LENGTHS_2021.csv"),stringsAsFactors=FALSE)
	BFISH_CAM_S_2021 = read.csv(paste0(proj.dir,"Data/CAM_SAMPLE_TIME_2021.csv"),stringsAsFactors=FALSE)
		# re-format DROP_TIME & convert to UTC
		BFISH_CAM_S_2021$DROP_TIME = sapply(BFISH_CAM_S_2021$DROP_TIME,reformat_new_drop_time)
		# add column for OFFICIAL_TEMP_C
		BFISH_CAM_S_2021$OFFICIAL_TEMP_C = BFISH_CAM_S_2021$TDR_temp_mean
		# re-format DROP_DATE
		BFISH_CAM_S_2021$DROP_DATE = sapply(BFISH_CAM_S_2021$DROP_DATE,reformat_drop_date)


	BFISH_D = read.csv(paste0(proj.dir,"Data/CRF_DRIFT.csv"),stringsAsFactors=FALSE)
	BFISH_S = read.csv(paste0(proj.dir,"Data/CRF_SAMPLE.csv"),stringsAsFactors=FALSE)
	BFISH_C = read.csv(paste0(proj.dir,"Data/CRF_CATCH.csv"),stringsAsFactors=FALSE)

	BFISH_D_2021 = read.csv(paste0(proj.dir,"Data/CRF_DRIFT_2021.csv"),stringsAsFactors=FALSE)
		# re-format drift times	
		BFISH_D_2021$DRIFT_START_TIME = sapply(BFISH_D_2021$DRIFT_START_TIME,reformat_drift_time)
		BFISH_D_2021$DRIFT_END_TIME = sapply(BFISH_D_2021$DRIFT_END_TIME,reformat_drift_time)
	BFISH_S_2021 = read.csv(paste0(proj.dir,"Data/CRF_SAMPLE_2021.csv"),stringsAsFactors=FALSE)
		# re-format SAMPLE_DATE
		BFISH_S_2021$SAMPLE_DATE = sapply(BFISH_S_2021$SAMPLE_DATE,reformat_sample_date)
	
	BFISH_C_2021 = read.csv(paste0(proj.dir,"Data/CRF_CATCH_2021.csv"),stringsAsFactors=FALSE)
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

		print(missing_cols)
		return(updated)
	}

#_____________________________________________________________________________________________________________________________
# update data
	update.BFISH_CAM_COUNT = update(old=BFISH_CAM_COUNT,new=BFISH_CAM_COUNT_2021)
	update.BFISH_CAM_LENGTHS = update(old=BFISH_CAM_LENGTHS,new=BFISH_CAM_LENGTHS_2021)
	update.BFISH_CAM_S = update(old=BFISH_CAM_S,new=BFISH_CAM_S_2021)
	# re-format VESSEL
	update.BFISH_CAM_S$VESSEL = ifelse(update.BFISH_CAM_S$VESSEL=="Oscar Sette","Oscar Elton Sette",update.BFISH_CAM_S$VESSEL)

	update.BFISH_D = update(old=BFISH_D,new=BFISH_D_2021)
	update.BFISH_S = update(old=BFISH_S,new=BFISH_S_2021)
	update.BFISH_C = update(old=BFISH_C,new=BFISH_C_2021)

#_____________________________________________________________________________________________________________________________
# write-out
	write.csv(update.BFISH_CAM_COUNT,file=paste0(proj.dir,"Data/2021_CAM_MAXN.csv"))
	write.csv(update.BFISH_CAM_LENGTHS,file=paste0(proj.dir,"Data/2021_CAM_LENGTHS.csv"))
	write.csv(update.BFISH_CAM_S,file=paste0(proj.dir,"Data/2021_CAM_SAMPLE_TIME.csv"))

	write.csv(update.BFISH_D,file=paste0(proj.dir,"Data/2021_CRF_DRIFT.csv"))
	write.csv(update.BFISH_S,file=paste0(proj.dir,"Data/2021_CRF_SAMPLE.csv"))
	write.csv(update.BFISH_C,file=paste0(proj.dir,"Data/2021_CRF_CATCH.csv"))
	