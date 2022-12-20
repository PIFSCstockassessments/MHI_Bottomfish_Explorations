

# Nicholas Ducharme-Barth
# 12/19/2022
# Format BFISH camera data
# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

# This is a sensitivity analysis to see how much of a difference it makes to PSU level biomass if: 
# a) biomass expansion is done at the ssu level and then averaged OR
# b) counts are averaged at the PSU level and then biomass expansion is done with total psu length composition

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(suncalc)
	library(ggplot2)
	library(ggthemes)

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
# determine which camera drops were unable to record lengths for some species
	# bring in conversion factor data
		species_dt = fread(paste0(proj.dir,"Data/SPECIES_TABLE.csv")) %>%
				 .[SPECIES_CD %in% c("ETCO","ETCA","PRSI","PRFI","PRZO","HYQU","APRU")] %>%
				 unique(.)

		count_a = fread(paste0(proj.dir,"Data/",data_flag,"CAM_MAXN.csv")) %>%
			  .[,.(MAXN=sum(MAXN)),by=.(DROP_CD,SPECIES_CD)] %>%
			  .[SPECIES_CD %in% c("ETCO","ETCA","PRSI","PRFI","PRZO","HYQU","APRU")]

		lengths_a = fread(paste0(proj.dir,"Data/",data_flag,"CAM_LENGTHS.csv")) %>%
			      .[,LENGTH_CM:=round(MEAN_MM/10)] %>%
				  .[,.(DROP_CD,SPECIES_CD,LENGTH_CM)] %>%
				  .[,.N,by=.(DROP_CD,SPECIES_CD,LENGTH_CM)] %>%
				  .[,TOTAL_N:=sum(N),by=.(DROP_CD,SPECIES_CD)] %>%
				  .[,PROP_N:=N/TOTAL_N] %>%
				  .[,.(DROP_CD,SPECIES_CD,LENGTH_CM,PROP_N)]

		bfish_a = merge(lengths_a,count_a,by=c("DROP_CD","SPECIES_CD"),all.x=TRUE,all.y=TRUE) %>%
			  	  merge(.,species_dt[,.(SPECIES_CD,A,B,GCF)],by="SPECIES_CD") %>%
			  .[,N:=MAXN*PROP_N] %>%
			  .[,Biomass:=A * (LENGTH_CM^B)] %>%
			  .[,KG:=Biomass*N] %>%
			  .[,Length_category:=as.character(NA)] %>%
			  .[LENGTH_CM>=29,Length_category:="Exploitable"] %>%
			  .[LENGTH_CM<29,Length_category:="Un-exploitable"] %>%
			  .[is.na(LENGTH_CM),Length_category:="Unknown"]
	# remove drops where length is missing
		missing_lengths = unique(bfish_a[Length_category == "Unknown"]$DROP_CD)

#_____________________________________________________________________________________________________________________________
# bring in camera data
	# PSU specific information
		PSU_table = fread(paste0(proj.dir,"Data/BFISH PSU lookup table.csv")) %>%
				.[,.(PSU,Island,lon_deg,lat_deg,STRATA,STRATA_2020,Depth_MEDIAN_m,med_slp,med_acr,BS_pct_over_136j,pctHB,pctHS)] %>%
				.[,substrate:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][1])] %>%
				.[,slope:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][2])]
	
	# remove these drops/PSUs
		dark_drops = fread(paste0(proj.dir,"Data/",data_flag,"CAM_MAXN.csv"))[SPECIES_CD == "DARK"]$DROP_CD

	# drop-specific information
		BFISH_CAM_S = fread(paste0(proj.dir,"Data/",data_flag,"CAM_SAMPLE_TIME.csv")) %>%
			  .[,.(DROP_CD,DROP_DATE,DROP_TIME,VESSEL,GEAR_ID,PSU,OBS_LON,OBS_LAT,OFFICIAL_DEPTH_M,OFFICIAL_TEMP_C)] %>%
			  .[,DROP_TIME_HST:=sapply(DROP_TIME,convert_drop_time)] %>%
			  .[,SAMPLE_DATE:=as.POSIXct(as.character(DROP_DATE),format=c("%Y%m%d"))] %>%
			  .[,YEAR:=format(SAMPLE_DATE,format="%Y")] %>%
			  .[,MONTH:=format(SAMPLE_DATE,format="%m")] %>%
			  .[,DAY:=format(SAMPLE_DATE,format="%d")] %>%
			  .[,JD:=as.numeric(format(SAMPLE_DATE,format="%j"))] %>%
			  .[,SEASON:=ifelse(as.numeric(MONTH)>6,"Fall","Spring")]
	
				# find closest PSU
					close_psu_vec = rep(NA,nrow(BFISH_CAM_S))
					for(i in seq_along(close_psu_vec))
					{
						tmp_lon = BFISH_CAM_S$OBS_LON[i]
						tmp_lat = BFISH_CAM_S$OBS_LAT[i]

						tmp_dist = geosphere::distHaversine(c(tmp_lon,tmp_lat),as.matrix(PSU_table[,.(lon_deg,lat_deg)]))
						close_psu_vec[i] = PSU_table$PSU[which(tmp_dist==min(tmp_dist))]
						rm(list=c("tmp_lon","tmp_lat","tmp_dist"))
					}

			  BFISH_CAM_S = BFISH_CAM_S %>%
			  .[,actual_psu:=close_psu_vec] %>%
			  .[,.(PSU,actual_psu,DROP_CD,SAMPLE_DATE,YEAR,JD,SEASON,DROP_TIME_HST,VESSEL,OBS_LON,OBS_LAT)] %>%
			  .[,sampling_unit:=paste0(YEAR,"_",SEASON,"_",actual_psu)] %>%
			  .[!(DROP_CD %in% dark_drops)] %>%
			  .[!(DROP_CD %in% missing_lengths)] %>%
			  .[,N:=.N,by=sampling_unit] %>%
			  .[,.(DROP_CD,sampling_unit)]
	
#_____________________________________________________________________________________________________________________________
# get biomass per sampling unit using method a)
	bfish_a = bfish_a %>%
			  .[,.(N=sum(N,na.rm=TRUE),KG=sum(KG,na.rm=TRUE)),by=.(DROP_CD,SPECIES_CD,Length_category)] %>%
			  # keep only exploitable "sized" biomass
			  .[Length_category=="Exploitable"] %>%
			  .[,.(DROP_CD,SPECIES_CD,N,KG)] %>%
			  .[DROP_CD %in% BFISH_CAM_S$DROP_CD] %>%
			  .[,.(DROP_CD,SPECIES_CD,KG)] %>%
			  # keep only biomass measurement
			  dcast(.,DROP_CD~SPECIES_CD,value.var="KG",fill=0,fun.aggregate=sum) %>%
			  .[,.(DROP_CD,ETCO,ETCA,PRSI,PRFI,PRZO,HYQU,APRU)]
	
	bfish_a_zero = data.table(DROP_CD=BFISH_CAM_S[!(DROP_CD %in% bfish_a$DROP_CD)]$DROP_CD) %>%
				   .[,ETCO:=0] %>%
				   .[,ETCA:=0] %>%
				   .[,PRSI:=0] %>%
				   .[,PRFI:=0] %>%
				   .[,PRZO:=0] %>%
				   .[,HYQU:=0] %>%
				   .[,APRU:=0]

	weight_per_psu_a = rbind(bfish_a,bfish_a_zero) %>%
					   merge(.,BFISH_CAM_S,by="DROP_CD") %>%
					   .[,.(ETCO=mean(ETCO),ETCA=mean(ETCA),PRSI=mean(PRSI),PRFI=mean(PRFI),PRZO=mean(PRZO),HYQU=mean(HYQU),APRU=mean(APRU)),by=sampling_unit] %>%
					   melt(.,id.vars="sampling_unit") %>%
					   setnames(.,"value","A")



#_____________________________________________________________________________________________________________________________
# get biomass per sampling unit using method b)
		count_b_nzero = fread(paste0(proj.dir,"Data/",data_flag,"CAM_MAXN.csv")) %>%
			  .[,.(MAXN=sum(MAXN)),by=.(DROP_CD,SPECIES_CD)] %>%
			  .[DROP_CD %in% BFISH_CAM_S$DROP_CD] %>%
			  .[SPECIES_CD %in% c("ETCO","ETCA","PRSI","PRFI","PRZO","HYQU","APRU")] %>%
			  dcast(.,DROP_CD~SPECIES_CD,value.var="MAXN",fill=0,fun.aggregate=sum) %>%
			  .[,.(DROP_CD,ETCO,ETCA,PRSI,PRFI,PRZO,HYQU,APRU)]

		count_b_zero = data.table(DROP_CD=BFISH_CAM_S[!(DROP_CD %in% count_b_nzero$DROP_CD)]$DROP_CD) %>%
				   .[,ETCO:=0] %>%
				   .[,ETCA:=0] %>%
				   .[,PRSI:=0] %>%
				   .[,PRFI:=0] %>%
				   .[,PRZO:=0] %>%
				   .[,HYQU:=0] %>%
				   .[,APRU:=0]
		
		count_b = rbind(count_b_nzero,count_b_zero) %>%
				  merge(.,BFISH_CAM_S,by="DROP_CD") %>%
				  .[,.(ETCO=mean(ETCO),ETCA=mean(ETCA),PRSI=mean(PRSI),PRFI=mean(PRFI),PRZO=mean(PRZO),HYQU=mean(HYQU),APRU=mean(APRU)),by=sampling_unit] %>%
				  melt(.,id.vars="sampling_unit") %>%
				  setnames(.,c("variable","value"),c("SPECIES_CD","MAXN"))				  

		lengths_b = fread(paste0(proj.dir,"Data/",data_flag,"CAM_LENGTHS.csv")) %>%
			      .[,LENGTH_CM:=round(MEAN_MM/10)] %>%
				  .[,.(DROP_CD,SPECIES_CD,LENGTH_CM)] %>%
				  .[,.N,by=.(DROP_CD,SPECIES_CD,LENGTH_CM)] %>%
				  .[DROP_CD %in% BFISH_CAM_S$DROP_CD] %>%
				  merge(.,BFISH_CAM_S,by="DROP_CD") %>%
				  .[,TOTAL_N:=sum(N),by=.(sampling_unit,SPECIES_CD)] %>%
				  .[,PROP_N:=N/TOTAL_N] %>%
				  .[,.(sampling_unit,SPECIES_CD,LENGTH_CM,PROP_N)] %>%
				  .[order(sampling_unit,SPECIES_CD,LENGTH_CM)]


		bfish_b = merge(lengths_b,count_b,by=c("sampling_unit","SPECIES_CD"),all.x=TRUE,all.y=TRUE) %>%
			  	  merge(.,species_dt[,.(SPECIES_CD,A,B,GCF)],by="SPECIES_CD") %>%
			  .[,N:=MAXN*PROP_N] %>%
			  .[,Biomass:=A * (LENGTH_CM^B)] %>%
			  .[,KG:=Biomass*N] %>%
			  .[,Length_category:=as.character(NA)] %>%
			  .[LENGTH_CM>=29,Length_category:="Exploitable"] %>%
			  .[LENGTH_CM<29,Length_category:="Un-exploitable"] %>%
			  .[is.na(LENGTH_CM),Length_category:="Unknown"]

	bfish_b = bfish_b %>%
			  .[,.(N=sum(N,na.rm=TRUE),KG=sum(KG,na.rm=TRUE)),by=.(sampling_unit,SPECIES_CD,Length_category)] %>%
			  # keep only exploitable "sized" biomass
			  .[Length_category=="Exploitable"] %>%
			  .[,.(sampling_unit,SPECIES_CD,KG)] %>%
			  # keep only biomass measurement
			  dcast(.,sampling_unit~SPECIES_CD,value.var="KG",fill=0,fun.aggregate=sum) %>%
			  .[,.(sampling_unit,ETCO,ETCA,PRSI,PRFI,PRZO,HYQU,APRU)]
	
	bfish_b_zero = data.table(sampling_unit=unique(BFISH_CAM_S[!(sampling_unit %in% bfish_b$sampling_unit)]$sampling_unit)) %>%
				   .[,ETCO:=0] %>%
				   .[,ETCA:=0] %>%
				   .[,PRSI:=0] %>%
				   .[,PRFI:=0] %>%
				   .[,PRZO:=0] %>%
				   .[,HYQU:=0] %>%
				   .[,APRU:=0]

	weight_per_psu_b = rbind(bfish_b,bfish_b_zero) %>%
					   melt(.,id.vars="sampling_unit") %>%
					   setnames(.,"value","B")


#_____________________________________________________________________________________________________________________________
# combine and plot
	weight_per_psu = merge(weight_per_psu_a,weight_per_psu_b,by=c("sampling_unit","variable")) %>%
					.[,diff:=abs(B-A)]

	p = weight_per_psu %>%
		.[A>0 | B> 0] %>%
		ggplot() + 
		xlim(0,50) + 
		ylim(0,50) +
		facet_wrap(~variable) +
		xlab("Biomass (KG) method A") +
		ylab("Biomass (KG) method B") +
		geom_abline(slope=1,intercept=0,linetype="dashed",color="gray70") +
		geom_point(aes(x=A,y=B,fill=diff),size=2,shape=21) +
		theme_few(base_size = 20) +
		viridis::scale_fill_viridis("Difference\nbetween\nMethods B\nand A",begin = 0.1,end = 0.8,direction = 1,option = "H")
	p

	weight_per_psu[B!=A]

	BFISH_CAM_S[sampling_unit == "2017_Spring_18"]$DROP_CD
	bfish_a[DROP_CD %in% BFISH_CAM_S[sampling_unit == "2017_Spring_18"]$DROP_CD]
	count_a[DROP_CD %in% BFISH_CAM_S[sampling_unit == "2017_Spring_18"]$DROP_CD&SPECIES_CD=="APRU"]
	lengths_a[DROP_CD %in% BFISH_CAM_S[sampling_unit == "2017_Spring_18"]$DROP_CD&SPECIES_CD=="APRU"]

	bfish_b[sampling_unit == "2017_Spring_18"]
	count_b[sampling_unit == "2017_Spring_18"&SPECIES_CD=="APRU"]
	lengths_b[sampling_unit == "2017_Spring_18"&SPECIES_CD=="APRU"]


# doesn't look like there is a major difference between two approaches
# when there is it looks like option B can result in greater biomass (see above for Lehi)

# the key issue appears to be how to deal with PSUs where one of the SSUs get dropped (i.e., due to missing length measurements or darkness)
# do you throw out the whole PSU or do some PSUs only have one SSU in them ?