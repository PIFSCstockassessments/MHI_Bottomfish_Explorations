

# Nicholas Ducharme-Barth
# 09/23/2022
# Match BFISH to PacIOOS surface wind data
# Code is heavily influenced by John Syslo's code 'matching_windu_psu.R' and 'prep wind for r.R'
# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(ncdf4)
	library(RANN)
	library(ggplot2)
	library(ggthemes)

#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"

	data_flag = "2021_" # includes data through 2021
#_____________________________________________________________________________________________________________________________
# 1) bring in bfsish sampling data
	load(file=paste0(proj.dir,"Data/",data_flag,"bfish_combined_long_dt.RData"))
	# collapse weight_kg across bait_type
		sample_weight_dt = bfish_combined_long_dt[,.(weight_kg=sum(weight_kg)),by=.(sample_id,species_cd)]
		unique_dt = copy(bfish_combined_long_dt) %>% .[,bait_type:=NULL] %>% .[,weight_kg:=NULL] %>% unique(.)
		bfish_combined_long_dt = merge(unique_dt,sample_weight_dt,by=c("sample_id","species_cd"))
	
	# subset to species
	bfish_df = bfish_combined_long_dt %>% .[species_cd %in% c("prfi","etco","etca","prsi","przo","hyqu","apru")] %>% as.data.frame(.)

	# remove sample with large lehi observation
	sample_dt = as.data.table(bfish_df) %>%
			   .[,year:=as.numeric(as.character(year))] %>%
			   .[,season:=toupper(substr(season,1,1))] %>%
			   .[,wind_file:=paste0(year,"_",season)] %>%
			   .[!(year%in%c(2016,2017)),wind_file:=as.character(year)] %>%
			   .[,isodate:=as.numeric(gsub("-","",date))] %>%
   			   .[,.(sample_id,isodate,time,lon,lat,wind_file)] %>%
   			   unique(.)



#_____________________________________________________________________________________________________________________________
# 2) iterate across years and match with wind
# match based on closest spatial (lon/lat) and temporal (time) wind value
# time is hours since 2010-05-14 00-00-00.000 UTC for PacIOOS data
# wind velocity by component (u and v) at 10m altitude and given in m/s
# https://www.pacioos.hawaii.edu/weather/model-wind-hawaii/#about
# https://www.eol.ucar.edu/content/wind-direction-quick-reference
	
	u_windfile = sort(unique(sample_dt$wind_file))
	sample_wind_dt.list = as.list(rep(NA,length(u_windfile)))
	A = proc.time()
	for(i in 1:length(u_windfile))
	{
		# A = proc.time()

		# tmp_nc = nc_open(paste0(proj.dir,"Data/Wind/wrf_hi_WRF_Hawaii_Regional_Atmospheric_Model_best_TEST.nc"))
		# 	time_vec = ncvar_get(tmp_nc, "time")
		# nc_close(tmp_nc) 

		tmp_nc = nc_open(paste0(proj.dir,"Data/Wind/wrf_hi_WRF_Hawaii_Regional_Atmospheric_Model_best_",u_windfile[i],".nc"))
			lat_vec = ncvar_get(tmp_nc, "lat")
			lon_vec = ncvar_get(tmp_nc, "lon")
			time_vec = ncvar_get(tmp_nc, "time")
			wind_u_array = ncvar_get(tmp_nc, "Uwind")
			wind_v_array = ncvar_get(tmp_nc, "Vwind")
		nc_close(tmp_nc) 

		wind_u_dt = as.data.table(wind_u_array) %>%
					.[,V1:=lon_vec[V1]] %>%
					.[,V2:=lat_vec[V2]] %>%
					.[,V3:=time_vec[V3]] %>%
					setnames(.,c("V1","V2","V3","value"),c("lon","lat","time_utc","wind_u")) %>%
					.[lon >= range(sample_dt$lon)[1] & lon <= range(sample_dt$lon)[2]] %>%
					.[lat >= range(sample_dt$lat)[1] & lat <= range(sample_dt$lat)[2]]
		wind_v_dt = as.data.table(wind_v_array) %>%
					.[,V1:=lon_vec[V1]] %>%
					.[,V2:=lat_vec[V2]] %>%
					.[,V3:=time_vec[V3]] %>%
					setnames(.,c("V1","V2","V3","value"),c("lon","lat","time_utc","wind_v")) %>%
					.[lon >= range(sample_dt$lon)[1] & lon <= range(sample_dt$lon)[2]] %>%
					.[lat >= range(sample_dt$lat)[1] & lat <= range(sample_dt$lat)[2]]					
		wind_dt = merge(wind_u_dt,wind_v_dt,by=c("lon","lat","time_utc")) %>%
					.[,speed_ms:=sqrt(wind_u^2+wind_v^2 )]	%>%
					.[,speed_kt:=speed_ms*1.943844] %>%
					.[,direction:=(270-atan2(wind_v,wind_u)*180/pi) %% 360 ] %>%		
					.[,time_utc:=as.POSIXct(time_utc*60*60, origin = "2010-05-14", tz = "UTC")] %>%
					.[,time_hst:=as.POSIXct(format(time_utc, tz="HST",usetz=TRUE))] %>%
			  		.[,isodate:=as.numeric(format(time_hst,format="%Y%m%d"))] %>%
					.[,time:=as.numeric(as.character(format(time_hst,format="%H")))] %>%
					.[,st_id:=paste0(lon,"_",lat,"_",isodate,"_",time)]


		u_lon = unique(wind_dt$lon)
		u_lat = unique(wind_dt$lat)
		u_isodate = unique(wind_dt$isodate)
		u_time = unique(wind_dt$time)

		tmp_sample_dt = sample_dt[wind_file==u_windfile[i]]

		nn_lon = u_lon[nn2(u_lon,tmp_sample_dt$lon,k=1,searchtype="priority")$nn.idx]
		nn_lat = u_lat[nn2(u_lat,tmp_sample_dt$lat,k=1,searchtype="priority")$nn.idx]
		nn_isodate = u_isodate[nn2(u_isodate,tmp_sample_dt$isodate,k=1,searchtype="priority")$nn.idx]
		nn_time = u_time[nn2(u_time,tmp_sample_dt$time,k=1,searchtype="priority")$nn.idx]
		nn_id = paste0(nn_lon,"_",nn_lat,"_",nn_isodate,"_",nn_time)

		# tmp_sample_dt
		tmp_match_dt = wind_dt[match(nn_id,wind_dt$st_id),.(lon,lat,isodate,time,speed_ms,speed_kt,direction)]
		# tmp_match_dt
		sample_wind_dt.list[[i]] = cbind(tmp_sample_dt,tmp_match_dt[,.(speed_ms,speed_kt,direction)]) %>%
								   .[,direction_d:=as.character(NA)] %>%
								   .[direction>=337.5 | direction<22.5,direction_d :="N"] %>%
								   .[direction>=22.5 & direction<67.5,direction_d :="NE"] %>%
								   .[direction>=67.5 & direction<112.5,direction_d :="E"] %>%
								   .[direction>=112.5 & direction<157.5,direction_d :="SE"] %>%
								   .[direction>=157.5 & direction<202.5,direction_d :="S"] %>%
								   .[direction>=202.5 & direction<247.5,direction_d :="SW"] %>%
								   .[direction>=247.5 & direction<292.5,direction_d :="W"] %>%
								   .[direction>=292.5 & direction<337.5,direction_d :="NW"] %>%
								   .[,direction_d:=factor(direction_d,levels=c("N","NE","E","SE","S","SW","W","NW"))]

		# p = sample_wind_dt.list[[i]] %>%
		# 	ggplot() +
		# 	geom_point(aes(x=lon,y=lat,color=speed_kt)) +
		# 	viridis::scale_color_viridis("Wind speed (kt)",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE)

		# p

		# B = proc.time()
		# B-A

		rm(list=c("tmp_nc","lat_vec","lon_vec","time_vec","wind_u_array","wind_v_array","wind_u_dt","wind_v_dt","wind_dt","u_lon","u_lat","u_isodate","u_time","tmp_sample_dt","nn_lon","nn_lat","nn_isodate","nn_time","nn_id","tmp_match_dt"))			
	}
	B = proc.time()
	B-A

	sample_wind_dt = rbindlist(sample_wind_dt.list)
	save(sample_wind_dt,file=paste0(proj.dir,"Data/sample_wind_dt.RData"))

