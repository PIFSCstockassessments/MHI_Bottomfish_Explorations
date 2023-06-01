

# Nicholas Ducharme-Barth
# 2023/06/01
# Match BFISH to PacIOOS surface wind data
# Code is heavily influenced by John Syslo's code 'matching_windu_psu.R' and 'prep wind for r.R'
# Copyright (c) 2023 Nicholas Ducharme-Barth
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

	data_flag = "2022_" # includes data through 2022
#_____________________________________________________________________________________________________________________________
# 1) bring in bfsish sampling data
	load(file=paste0(proj.dir,"Data/",data_flag,"02.bfish_combined_long_dt.RData"))

	bfish_df = bfish_combined_long_dt %>% .[species_cd %in% c("prfi","etco","etca","prsi","przo","hyqu","apru")] %>% as.data.frame(.)

	sample_dt = as.data.table(bfish_df) %>%
			   .[,year:=as.numeric(as.character(year))] %>%
			   .[,season:=toupper(substr(season,1,1))] %>%
			   .[,wind_file:=paste0(year,"_",season)] %>%
			   .[!(year%in%c(2016,2017)),wind_file:=as.character(year)] %>%
			   .[,isodate:=as.numeric(gsub("-","",date))] %>%
   			   .[,.(design_sampling_unit,isodate,time,lon,lat,wind_file)] %>%
   			   unique(.)



#_____________________________________________________________________________________________________________________________
# 2) iterate across years and match with wind
# match based on closest spatial (lon/lat) and temporal (time) wind value
# time is hours since 2010-05-14 12-00-00 UTC for PacIOOS data
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
					.[,time_utc:=as.POSIXct(time_utc*60*60, origin = "2010-05-14 12:00:00", tz = "UTC")] %>%
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
		# for 7 missing observations from 2017-11-15 (PacIOOS has no predictions for 2017-11-16 UTC)
		if(u_windfile[i] == "2017_F")
		{
			missing_id = which(is.na(tmp_match_dt$speed_ms))
			nn_id[missing_id] = paste0(nn_lon[missing_id],"_",nn_lat[missing_id],"_",nn_isodate[missing_id],"_14")
			tmp_match_dt = wind_dt[match(nn_id,wind_dt$st_id),.(lon,lat,isodate,time,speed_ms,speed_kt,direction)]
			rm(list=c("missing_id"))
		}

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
								   .[,direction_d:=factor(direction_d,levels=c("N","NE","E","SE","S","SW","W","NW"))] %>%
								   .[,speed_mph:=speed_kt*1.15077945] %>%
								   .[,beaufort:=cut(speed_ms,breaks=c(0,0.5,1.6,3.4,5.5,8,10.8,13.9,17.2,20.8,24.5,28.5,32.7,100),labels=FALSE)-1]								   
		# calculate difference in distance and time between sample and wind observation
			sample_wind_dt.list[[i]]$temporal_diff_hours = abs(((as.numeric(as.POSIXct(as.character(tmp_sample_dt$isodate),format="%Y%m%d")) + tmp_sample_dt$time*60*60) -	(as.numeric(as.POSIXct(as.character(tmp_match_dt$isodate),format="%Y%m%d")) + tmp_match_dt$time*60*60))/(60*60))
			sample_wind_dt.list[[i]]$spatial_diff_km = geosphere::distHaversine(as.matrix(tmp_sample_dt[,.(lon,lat)]),as.matrix(tmp_match_dt[,.(lon,lat)]))/1000

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
	save(sample_wind_dt,file=paste0(proj.dir,"Data/",data_flag,"sample_wind_dt.RData"))


#_____________________________________________________________________________________________________________________________
# 3) make some summary plots
	load(file=paste0(proj.dir,"Data/",data_flag,"sample_wind_dt.RData"))
	plot_dir = paste0(proj.dir,"Plot/BFISH_CPUE/")

	wind_plot_dt = merge(bfish_combined_long_dt,sample_wind_dt[,.(design_sampling_unit,speed_ms,speed_kt,speed_mph,beaufort,direction,direction_d)],by="design_sampling_unit")

	p = copy(wind_plot_dt) %>%
		.[,bin:=ifelse(weight_kg>0,1,0)] %>%
		.[,wind_floor:=floor(speed_kt)] %>%
		.[,.(encounter=mean(bin),.N),by=.(species_cd,gear_type,wind_floor)] %>%
		.[,species_cd:=factor(species_cd,levels=c("prfi","etca","etco","prsi","przo","hyqu","apru"),labels=c("Opakapaka", "Ehu", "Onaga", "Kalekale", "Gindai", "Hapuupuu", "Lehi"))] %>%
		ggplot() +
		ylab("Empirical encounter rate") +
		xlab("Wind speed (kt)") +
		facet_wrap(~species_cd,scales="free_y") +
		ylim(0,NA) +
		geom_hline(yintercept=0) +
		geom_point(aes(x=wind_floor,y=encounter,fill=gear_type,size=N),shape=21) +
		geom_smooth(aes(x=wind_floor,y=encounter,color=gear_type,weight=N),se=FALSE) +
		theme_few(base_size=20) +
     	viridis::scale_color_viridis("Gear\ntype",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) + 
     	viridis::scale_fill_viridis("Gear\ntype",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
			
		ggsave(filename=paste0(data_flag,"wind_encounter_species.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	p = copy(wind_plot_dt) %>%
		.[,.(weight_kg=sum(weight_kg),speed_kt=mean(speed_kt)),by=.(design_sampling_unit,gear_type)] %>%
		.[,bin:=ifelse(weight_kg>0,1,0)] %>%
		.[,wind_floor:=floor(speed_kt)] %>%
		.[,.(encounter=mean(bin),.N),by=.(gear_type,wind_floor)] %>%
		ggplot() +
		ylab("Empirical encounter rate") +
		xlab("Wind speed (kt)") +
		ylim(0,NA) +
		geom_hline(yintercept=0) +
		geom_point(aes(x=wind_floor,y=encounter,fill=gear_type,size=N),shape=21) +
		geom_smooth(aes(x=wind_floor,y=encounter,color=gear_type,weight=N),se=FALSE) +
		theme_few(base_size=20) +
     	viridis::scale_color_viridis("Gear\ntype",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) + 
     	viridis::scale_fill_viridis("Gear\ntype",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)		
		ggsave(filename=paste0(data_flag,"wind_encounter.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	p = copy(wind_plot_dt) %>%
		.[,.(weight_kg=sum(weight_kg),speed_kt=mean(speed_kt)),by=.(design_sampling_unit,gear_type,island)] %>%
		.[,bin:=ifelse(weight_kg>0,1,0)] %>%
		.[,wind_floor:=floor(speed_kt)] %>%
		.[,.(encounter=mean(bin),.N),by=.(gear_type,wind_floor,island)] %>%
		ggplot() +
		facet_wrap(~island,scales="free_y") +
		ylab("Empirical encounter rate") +
		xlab("Wind speed (kt)") +
		ylim(0,NA) +
		geom_hline(yintercept=0) +
		geom_point(aes(x=wind_floor,y=encounter,fill=gear_type,size=N),shape=21) +
		geom_smooth(aes(x=wind_floor,y=encounter,color=gear_type,weight=N),se=FALSE) +
		theme_few(base_size=20) +
     	viridis::scale_color_viridis("Gear\ntype",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) + 
     	viridis::scale_fill_viridis("Gear\ntype",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)		
		ggsave(filename=paste0(data_flag,"wind_encounter_island.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	p = copy(wind_plot_dt) %>%
		.[,.(weight_kg=sum(weight_kg),speed_kt=mean(speed_kt)),by=.(design_sampling_unit,gear_type)] %>%
		.[weight_kg>0] %>%
		.[,wind_floor:=floor(speed_kt)] %>%
		.[,.(positive_catch=mean(weight_kg),.N),by=.(gear_type,wind_floor)] %>%
		ggplot() +
		facet_wrap(~gear_type,scales="free_y") +
		ylab("Positive catch (kg)") +
		xlab("Wind speed (kt)") +
		ylim(0,NA) +
		geom_hline(yintercept=0) +
		geom_point(aes(x=wind_floor,y=positive_catch,fill=gear_type,size=N),shape=21) +
		geom_smooth(aes(x=wind_floor,y=positive_catch,color=gear_type,weight=N),se=FALSE) +
		theme_few(base_size=20) +
     	viridis::scale_color_viridis("Gear\ntype",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) + 
     	viridis::scale_fill_viridis("Gear\ntype",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)		
		ggsave(filename=paste0(data_flag,"wind_poscatch.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)


	p = copy(wind_plot_dt) %>%
		.[,.(weight_kg=sum(weight_kg),speed_kt=mean(speed_kt)),by=.(design_sampling_unit,gear_type,island)] %>%
		.[weight_kg>0] %>%
		.[,wind_floor:=floor(speed_kt)] %>%
		.[,.(positive_catch=mean(weight_kg),.N),by=.(gear_type,wind_floor,island)] %>%
		ggplot() +
		facet_grid(gear_type~island,scales="free_y") +
		ylab("Positive catch (kg)") +
		xlab("Wind speed (kt)") +
		ylim(0,NA) +
		geom_hline(yintercept=0) +
		geom_point(aes(x=wind_floor,y=positive_catch,fill=gear_type,size=N),shape=21) +
		geom_smooth(aes(x=wind_floor,y=positive_catch,color=gear_type,weight=N),se=FALSE) +
		theme_few(base_size=20) +
     	viridis::scale_color_viridis("Gear\ntype",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) + 
     	viridis::scale_fill_viridis("Gear\ntype",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)		
		ggsave(filename=paste0(data_flag,"wind_poscatch_island.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	p = copy(sample_wind_dt) %>%
		merge(.,as.data.table(bfish_df)[,.(design_sampling_unit,gear_type)]) %>%
		.[,time_floor:=floor(time)] %>%
		.[,time_floor:=factor(time_floor)] %>%
		ggplot() +
		ylab("Wind speed (kt)") +
		xlab("Time of day (24h)") +
		ylim(0,NA) +
		geom_hline(yintercept=0) +
		geom_boxplot(aes(x=time_floor,y=speed_kt,fill=gear_type)) +
		theme_few(base_size=20) +
     	viridis::scale_fill_viridis("Gear\ntype",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)		
		ggsave(filename=paste0(data_flag,"wind_time.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	p = copy(sample_wind_dt) %>%
		merge(.,as.data.table(bfish_df)[,.(design_sampling_unit,island)]) %>%
		ggplot() +
		ylab("Latitude") +
		xlab("Longitude") +
		facet_wrap(~island,scales="free") +
		geom_point(aes(x=lon,y=lat,fill=speed_kt),shape=21,size=2) +
		theme_few(base_size=20) +
     	viridis::scale_fill_viridis("Wind\nspeed (kt)",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE)		
		ggsave(filename=paste0(data_flag,"wind_lonlat.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	p = wind_plot_dt %>%
		.[gear_type=="research_fishing"] %>%
		.[,.(design_sampling_unit,obs_wind,speed_kt,direction_d,island,time)] %>%
		unique(.) %>%
		ggplot() +
		xlab("Estimated observed wind speed (kt)") +
		ylab("PacIOOS wind speed (kt)") +
		ylim(0,NA) +
		geom_hline(yintercept=0) +
		geom_abline(intercept=0,slope=1,linetype="dashed",size=1.25,color="black") +
		geom_point(aes(x=jitter(obs_wind,factor=2),y=speed_kt,fill=time),shape=21,size=3) +
		theme_few(base_size=20) +
 		viridis::scale_fill_viridis("Time\nof day\n(24 h)",begin = 0.1,end = 0.8,direction = 1,option = "H")
		ggsave(filename=paste0(data_flag,"wind_rf_pacioos_v_obs.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	p = wind_plot_dt %>%
		.[gear_type=="research_fishing"] %>%
		.[,.(design_sampling_unit,obs_wind,speed_kt,direction_d,island,time)] %>%
		unique(.) %>%
		ggplot() +
		xlab("Estimated observed wind speed (kt)") +
		ylab("PacIOOS wind speed (kt)") +
		ylim(0,NA) +
		geom_hline(yintercept=0) +
		geom_abline(intercept=0,slope=1,linetype="dashed",size=1.25,color="black") +
		geom_point(aes(x=jitter(obs_wind,factor=2),y=speed_kt,fill=direction_d),shape=21,size=3) +
		theme_few(base_size=20) +
 		viridis::scale_fill_viridis("Direction",begin = 0.1,end = 0.8,discrete=TRUE,direction = 1,option = "H")
		ggsave(filename=paste0(data_flag,"wind_rf_pacioos_v_obs_direction.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	p = wind_plot_dt %>%
		.[gear_type=="research_fishing"] %>%
		.[,.(design_sampling_unit,obs_wind,speed_kt,direction_d,island,time,lon,lat)] %>%
		.[,diff:=speed_kt-obs_wind] %>%
		.[,lon_bin:=floor(lon/0.125)*0.125+0.06125] %>%
		.[,lat_bin:=floor(lat/0.125)*0.125+0.06125] %>%
		.[,bin_id := paste0(lon_bin,"_",lat_bin)] %>%
		.[,.(diff=mean(diff),lon=mean(lon_bin),lat=mean(lat_bin),.N),by=bin_id] %>%
		unique(.) %>%
		ggplot() +
		ggtitle("Wind speed (kt) residual") +
		xlab("Longitude") +
		ylab("Latitude") +
		geom_point(aes(x=lon,y=lat,fill=diff,size=N),shape=21) +
		theme_few(base_size=20) +
 		scale_fill_gradient2("Difference\nPacIOOS\nv.\n'observed'",low = "blue",mid = "white",high = "red")
		ggsave(filename=paste0(data_flag,"wind_rf_agg_wind_spatial_residual.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	p = wind_plot_dt %>%
		.[gear_type=="research_fishing"] %>%
		.[,.(design_sampling_unit,obs_wind,speed_kt,direction_d,island,time,lon,lat)] %>%
		.[,diff:=speed_kt-obs_wind] %>%
		unique(.) %>%
		ggplot() +
		ggtitle("Wind speed (kt) residual") +
		xlab("Time of day (24h)") +
		ylab("Difference PacIOOS v.'observed'") +
		geom_hline(yintercept=0) +
		geom_point(aes(x=time,y=diff),alpha=0.5) +
		geom_smooth(aes(x=time,y=diff)) +
		theme_few(base_size=20)
		ggsave(filename=paste0(data_flag,"wind_rf_time_residual.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	p = wind_plot_dt %>%
		.[gear_type=="research_fishing"] %>%
		.[,.(design_sampling_unit,obs_wind,speed_kt,direction_d,island,time,lon,lat)] %>%
		.[,diff:=speed_kt-obs_wind] %>%
		unique(.) %>%
		ggplot() +
		ggtitle("Wind speed (kt) residual") +
		xlab("Island group") +
		ylab("Difference PacIOOS v.'observed'") +
		geom_hline(yintercept=0) +
		geom_boxplot(aes(x=island,y=diff,fill=island)) +
		theme_few(base_size=20) +
 		viridis::scale_fill_viridis("Island group",begin = 0.1,end = 0.8,discrete=TRUE,direction = 1,option = "H")
		ggsave(filename=paste0(data_flag,"wind_rf_island_residual.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	p = wind_plot_dt %>%
		.[gear_type=="research_fishing"] %>%
		.[,.(design_sampling_unit,obs_wind,speed_kt,direction_d,island,time,lon,lat)] %>%
		.[,diff:=speed_kt-obs_wind] %>%
		unique(.) %>%
		ggplot() +
		ggtitle("Wind speed (kt) residual") +
		xlab("Direction") +
		ylab("Difference PacIOOS v.'observed'") +
		geom_hline(yintercept=0) +
		geom_boxplot(aes(x=direction_d,y=diff,fill=direction_d)) +
		theme_few(base_size=20) +
 		viridis::scale_fill_viridis("Wind direction",begin = 0.1,end = 0.8,discrete=TRUE,direction = 1,option = "H")
		ggsave(filename=paste0(data_flag,"wind_rf_direction_residual.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	p = wind_plot_dt %>%
		.[gear_type=="research_fishing"] %>%
		.[,.(design_sampling_unit,obs_wind,speed_kt,direction_d,island,time,lon,lat,year,month)] %>%
		.[,diff:=speed_kt-obs_wind] %>%
		unique(.) %>%
		ggplot() +
		ggtitle("Wind speed (kt) residual") +
		xlab("Month") +
		ylab("Difference PacIOOS v.'observed'") +
		geom_hline(yintercept=0) +
		geom_boxplot(aes(x=month,y=diff,fill=month)) +
		theme_few(base_size=20) +
 		viridis::scale_fill_viridis("Month",begin = 0.1,end = 0.8,discrete=TRUE,direction = 1,option = "H")
		ggsave(filename=paste0(data_flag,"wind_rf_month_residual.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	p = wind_plot_dt %>%
		.[gear_type=="research_fishing"] %>%
		.[,.(design_sampling_unit,obs_wind,speed_kt,direction_d,island,time,lon,lat,year,month)] %>%
		.[,diff:=speed_kt-obs_wind] %>%
		unique(.) %>%
		ggplot() +
		ggtitle("Wind speed (kt) residual") +
		xlab("Year") +
		ylab("Difference PacIOOS v.'observed'") +
		geom_hline(yintercept=0) +
		geom_boxplot(aes(x=year,y=diff,fill=year)) +
		theme_few(base_size=20) +
 		viridis::scale_fill_viridis("Year",begin = 0.1,end = 0.8,discrete=TRUE,direction = 1,option = "H")
		ggsave(filename=paste0(data_flag,"wind_rf_year_residual.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

