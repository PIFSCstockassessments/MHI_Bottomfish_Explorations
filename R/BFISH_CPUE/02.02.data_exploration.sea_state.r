

# Nicholas Ducharme-Barth
# 12/21/2022
# Data exploration of effect of observed wind, wave (e.g., sea state) on research fishing catches
# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.


#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
    library(ggplot2)
    library(ggthemes)


#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"
   	plot_dir = paste0(proj.dir,"Plot/BFISH_CPUE/")

#_____________________________________________________________________________________________________________________________
# define data_flag
	# data_flag = "" # only loads data up through 2020
	data_flag = "2021_" # includes data through 2021

#_____________________________________________________________________________________________________________________________
# load wind data
	load(file=paste0(proj.dir,"Data/sample_wind_dt.RData"))

#_____________________________________________________________________________________________________________________________
# load 02 data
	load(file=paste0(proj.dir,"Data/",data_flag,"02.bfish_combined_long_dt.RData"))

#_____________________________________________________________________________________________________________________________
# plot data exploration - research fishing
    mode = function(x)
    # grabs the first item in the vector in the case of a tie
    {
        tmp_class = class(x)[1]
        if(!(tmp_class %in% c("integer","numeric")))
        {
            x = as.character(x)
        }
        unique_x = unique(x)
        return(unique_x[which.max(tabulate(match(x, unique_x)))])
    }

    sea_state_dt = bfish_combined_long_dt %>%
                         .[,.(design_sampling_unit,gear_type,lon,lat,obs_wind,obs_wave,obs_current,weight_kg)] %>%
                         .[gear_type == "research_fishing"] %>%
                         .[obs_current=="G ",obs_current:="G"] %>%
                         .[obs_current=="O",obs_current:="0"] %>%
                         .[,.(lon=mean(lon),lat=mean(lat),obs_wind=mean(obs_wind),obs_wave=mean(obs_wave),obs_current=mode(obs_current),weight_kg=sum(weight_kg)),by=design_sampling_unit] %>%
                         unique(.) 
       
    p = sea_state_dt %>%
        .[,.(obs_wave,weight_kg)] %>%
        .[,bin_kg:=ifelse(weight_kg>0,1,0)] %>%
        .[,.(bin_kg=mean(bin_kg),.N),by=obs_wave] %>%
        merge(.,sea_state_dt[weight_kg>0,.(weight_kg=mean(weight_kg)),by=obs_wave],by="obs_wave",all.x=TRUE) %>%
        ggplot() +
        ylim(0,NA) +
        geom_hline(yintercept = 0) +
        xlab("Observed wave (ft)") +
        ylab("Encounter rate") +  
        geom_point(aes(x=obs_wave,y=bin_kg,size=N,fill=weight_kg),shape=21) +
        theme_few(base_size=20) +
     	viridis::scale_fill_viridis("Positive catch",begin = 0.1,end = 0.8,direction = 1,option = "H")
		ggsave(filename=paste0("obs_wave_encounter.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
    
    p = sea_state_dt %>%
        .[,.(obs_wind,weight_kg)] %>%
        .[,bin_kg:=ifelse(weight_kg>0,1,0)] %>%
        .[,.(bin_kg=mean(bin_kg),.N),by=obs_wind] %>%
        merge(.,sea_state_dt[weight_kg>0,.(weight_kg=mean(weight_kg)),by=obs_wind],by="obs_wind",all.x=TRUE) %>%
        ggplot() +
        ylim(0,NA) +
        geom_hline(yintercept = 0) +
        xlab("Observed wind (kt)") +
        ylab("Encounter rate") + 
        geom_point(aes(x=obs_wind,y=bin_kg,size=N,fill=weight_kg),shape=21) +
        theme_few(base_size=20) +
     	viridis::scale_fill_viridis("Positive catch",begin = 0.1,end = 0.8,direction = 1,option = "H")
		ggsave(filename=paste0("obs_wind_encounter.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

    p = sea_state_dt %>%
        .[,.(obs_current,weight_kg)] %>%
        .[,bin_kg:=ifelse(weight_kg>0,1,0)] %>%
        .[,.(bin_kg=mean(bin_kg),.N),by=obs_current] %>%
        merge(.,sea_state_dt[weight_kg>0,.(weight_kg=mean(weight_kg)),by=obs_current],by="obs_current",all.x=TRUE) %>%
        ggplot() +
        ylim(0,NA) +
        geom_hline(yintercept = 0) +
        xlab("Observed current") +
        ylab("Encounter rate") + 
        geom_point(aes(x=obs_current,y=bin_kg,size=N,fill=weight_kg),shape=21) +
        theme_few(base_size=20) +
     	viridis::scale_fill_viridis("Positive catch",begin = 0.1,end = 0.8,direction = 1,option = "H")
		ggsave(filename=paste0("obs_current_encounter.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

  tmp_dt = sea_state_dt %>%
            .[weight_kg>0] %>%
            .[,lon_bin:=floor(lon/0.125)*0.125+0.06125] %>%
		    .[,lat_bin:=floor(lat/0.125)*0.125+0.06125] %>%
            .[,bin_id := paste0(lon_bin,"_",lat_bin)] %>%
            .[,.(weight_kg=mean(weight_kg)),by=bin_id]

    p = sea_state_dt %>%
        .[,.(lon,lat,obs_wave,obs_wind,weight_kg)] %>%
        .[,bin_kg:=ifelse(weight_kg>0,1,0)] %>%
        .[,lon_bin:=floor(lon/0.125)*0.125+0.06125] %>%
		.[,lat_bin:=floor(lat/0.125)*0.125+0.06125] %>%
        .[,bin_id := paste0(lon_bin,"_",lat_bin)] %>%
        .[,.(lon=mean(lon_bin),lat=mean(lat_bin),bin_kg=mean(bin_kg),obs_wave=mean(obs_wave),obs_wind=mean(obs_wind),.N),by=bin_id] %>%
        merge(.,tmp_dt,by="bin_id",all.x=TRUE) %>%
        melt(.,id.vars=c("bin_id","lon","lat","N")) %>%
        .[,variable:=factor(variable,levels=c("bin_kg","weight_kg","obs_wind","obs_wave"),labels=c("Encounter rate","Positive catch","Obs. wind","Obs. wave"))] %>%
        .[,value:=scale(value),by=variable] %>%
        ggplot() +
        facet_wrap(~variable) +
        xlab("Longitude") +
        ylab("Latitude") +  
        geom_point(aes(x=lon,y=lat,size=N,fill=value),shape=21) +
        theme_few(base_size=20) +
   		scale_fill_gradient2("Relative value",low = "blue",mid = "white",high = "red")
		ggsave(filename=paste0("rf_spatial_encounter_catch_wind_wave.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

 p = sea_state_dt %>%
        .[,.(lon,lat,obs_current)] %>%
        .[,lon_bin:=floor(lon/0.125)*0.125+0.06125] %>%
		.[,lat_bin:=floor(lat/0.125)*0.125+0.06125] %>%
        .[,bin_id := paste0(lon_bin,"_",lat_bin)] %>%
        .[,.(lon=mean(lon_bin),lat=mean(lat_bin),obs_current=mode(obs_current),.N),by=bin_id] %>%
        ggplot() +
        xlab("Longitude") +
        ylab("Latitude") +  
        geom_point(aes(x=lon,y=lat,size=N,fill=obs_current),shape=21) +
        theme_few(base_size=20) +
     	viridis::scale_fill_viridis("Obs.\ncurrent",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
		ggsave(filename=paste0("rf_spatial_current.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)		

#_____________________________________________________________________________________________________________________________
# plot data exploration - camera
    camera_dt = bfish_combined_long_dt %>%
                         .[,.(design_sampling_unit,gear_type,lon,lat,weight_kg)] %>%
                         .[gear_type == "camera"] %>%
                         .[,.(lon=mean(lon),lat=mean(lat),weight_kg=sum(weight_kg)),by=design_sampling_unit] %>%
                         unique(.)

    tmp_dt = camera_dt %>%
            .[weight_kg>0] %>%
            .[,lon_bin:=floor(lon/0.125)*0.125+0.06125] %>%
		    .[,lat_bin:=floor(lat/0.125)*0.125+0.06125] %>%
            .[,bin_id := paste0(lon_bin,"_",lat_bin)] %>%
            .[,.(weight_kg=mean(weight_kg)),by=bin_id]

    p = camera_dt %>%
        .[,.(lon,lat,weight_kg)] %>%
        .[,bin_kg:=ifelse(weight_kg>0,1,0)] %>%
        .[,lon_bin:=floor(lon/0.125)*0.125+0.06125] %>%
		.[,lat_bin:=floor(lat/0.125)*0.125+0.06125] %>%
        .[,bin_id := paste0(lon_bin,"_",lat_bin)] %>%
        .[,.(lon=mean(lon_bin),lat=mean(lat_bin),bin_kg=mean(bin_kg),.N),by=bin_id] %>%
        merge(.,tmp_dt,by="bin_id",all.x=TRUE) %>%
        melt(.,id.vars=c("bin_id","lon","lat","N")) %>%
        .[,variable:=factor(variable,levels=c("bin_kg","weight_kg"),labels=c("Encounter rate","Positive catch"))] %>%
        .[,value:=scale(value),by=variable] %>%
        ggplot() +
        facet_wrap(~variable) +
        xlab("Longitude") +
        ylab("Latitude") +  
        geom_point(aes(x=lon,y=lat,size=N,fill=value),shape=21) +
        theme_few(base_size=20) +
   		scale_fill_gradient2("Relative value",low = "blue",mid = "white",high = "red")
		ggsave(filename=paste0("camera_spatial_encounter_catch.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)