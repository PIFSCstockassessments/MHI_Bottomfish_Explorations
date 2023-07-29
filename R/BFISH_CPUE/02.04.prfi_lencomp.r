

# Nicholas Ducharme-Barth
# 2023/06/01
# calc length composition effective sample size for opakapaka
# plot length comp at depth by gear
# Copyright (c) 2023 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(ggplot2)

#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"
	plot_dir = paste0(proj.dir,"Plot/BFISH_CPUE/")
#_____________________________________________________________________________________________________________________________
# define data_flag
	# data_flag = "" # only loads data up through 2020
	data_flag = "2022" # includes data through 2022

#_____________________________________________________________________________________________________________________________
# plot length comps by gear and depth
	load(file=paste0(proj.dir,"Data/",data_flag,"_01.bfish_combined_long_dt.all_lengths.RData"))
	bfish_dt = bfish_combined_long_dt
	bfish_cam_lengths = fread(paste0(proj.dir,"Data/",data_flag,"_CAM_LENGTHS.csv")) %>%
			      .[,LENGTH_CM:=round(MEAN_MM/10)] %>%
				  .[,.(DROP_CD,SPECIES_CD,LENGTH_CM)] %>%
				  .[SPECIES_CD %in% "PRFI" & DROP_CD %in% bfish_combined_long_dt$model_sampling_unit] %>%
				  setnames(.,c("DROP_CD","SPECIES_CD","LENGTH_CM"),c("model_sampling_unit","species_cd","length_cm")) %>%
				  merge(.,unique(bfish_combined_long_dt[,.(model_sampling_unit,depth)]),by=c("model_sampling_unit"),all.x=TRUE) %>%
				  .[,gear:="camera"]

	bfish_rf_lengths = fread(paste0(proj.dir,"Data/",data_flag,"_CRF_CATCH.csv")) %>%
			  .[SPECIES_CD == "PRFI"& SAMPLE_ID %in% bfish_combined_long_dt$model_sampling_unit] %>%
			  .[,LENGTH_CM:=round(as.numeric(LENGTH_CM))] %>%
			  .[,.(SAMPLE_ID,SPECIES_CD,LENGTH_CM)] %>%
			  na.omit(.) %>%
			  setnames(.,c("SAMPLE_ID","SPECIES_CD","LENGTH_CM"),c("model_sampling_unit","species_cd","length_cm")) %>%
			  merge(.,unique(bfish_combined_long_dt[,.(model_sampling_unit,depth)]),by=c("model_sampling_unit"),all.x=TRUE) %>%
				.[,gear:="research fishing"]
	
	bfish_lengths = rbind(bfish_cam_lengths,bfish_rf_lengths) %>%
					.[,depth_bin_50:=floor(depth/50)*50] %>%
					.[,depth_bin_100:=floor(depth/100)*100] %>%
					.[,length_bin_10:=floor(length_cm/10)*10]
    rf_sample = unique(bfish_combined_long_dt[gear_type=="research_fishing"]$model_sampling_unit)
    rm(list=c("bfish_combined_long_dt"))                

#_____________________________________________________________________________________________________________________________
# loop over different options and get opakapaka sample sizes (assume 1 for each model sampling event)
	lencomp_dt.list = sample_size_dt.list = as.list(rep(NA,10))
    
	iter=1
	for(all_lengths in c(TRUE,FALSE))
	{
		for(d_t in 1:5)
		{
			if(all_lengths)
			{
				load(file=paste0(proj.dir,"Data/",data_flag,"_0",d_t,".bfish_combined_long_dt.all_lengths.RData"))
			} else {
				load(file=paste0(proj.dir,"Data/",data_flag,"_0",d_t,".bfish_combined_long_dt.RData"))
			}
			
			sample_size_dt.list[[iter]] = bfish_combined_long_dt %>%
			.[species_cd=="prfi"&weight_kg>0] %>%
			.[,.(input_sample_size=.N),by=.(year,gear_type)] %>%
			.[order(year,gear_type)] %>%
			.[,data_treatment:=d_t] %>%
			.[,all_lengths:=as.character(all_lengths)] %>%
			.[,.(data_treatment,all_lengths,year,gear_type,input_sample_size)]

            # get appropriate camera and research fishing samples
            if(all_lengths)
            {
                drop_cd = unique(fread(paste0(proj.dir,"Data/",data_flag,"_0",d_t,".camera_dt.all_lengths.drop_cd.csv"))$DROP_CD)
            } else {
                drop_cd = unique(fread(paste0(proj.dir,"Data/",data_flag,"_0",d_t,".camera_dt.all_lengths.drop_cd.csv"))$DROP_CD)
            }
            model_samp_vec = unique(c(drop_cd,rf_sample))

            lencomp_dt.list[[iter]] = merge(unique(bfish_dt[model_sampling_unit %in% model_samp_vec,.(model_sampling_unit,year)]),bfish_lengths[,.( model_sampling_unit,gear,length_cm,depth)],by="model_sampling_unit",all.x=TRUE) %>%
                                      .[,data_treatment:=d_t] %>%
			                          .[,all_lengths:=as.character(all_lengths)] %>%
                                      .[,.(data_treatment,all_lengths,year,gear,length_cm,depth)]
            if(!all_lengths)
            {
                lencomp_dt.list[[iter]] = lencomp_dt.list[[iter]][length_cm>=29]
            }

		rm(list=c("bfish_combined_long_dt","drop_cd"))
		iter = iter+1
		}
	}

	sample_size_dt = rbindlist(sample_size_dt.list)
    lencomp_dt = na.omit(rbindlist(lencomp_dt.list))
	save(sample_size_dt,file=paste0(proj.dir,"Data/",data_flag,".sample_size_dt.RData"))
	save(lencomp_dt,file=paste0(proj.dir,"Data/",data_flag,".lencomp_dt.RData"))

#_____________________________________________________________________________________________________________________________
# plot length comps by gear and depth

	p = lencomp_dt[gear=="camera"&all_lengths==TRUE] %>%
	ggplot() +
	xlab("Length (cm)") +
	ylab("Density") +
	geom_hline(yintercept=0) +
	geom_vline(xintercept=29,linetype="dashed") +
	xlim(0,NA) +
	geom_density(aes(x=length_cm,y=after_stat(ndensity),color=as.character(data_treatment),group=as.character(data_treatment)),linewidth=2) +
	# geom_density(aes(x=length_cm,linetype=gear,color=depth_bin,group=paste0(gear,".",depth_bin))) +
	viridis::scale_color_viridis("Data treatment",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
	theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
				panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
				panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
				strip.background =element_rect(fill="white"),
				legend.key = element_rect(fill = "white"))
	ggsave(filename=paste0(data_flag,"_lf.dens.camera.datatreatment.png"), plot = p, device = "png", path = plot_dir,
						scale = 0.75, width = 12, height = 8, units = c("in"),
						dpi = 300, limitsize = TRUE)
			  
	p = bfish_lengths %>%
	ggplot() +
	xlab("Length (cm)") +
	ylab("Density") +
	geom_hline(yintercept=0) +
	geom_vline(xintercept=29,linetype="dashed") +
	xlim(0,NA) +
	geom_density(aes(x=length_cm,y=after_stat(ndensity),color=depth_bin_50,group=depth_bin_50),linewidth=2) +
	# geom_density(aes(x=length_cm,linetype=gear,color=depth_bin,group=paste0(gear,".",depth_bin))) +
	viridis::scale_color_viridis("Depth (m)",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
	theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
				panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
				panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
				strip.background =element_rect(fill="white"),
				legend.key = element_rect(fill = "white"))
	ggsave(filename=paste0(data_flag,"_lf.dens.depth50.png"), plot = p, device = "png", path = plot_dir,
						scale = 0.75, width = 12, height = 8, units = c("in"),
						dpi = 300, limitsize = TRUE)

	p = bfish_lengths %>%
	ggplot() +
	xlab("Length (cm)") +
	ylab("Density") +
	geom_hline(yintercept=0) +
	geom_vline(xintercept=29,linetype="dashed") +
	xlim(0,NA) +
	geom_density(aes(x=length_cm,y=after_stat(ndensity),color=depth_bin_100,group=depth_bin_100),linewidth=2) +
	# geom_density(aes(x=length_cm,linetype=gear,color=depth_bin,group=paste0(gear,".",depth_bin))) +
	viridis::scale_color_viridis("Depth (m)",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
	theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
				panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
				panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
				strip.background =element_rect(fill="white"),
				legend.key = element_rect(fill = "white"))
		ggsave(filename=paste0(data_flag,"_lf.dens.depth100.png"), plot = p, device = "png", path = plot_dir,
						scale = 0.75, width = 12, height = 8, units = c("in"),
						dpi = 300, limitsize = TRUE)

	p = bfish_lengths %>%
	ggplot() +
	xlab("Length (cm)") +
	ylab("Density") +
	geom_hline(yintercept=0) +
	geom_vline(xintercept=29,linetype="dashed") +
	xlim(0,NA) +
	geom_density(aes(x=length_cm,y=after_stat(ndensity),linetype=gear,group=gear),linewidth=2) +
	# geom_density(aes(x=length_cm,linetype=gear,color=depth_bin,group=paste0(gear,".",depth_bin))) +
	# viridis::scale_color_viridis("Depth",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
	theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
				panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
				panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
				strip.background =element_rect(fill="white"),
				legend.key = element_rect(fill = "white"))
	ggsave(filename=paste0(data_flag,"_lf.dens.gear.png"), plot = p, device = "png", path = plot_dir,
						scale = 0.75, width = 12, height = 8, units = c("in"),
						dpi = 300, limitsize = TRUE)

	p = bfish_lengths %>%
	ggplot() +
	xlab("Length (cm)") +
	ylab("Density") +
	geom_hline(yintercept=0) +
	geom_vline(xintercept=29,linetype="dashed") +
	xlim(0,NA) +
	facet_wrap(~gear,nrow=2) +
	geom_density(aes(x=length_cm,y=after_stat(ndensity),linetype=gear,color=depth_bin_50,group=paste0(gear,".",depth_bin_50)),linewidth=2) +
	viridis::scale_color_viridis("Depth (m)",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
	theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
				panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
				panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
				strip.background =element_rect(fill="white"),
				legend.key = element_rect(fill = "white"))
	ggsave(filename=paste0(data_flag,"_lf.dens.gear_depth50.png"), plot = p, device = "png", path = plot_dir,
						scale = 0.75, width = 12, height = 8, units = c("in"),
						dpi = 300, limitsize = TRUE)

	p = bfish_lengths %>%
	ggplot() +
	xlab("Length (cm)") +
	ylab("Number of samples") +
	geom_hline(yintercept=0) +
	geom_vline(xintercept=29,linetype="dashed") +
	xlim(0,NA) +
	facet_wrap(~gear,nrow=2) +
	geom_histogram(aes(x=length_cm,fill=depth_bin_50,group=paste0(gear,".",depth_bin_50)),binwidth=2) +
	viridis::scale_fill_viridis("Depth (m)",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
	theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
				panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
				panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
				strip.background =element_rect(fill="white"),
				legend.key = element_rect(fill = "white"))
	ggsave(filename=paste0(data_flag,"_lf.hist.gear_depth50.png"), plot = p, device = "png", path = plot_dir,
						scale = 0.75, width = 12, height = 8, units = c("in"),
						dpi = 300, limitsize = TRUE)

	p = bfish_lengths %>%
	.[order(length_cm,depth)] %>%
	ggplot() +
	xlab("Depth (m)") +
	ylab("Number of samples") +
	geom_hline(yintercept=0) +
	xlim(0,NA) +
	facet_wrap(~gear,nrow=2) +
	geom_histogram(aes(x=depth,fill=length_bin_10,group=paste0(gear,".",length_bin_10)),binwidth=10) +
	viridis::scale_fill_viridis("Length (cm)",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=FALSE) +
	theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
				panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
				panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
				strip.background =element_rect(fill="white"),
				legend.key = element_rect(fill = "white"))
	ggsave(filename=paste0(data_flag,"_lf.hist.gear_len10.png"), plot = p, device = "png", path = plot_dir,
						scale = 0.75, width = 12, height = 8, units = c("in"),
						dpi = 300, limitsize = TRUE)


