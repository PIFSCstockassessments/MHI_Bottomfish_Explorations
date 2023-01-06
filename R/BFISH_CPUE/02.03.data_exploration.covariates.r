
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
# 1) bring in data
    data_flag = 2021
    species = "mv"
    data_treatment = "01"
    lehi_filter = TRUE
	load(file=paste0(proj.dir,"Data/",data_flag,"_",data_treatment,".bfish_combined_long_dt.RData"))

	# subset to species
    if(species == "mv")
    { 
        target_species = c("prfi","etca","etco","prsi","przo","hyqu","apru")
    } else{
        target_species = species
    }
	bfish_df = bfish_combined_long_dt %>% .[species_cd %in% target_species] %>% as.data.frame(.)
	
	# remove sample with large lehi observation
    if(lehi_filter)
    {
        bfish_df =  subset(bfish_df,design_sampling_unit!="2021_Fall_32293")
    }
	

#_____________________________________________________________________________________________________________________________
# covariate exploration: categorical
	explore_dt = as.data.table(bfish_df)
	target_vec = c("island","strata","strata_2020","substrate","slope","depth_strata","depth_strata_2020","complexity","hardness","year","season","month","gear_type","platform")
	
	for(i in 1:length(target_vec))
	{
		p = copy(explore_dt) %>%
					setnames(.,target_vec[i],"target") %>%
					.[,weight_kg:=ifelse(weight_kg>0,1,0)] %>%
					.[,.(mean = mean(weight_kg),sd=sd(weight_kg),.N),by=.(species_cd,target)] %>%
					.[,y:=mean-1.96*sd] %>%
					.[,yend:=mean+1.96*sd] %>%
					.[,y:=ifelse(y<0,0,y)] %>%
					ggplot() +
					ggtitle("1st component") +
					ylab("Empricial effect") +
					xlab("Variable") +
					facet_wrap(~species_cd,scales="free_y") +
					geom_segment(aes(x=target,xend=target,y=y,yend=yend)) +
					geom_point(aes(x=target,y=mean,fill=target,size=N),shape=21) +
					theme_few(base_size=20) +
					theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
					viridis::scale_fill_viridis(target_vec[i],discrete=TRUE,begin = 0.1,end = 0.8,direction = 1,option = "H")
					ggsave(filename=paste0(data_flag,".",species,".",data_treatment,".",lehi_filter,".",target_vec[i],".1st.png"), plot = p, device = "png", path = plot_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)

				p = copy(explore_dt) %>%
					setnames(.,target_vec[i],"target") %>%
					.[weight_kg>0] %>%
					.[,.(mean = mean(weight_kg),sd=sd(weight_kg),.N),by=.(species_cd,target)] %>%
					.[,y:=mean-1.96*sd] %>%
					.[,yend:=mean+1.96*sd] %>%
					.[,y:=ifelse(y<0,0,y)] %>%
					ggplot() +
					ggtitle("2nd component") +
					ylab("Empricial effect") +
					xlab("Variable") +
					facet_wrap(~species_cd,scales="free_y") +
					geom_segment(aes(x=target,xend=target,y=y,yend=yend)) +
					geom_point(aes(x=target,y=mean,fill=target,size=N),shape=21) +
					theme_few(base_size=20) +
					theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
					viridis::scale_fill_viridis(target_vec[i],discrete=TRUE,begin = 0.1,end = 0.8,direction = 1,option = "H")
					ggsave(filename=paste0(data_flag,".",species,".",data_treatment,".",lehi_filter,".",target_vec[i],".2nd.png"), plot = p, device = "png", path = plot_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)	
	}

	target_vec = c("jd","lon","lat","depth","time","lunar_phase","obs_wind","obs_wave")
	
	for(i in 1:length(target_vec))
	{
		p = copy(explore_dt) %>%
					setnames(.,target_vec[i],"target") %>%
					.[,weight_kg:=ifelse(weight_kg>0,1,0)] %>%
					ggplot() +
					ggtitle("1st component") +
					ylab("Empricial effect") +
					xlab("Variable") +
					facet_wrap(~species_cd,scales="free_y") +
					geom_rug(aes(x=target),alpha=0.25) +
					geom_smooth(aes(x=target,y=weight_kg)) +
					theme_few(base_size=20) +
					theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
					viridis::scale_fill_viridis(target_vec[i],discrete=TRUE,begin = 0.1,end = 0.8,direction = 1,option = "H")
					ggsave(filename=paste0(data_flag,".",species,".",data_treatment,".",lehi_filter,".",target_vec[i],".1st.png"), plot = p, device = "png", path = plot_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)

				p = copy(explore_dt) %>%
					setnames(.,target_vec[i],"target") %>%
					.[weight_kg>0] %>%
					ggplot() +
					ggtitle("2nd component") +
					ylab("Empricial effect") +
					xlab("Variable") +
					facet_wrap(~species_cd,scales="free_y") +
					geom_point(aes(x=target,y=weight_kg),alpha=0.25) +
					geom_smooth(aes(x=target,y=weight_kg)) +
					theme_few(base_size=20) +
					theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
					viridis::scale_fill_viridis(target_vec[i],discrete=TRUE,begin = 0.1,end = 0.8,direction = 1,option = "H")
					ggsave(filename=paste0(data_flag,".",species,".",data_treatment,".",lehi_filter,".",target_vec[i],".2nd.png"), plot = p, device = "png", path = plot_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)	
	}		


					p = copy(explore_dt) %>%
					setnames(.,"depth","target") %>%
					.[,weight_kg:=ifelse(weight_kg>0,1,0)] %>%
					ggplot() +
					ggtitle("1st component") +
					ylab("Empricial effect") +
					xlab("Variable") +
					coord_cartesian(ylim=c(0,NA)) +
					# facet_wrap(~species_cd,scales="free_y") +
					geom_rug(aes(x=target),alpha=0.25) +
					geom_smooth(aes(x=target,y=weight_kg,color=gear_type)) +
					theme_few(base_size=20) +
					theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
					viridis::scale_color_viridis("gear_type",discrete=TRUE,begin = 0.1,end = 0.8,direction = 1,option = "H")
					ggsave(filename=paste0(data_flag,".",species,".",data_treatment,".",lehi_filter,".","depth_by_gear_agg",".1st.png"), plot = p, device = "png", path = plot_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)
					
					p = copy(explore_dt) %>%
					setnames(.,"depth","target") %>%
					.[,weight_kg:=ifelse(weight_kg>0,1,0)] %>%
					ggplot() +
					ggtitle("1st component") +
					ylab("Empricial effect") +
					xlab("Variable") +
					coord_cartesian(ylim=c(0,NA)) +
					facet_wrap(~species_cd,scales="free_y") +
					geom_rug(aes(x=target),alpha=0.25) +
					geom_smooth(aes(x=target,y=weight_kg,color=gear_type)) +
					theme_few(base_size=20) +
					theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
					viridis::scale_color_viridis("gear_type",discrete=TRUE,begin = 0.1,end = 0.8,direction = 1,option = "H")
					ggsave(filename=paste0(data_flag,".",species,".",data_treatment,".",lehi_filter,".","depth_by_gear",".1st.png"), plot = p, device = "png", path = plot_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)



