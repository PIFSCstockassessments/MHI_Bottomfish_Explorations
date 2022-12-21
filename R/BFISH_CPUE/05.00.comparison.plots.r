
#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(ggplot2)
	library(ggthemes)


#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"
	dir_plot = paste0(proj.dir,"VAST/model_runs/comparisons/")
	dir.create(dir_plot,recursive = TRUE)

#_____________________________________________________________________________________________________________________________
# bring in indices
	design_dt = fread(paste0(proj.dir,"Data/2021_design_based_estimates.csv")) %>%
				.[,.(Model,Category,Time,Estimate,CV)]

	data_flag = "2021_" # includes data through 2021
	load(file=paste0(proj.dir,"Data/",data_flag,"bfish_combined_long_dt.RData"))
	# collapse weight_kg across bait_type
		sample_weight_dt = bfish_combined_long_dt[,.(weight_kg=sum(weight_kg)),by=.(sample_id,species_cd)]
		unique_dt = copy(bfish_combined_long_dt) %>% .[,bait_type:=NULL] %>% .[,weight_kg:=NULL] %>% unique(.)
		nominal_species_dt = merge(unique_dt,sample_weight_dt,by=c("sample_id","species_cd")) %>%
		.[sample_id!="20210831_183445"] %>%
		.[,.(Estimate=mean(weight_kg)),by=.(year,species_cd)] %>%
		.[,Model:="nominal"] %>%
		.[,CV:=NA] %>%
		setnames(.,c("year","species_cd"),c("Time","Category")) %>%
		.[,.(Model,Category,Time,Estimate,CV)]

		nominal_all_dt = merge(unique_dt,sample_weight_dt,by=c("sample_id","species_cd")) %>%
		.[sample_id!="20210831_183445"] %>%
		.[,.(Estimate=mean(weight_kg)),by=.(year)] %>%
		.[,Model:="nominal"] %>%
		.[,CV:=NA] %>%
		.[,Category:="all"] %>%
		setnames(.,c("year"),c("Time")) %>%
		.[,.(Model,Category,Time,Estimate,CV)]

		nominal_dt = rbind(nominal_all_dt,nominal_species_dt) %>%
					.[,Estimate:=Estimate/mean(Estimate),by=Category] %>%
					.[,Time:=as.numeric(as.character(Time))]

	

	mv_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_mv7_vanilla/Index.csv")) %>%
			.[,Category:=factor(Category,levels=paste0("Category_",1:7),labels=c("prfi","etca","etco","prsi","przo","hyqu","apru"))] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="mv"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]

	mvpl_dt = fread(paste0(proj.dir,"VAST/model_runs/pl-dgamma_mv7_vanilla/Index.csv")) %>%
			.[,Category:=factor(Category,levels=paste0("Category_",1:7),labels=c("prfi","etca","etco","prsi","przo","hyqu","apru"))] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="mvpl"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]

	prfi_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_prfi/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="prfi")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]
	etca_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_etca/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="etca")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]
	etco_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_etco/simple/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="etco")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]
	prsi_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_prsi/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="prsi")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]
	przo_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_przo/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="przo")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]	
	hyqu_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_hyqu/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="hyqu")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]
	apru_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_apru/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="apru")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]
	apru_nof_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_apru_nofilter/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="apru")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single - no filter"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]

	mv_agg_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_mv7_vanilla/agg/Index.csv")) %>%
			.[Category=="Category_7"] %>%
			.[,Category:="all"] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="mv"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]

	mvpl_agg_dt = fread(paste0(proj.dir,"VAST/model_runs/pl-dgamma_mv7_vanilla/agg/Index.csv")) %>%
			.[Category=="Category_7"] %>%
			.[,Category:="all"] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="mvpl"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]

index_dt = rbind(mv_dt,prfi_dt,etca_dt,etco_dt,prsi_dt,przo_dt,hyqu_dt,apru_dt,apru_nof_dt)
index_2_dt = rbind(mv_dt,mvpl_dt,prfi_dt,etca_dt,etco_dt,prsi_dt,przo_dt,hyqu_dt,apru_dt,apru_nof_dt)


p1 = index_dt %>%
	ggplot() +
	ylim(0,NA) +
	facet_wrap(~Category) +
	geom_hline(yintercept=0) +
	geom_path(aes(x=Time,y=Estimate/1000000,group=Model,color=Model)) +
	ylab("Estimate (million kg)") +
	theme_few(base_size = 20)
    ggsave(filename=paste0("p1.png"), plot = p1, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)


p2 = index_dt %>%
	ggplot() +
	ylim(0,NA) +
	facet_wrap(~Category,scales="free_y") +
	geom_hline(yintercept=0) +
	geom_path(aes(x=Time,y=Estimate/1000000,group=Model,color=Model)) +
	ylab("Estimate (million kg)") +
	theme_few(base_size = 20)
    ggsave(filename=paste0("p2.png"), plot = p2, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

p3 =  index_dt %>%
	 .[Model!="single - no filter",.(Estimate=sum(Estimate)),by=.(Model,Time)] %>%
	 ggplot() +
	 ylim(0,NA) +
	 geom_hline(yintercept=0) +
	 geom_path(aes(x=Time,y=Estimate/1000000,group=Model,color=Model)) +
	ylab("Estimate (million kg)") +
	 theme_few(base_size = 20)
    ggsave(filename=paste0("p3.png"), plot = p3, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)


p4 =  rbind(mv_agg_dt,mvpl_agg_dt) %>% 
	  .[,l95:=exp(log(Estimate/1000000)-2*CV)] %>%
	  .[,u95:=exp(log(Estimate/1000000)+2*CV)] %>%
	 ggplot() +
	 ylim(0,NA) +
	 geom_hline(yintercept=0) +
	 geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	 geom_path(aes(x=Time,y=Estimate/1000000,group=Model,color=Model)) +
	ylab("Estimate (million kg)") +
	 theme_few(base_size = 20)
    ggsave(filename=paste0("p4.png"), plot = p4, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)


p5 = rbind(mv_dt,mvpl_dt) %>% 
	  .[,l95:=exp(log(Estimate/1000000)-2*CV)] %>%
	  .[,u95:=exp(log(Estimate/1000000)+2*CV)] %>%
	ggplot() +
	ylim(0,NA) +
	facet_wrap(~Category,scales="free_y") +
	geom_hline(yintercept=0) +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	geom_path(aes(x=Time,y=Estimate/1000000,group=Model,color=Model)) +
	ylab("Estimate (million kg)") +
	theme_few(base_size = 20)
    ggsave(filename=paste0("p5.png"), plot = p5, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

p6 =  index_2_dt %>%
	 .[Model=="single",.(Estimate=sum(Estimate)),by=.(Model,Time)] %>%
	 .[,Category:="all"] %>%
	 .[,CV:=NA] %>%
	 .[,.(Model,Category,Time,Estimate,CV)] %>%
	 rbind(.,mv_agg_dt,mvpl_agg_dt) %>%
	 .[,l95:=exp(log(Estimate/1000000)-2*CV)] %>%
	 .[,u95:=exp(log(Estimate/1000000)+2*CV)] %>%
	 ggplot() +
	 ylim(0,NA) +
	 geom_hline(yintercept=0) +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	 geom_path(aes(x=Time,y=Estimate/1000000,group=Model,color=Model)) +
	ylab("Estimate (million kg)") +
	 theme_few(base_size = 20)
    ggsave(filename=paste0("p6.png"), plot = p6, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

p7 = index_2_dt %>%
	.[Model=="single"] %>%
	 .[,.(Model,Category,Time,Estimate,CV)] %>%
	 rbind(.,mv_dt,mvpl_dt) %>%
	 .[,l95:=exp(log(Estimate/1000000)-2*CV)] %>%
	 .[,u95:=exp(log(Estimate/1000000)+2*CV)] %>%
	ggplot() +
	ylim(0,NA) +
	facet_wrap(~Category,scales="free_y") +
	geom_hline(yintercept=0) +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	geom_path(aes(x=Time,y=Estimate/1000000,group=Model,color=Model)) +
	ylab("Estimate (million kg)") +
	theme_few(base_size = 20)
    ggsave(filename=paste0("p7.png"), plot = p7, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

p8 =  index_2_dt %>%
	 .[Model=="single",.(Estimate=sum(Estimate)),by=.(Model,Time)] %>%
	 .[,Category:="all"] %>%
	 .[,CV:=NA] %>%
	 .[,.(Model,Category,Time,Estimate,CV)] %>%
	 rbind(.,mv_agg_dt,mvpl_agg_dt) %>%
	 .[,Scaled_Estimate:=(Estimate/mean(Estimate)),by=Model] %>%
	 .[,l95:=(exp(log(Scaled_Estimate)-2*CV))] %>%
	 .[,u95:=(exp(log(Scaled_Estimate)+2*CV))] %>%
	 ggplot() +
	 ylim(0,NA) +
	 geom_hline(yintercept=0) +
	 geom_hline(yintercept=1,linetype="dashed") +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	 geom_path(aes(x=Time,y=Scaled_Estimate,group=Model,color=Model)) +
	geom_point(data=nominal_dt[Category=="all"],aes(x=Time,y=Estimate,group=Model),size=2) +
	ylab("Estimate (relative)") +
	 theme_few(base_size = 20)
    ggsave(filename=paste0("p8.png"), plot = p8, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

p9 = index_2_dt %>%
	.[Model=="single"] %>%
	 .[,.(Model,Category,Time,Estimate,CV)] %>%
	 rbind(.,mv_dt,mvpl_dt) %>%
	 	 .[,Scaled_Estimate:=(Estimate/mean(Estimate)),by=.(Model,Category)] %>%
	 .[,l95:=exp(log(Scaled_Estimate)-2*CV)] %>%
	 .[,u95:=exp(log(Scaled_Estimate)+2*CV)] %>%
	ggplot() +
	ylim(0,NA) +
	facet_wrap(~Category,scales="free_y") +
	geom_hline(yintercept=0) +
	 geom_hline(yintercept=1,linetype="dashed") +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	geom_path(aes(x=Time,y=Scaled_Estimate,group=Model,color=Model)) +
	geom_point(data=nominal_dt[Category!="all"],aes(x=Time,y=Estimate,group=Model),size=2) +
	ylab("Estimate (relative)") +
	theme_few(base_size = 20)
    ggsave(filename=paste0("p9.png"), plot = p9, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

p10 =  index_2_dt %>%
	 .[Model=="single",.(Estimate=sum(Estimate)),by=.(Model,Time)] %>%
	 .[,Category:="all"] %>%
	 .[,CV:=NA] %>%
	 .[,.(Model,Category,Time,Estimate,CV)] %>%
	 rbind(.,mv_agg_dt,mvpl_agg_dt) %>%
	 rbind(.,design_dt[Category=="all"]) %>%
	 .[,l95:=exp(log(Estimate/1000000)-2*CV)] %>%
	 .[,u95:=exp(log(Estimate/1000000)+2*CV)] %>%
	 ggplot() +
	 ylim(0,NA) +
	 geom_hline(yintercept=0) +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	 geom_path(aes(x=Time,y=Estimate/1000000,group=Model,color=Model)) +
	ylab("Estimate (million kg)") +
	 theme_few(base_size = 20)
    ggsave(filename=paste0("p10.png"), plot = p10, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

p11 = index_2_dt %>%
	.[Model=="single"] %>%
	 .[,.(Model,Category,Time,Estimate,CV)] %>%
	 rbind(.,mv_dt,mvpl_dt) %>%
	 rbind(.,design_dt[Category!="all"]) %>%
	 .[,l95:=exp(log(Estimate/1000000)-2*CV)] %>%
	 .[,u95:=exp(log(Estimate/1000000)+2*CV)] %>%
	ggplot() +
	ylim(0,NA) +
	facet_wrap(~Category,scales="free_y") +
	geom_hline(yintercept=0) +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	geom_path(aes(x=Time,y=Estimate/1000000,group=Model,color=Model)) +
	ylab("Estimate (million kg)") +
	theme_few(base_size = 20)
    ggsave(filename=paste0("p11.png"), plot = p11, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

p12 =  index_2_dt %>%
	 .[Model=="single",.(Estimate=sum(Estimate)),by=.(Model,Time)] %>%
	 .[,Category:="all"] %>%
	 .[,CV:=NA] %>%
	 .[,.(Model,Category,Time,Estimate,CV)] %>%
	 rbind(.,mv_agg_dt,mvpl_agg_dt) %>%
	 rbind(.,design_dt[Category=="all"]) %>%
	 .[,Scaled_Estimate:=(Estimate/mean(Estimate)),by=Model] %>%
	 .[,l95:=(exp(log(Scaled_Estimate)-2*CV))] %>%
	 .[,u95:=(exp(log(Scaled_Estimate)+2*CV))] %>%
	 ggplot() +
	 ylim(0,NA) +
	 geom_hline(yintercept=0) +
	 geom_hline(yintercept=1,linetype="dashed") +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	 geom_path(aes(x=Time,y=Scaled_Estimate,group=Model,color=Model)) +
	ylab("Estimate (relative)") +
	 theme_few(base_size = 20)
    ggsave(filename=paste0("p12.png"), plot = p12, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

p13 = index_2_dt %>%
	.[Model=="single"] %>%
	 .[,.(Model,Category,Time,Estimate,CV)] %>%
	 rbind(.,mv_dt,mvpl_dt) %>%
	 rbind(.,design_dt[Category!="all"]) %>%
	 	 .[,Scaled_Estimate:=(Estimate/mean(Estimate)),by=.(Model,Category)] %>%
	 .[,l95:=exp(log(Scaled_Estimate)-2*CV)] %>%
	 .[,u95:=exp(log(Scaled_Estimate)+2*CV)] %>%
	ggplot() +
	ylim(0,NA) +
	facet_wrap(~Category,scales="free_y") +
	geom_hline(yintercept=0) +
	 geom_hline(yintercept=1,linetype="dashed") +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	geom_path(aes(x=Time,y=Scaled_Estimate,group=Model,color=Model)) +
	ylab("Estimate (relative)") +
	theme_few(base_size = 20)
    ggsave(filename=paste0("p13.png"), plot = p13, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

p14 =  index_2_dt %>%
	 .[Model=="single",.(Estimate=sum(Estimate)),by=.(Model,Time)] %>%
	 .[,Category:="all"] %>%
	 .[,CV:=NA] %>%
	 .[,.(Model,Category,Time,Estimate,CV)] %>%
	 rbind(.,mv_agg_dt,mvpl_agg_dt) %>%
	 rbind(.,design_dt[Category=="all"]) %>%
	 .[,Scaled_Estimate:=(Estimate/mean(Estimate)),by=Model] %>%
	 .[,l95:=(exp(log(Scaled_Estimate)-2*CV))] %>%
	 .[,u95:=(exp(log(Scaled_Estimate)+2*CV))] %>%
	 ggplot() +
	 ylim(0,NA) +
	 geom_hline(yintercept=0) +
	 geom_hline(yintercept=1,linetype="dashed") +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	 geom_path(aes(x=Time,y=Scaled_Estimate,group=Model,color=Model)) +
	geom_point(data=nominal_dt[Category=="all"],aes(x=Time,y=Estimate,group=Model),size=2) +
	ylab("Estimate (relative)") +
	 theme_few(base_size = 20)
    ggsave(filename=paste0("p14.png"), plot = p14, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

p15 = index_2_dt %>%
	.[Model=="single"] %>%
	 .[,.(Model,Category,Time,Estimate,CV)] %>%
	 rbind(.,mv_dt,mvpl_dt) %>%
	 rbind(.,design_dt[Category!="all"]) %>%
	 	 .[,Scaled_Estimate:=(Estimate/mean(Estimate)),by=.(Model,Category)] %>%
	 .[,l95:=exp(log(Scaled_Estimate)-2*CV)] %>%
	 .[,u95:=exp(log(Scaled_Estimate)+2*CV)] %>%
	ggplot() +
	ylim(0,NA) +
	facet_wrap(~Category,scales="free_y") +
	geom_hline(yintercept=0) +
	 geom_hline(yintercept=1,linetype="dashed") +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	geom_path(aes(x=Time,y=Scaled_Estimate,group=Model,color=Model)) +
	geom_point(data=nominal_dt[Category!="all"],aes(x=Time,y=Estimate,group=Model),size=2) +
	ylab("Estimate (relative)") +
	theme_few(base_size = 20)
    ggsave(filename=paste0("p15.png"), plot = p15, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)