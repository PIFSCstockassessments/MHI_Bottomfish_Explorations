

# Nicholas Ducharme-Barth
# 12/27/2022
# Comparison plots

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
	dir_plot = paste0(proj.dir,"VAST/model_runs/comparison_plots/")
	dir.create(dir_plot,recursive = TRUE)

    deep7_code_vec = tolower(c("ETCA","APRU","PRSI","HYQU","PRFI","PRZO","ETCO"))
	deep7_name_vec = c("Ehu (ETCA)", "Lehi (APRU)", "Kalekale (PRSI)", "Hapu'upu'u (HYQU)", "'Opakapaka (PRFI)", "Gindai (PRZO)", "Onaga (ETCO)")
		

#_____________________________________________________________________________________________________________________________
# bring in indices
    # design based bfish index data is in kg so convert to lbs (millions)
	design_dt = fread(paste0(proj.dir,"Data/2021_design_based_estimates.csv")) %>%
				.[,.(Model,Category,Time,Estimate,CV)] %>%
                .[,Estimate:=Estimate*2.20462262185] %>%
				.[,Estimate:=Estimate/1000000] %>%
				.[,l95:=exp(log(Estimate)-2*sqrt(log(CV^2+1)))] %>%
				.[,u95:=exp(log(Estimate)+2*sqrt(log(CV^2+1)))] %>%
                .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                .[is.na(Category),Category:="Total"] %>%
                .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

#_____________________________________________________________________________________________________________________________
# plot 1) apru single species, poisson-link delta-gamma model
# data type: 05
# with and without exclusion of the large lehi sample

    # bring in the index_dt from each model fit
    # Note: index_dt is in lbs (millions)

    apru_filter_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-27/2021_pldg_apru_05_v_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="model - filter"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

    apru_nofilter_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-27/2021_pldg_apru_05_v_v_FALSE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="model - no filter"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

    p1 = rbind(design_dt,apru_filter_dt,apru_nofilter_dt) %>%
            .[Category == "Lehi (APRU)"] %>%
            ggplot() +
			ylim(0,NA) +
			ylab("Predicted biomass (millions lbs)") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("p1_apru_filter_comp.png"), plot = p1, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)


#_____________________________________________________________________________________________________________________________
# plot 2) compare single species to multispecies, poisson-link delta-gamma model
# data type: 05

    apru_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-27/2021_pldg_apru_05_v_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="single"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

    hyqu_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-27/2021_pldg_hyqu_05_v_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="single"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    przo_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-27/2021_pldg_przo_05_v_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="single"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    prsi_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-27/2021_pldg_prsi_05_v_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="single"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    etco_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-27/2021_pldg_etco_05_v_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="single"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    etca_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-27/2021_pldg_etca_05_v_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="single"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    prfi_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-27/2021_pldg_prfi_05_v_v_TRUE_7.5_FALSE_TRUE_pit_xval/index_dt.csv")) %>%
                     .[,Model:="single"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

    mv_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-29/2021_pldg_mv_05_v_v_TRUE_7.5_FALSE_TRUE_pit_xval/index_dt.csv")) %>%
                     .[,Model:="mv"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
	
	total_dt = rbind(apru_dt,hyqu_dt,prfi_dt,prsi_dt,etco_dt,etca_dt,przo_dt) %>%
			.[,.(Estimate=sum(Estimate)),by=.(Model,Time)] %>%
			.[,Category:="Total"] %>%
			.[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))] %>%
			.[,CV := NA] %>%
			.[,l95:=NA] %>%
			.[,u95:=NA] %>%
			.[,.(Model,Category,Time,Estimate,CV,l95,u95)]

	p2 = rbind(mv_dt,total_dt,apru_dt,hyqu_dt,prfi_dt,prsi_dt,etco_dt,etca_dt,przo_dt,design_dt) %>%
            ggplot() +
			ylim(0,NA) +
			ylab("Predicted biomass (millions lbs)") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("p2_single_vs_mv.png"), plot = p2, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)



#_____________________________________________________________________________________________________________________________
# plot 3) compare multivariate across data types 01 - 05, poisson-link delta-gamma model

    mv_05_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-29/2021_pldg_mv_05_v_v_TRUE_7.5_FALSE_TRUE_pit_xval/index_dt.csv")) %>%
                     .[,Model:="mv - 05"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    mv_04_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-28/2021_pldg_mv_04_v_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="mv - 04 "] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    mv_03_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-28/2021_pldg_mv_03_v_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="mv - 03"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
		
    mv_02_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-28/2021_pldg_mv_02_v_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="mv - 02"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
		
    mv_01_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-28/2021_pldg_mv_01_v_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="mv - 01"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
		
	p3 = rbind(mv_01_dt,mv_02_dt,mv_03_dt,mv_04_dt,mv_05_dt,design_dt) %>%
            ggplot() +
			ylim(0,NA) +
			ylab("Predicted biomass (millions lbs)") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("p3_data_type.png"), plot = p3, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)	

	p3 = rbind(mv_01_dt,mv_02_dt,mv_03_dt,mv_04_dt,mv_05_dt,design_dt) %>%
			.[,Estimate:=Estimate/mean(Estimate),by=.(Model,Category)] %>%
			.[,l95:=exp(log(Estimate)-2*sqrt(log(CV^2+1)))] %>%
			.[,u95:=exp(log(Estimate)+2*sqrt(log(CV^2+1)))] %>%
            ggplot() +
			ylim(0,NA) +
			ylab("Relative biomass") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_hline(yintercept=1,linetype="dashed") +
			geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("p3_data_type_relative.png"), plot = p3, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)	

#_____________________________________________________________________________________________________________________________
# plot 4) compare multi-species across error structures
# data type: 05
 	pldg_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-29/2021_pldg_mv_05_v_v_TRUE_7.5_FALSE_TRUE_pit_xval/index_dt.csv")) %>%
                     .[,Model:="pldg"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
 	lgdg_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-29/2021_lgdg_mv_05_v_v_TRUE_7.5_FALSE_TRUE_pit_xval/index_dt.csv")) %>%
                     .[,Model:="lgdg"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
 	pldl_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-29/2021_pldl_mv_05_v_v_TRUE_7.5_FALSE_TRUE_pit_xval/index_dt.csv")) %>%
                     .[,Model:="pldl"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
  	plgg_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-30/2021_plgg_mv_05_v_v_TRUE_7.5_FALSE_TRUE_pit_xval/index_dt.csv")) %>%
                     .[,Model:="plgg"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

	p4 = rbind(pldg_dt,lgdg_dt,pldl_dt,plgg_dt,design_dt) %>%
            ggplot() +
			ylim(0,NA) +
			ylab("Predicted biomass (millions lbs)") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_hline(yintercept=1,linetype="dashed") +
			geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("p4_error_structure.png"), plot = p4, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)   
  
  	p4 = rbind(pldg_dt,lgdg_dt,pldl_dt,plgg_dt,design_dt) %>%
			.[,Estimate:=Estimate/mean(Estimate),by=.(Model,Category)] %>%
			.[,l95:=exp(log(Estimate)-2*sqrt(log(CV^2+1)))] %>%
			.[,u95:=exp(log(Estimate)+2*sqrt(log(CV^2+1)))] %>%
            ggplot() +
			ylim(0,NA) +
			ylab("Relative biomass") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_hline(yintercept=1,linetype="dashed") +
			geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("p4_error_structure_relative.png"), plot = p4, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)   
  

