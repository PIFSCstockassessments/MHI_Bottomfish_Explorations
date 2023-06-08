

# Nicholas Ducharme-Barth
# 2023/06/04
# Comparison plots

# Copyright (c) 2023 Nicholas Ducharme-Barth
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
	design_dt = fread(paste0(proj.dir,"Data/2022_design_based_estimates.csv")) %>%
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
  

#_____________________________________________________________________________________________________________________________
# plot 5) compare multi-species with catchability covariates
# data type: 05
    dt_01 = fread(paste0(proj.dir,"VAST/model_runs/2022-12-29/2021_pldg_mv_05_v_v_TRUE_7.5_FALSE_TRUE_pit_xval/index_dt.csv")) %>%
                     .[,Model:="01 base"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    dt_02 = fread(paste0(proj.dir,"VAST/model_runs/2023-01-03/2021_pldg_mv_05_gear_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="02 gear"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    dt_03 = fread(paste0(proj.dir,"VAST/model_runs/2023-01-03/2021_pldg_mv_05_gearSP_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="03 gearSP"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    dt_04 = fread(paste0(proj.dir,"VAST/model_runs/2023-01-05/2021_pldg_mv_05_gearSP_v_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="04 gearSP (fine_scale)"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

    dt_05 = fread(paste0(proj.dir,"VAST/model_runs/2023-01-05/2021_pldg_mv_05_gearSP_depth_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="05 gearSP & depth3 (fine_scale)"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

    dt_06 = fread(paste0(proj.dir,"VAST/model_runs/2023-01-06/2021_pldg_mv_05_gearSP_depth.hardness_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="06 gearSP & depth3.hardness (fine_scale)"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

    dt_07 = fread(paste0(proj.dir,"VAST/model_runs/2023-01-05/2021_pldg_mv_05_gearSP_depth.complexity.hardness_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="07 gearSP & depth3.hardness.complexity (fine_scale)"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

    dt_08a = fread(paste0(proj.dir,"VAST/model_runs/2023-01-06/2021_pldg_mv_05_gearSP_depth3.complexity3.hardness3.SCALED_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="08a gearSP & depth3.hardness3.complexity3 (fine_scale)"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

    dt_08b = fread(paste0(proj.dir,"VAST/model_runs/2023-01-06/2021_pldg_mv_05_gear.vesselE_depth3.complexity3.hardness3.SCALED_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="08b gear.vesselE & depth3.hardness3.complexity3 (fine_scale)"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

    dt_08c = fread(paste0(proj.dir,"VAST/model_runs/2023-01-06/2021_pldg_mv_05_gearSP.vesselEP_depth.complexity.hardness.SCALED_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="08c gearSP.vesselEP & depth.hardness.complexity (fine_scale)"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]


    p5 = rbind(design_dt,dt_01,dt_02,dt_03,dt_04,dt_05,dt_06,dt_07,dt_08a,dt_08b,dt_08c) %>%
            .[Model=="design",Model:="00 design"] %>%
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
        ggsave(filename=paste0("p5_catchability_covariates.png"), plot = p5, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)                        

    p5 = rbind(design_dt,dt_01,dt_02,dt_03,dt_04,dt_05,dt_06,dt_07,dt_08a,dt_08b,dt_08c) %>%
            .[Model=="design",Model:="00 design"] %>%
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
        ggsave(filename=paste0("p5_catchability_covariates_relative.png"), plot = p5, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)   
   
   
    p5 = rbind(design_dt,dt_01,dt_02,dt_03,dt_04,dt_05,dt_06,dt_07,dt_08a,dt_08b,dt_08c) %>%
            .[Model=="design",Model:="00 design"] %>%
            ggplot() +
			ylim(0,NA) +
			ylab("Predicted biomass (millions lbs)") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("p5_catchability_covariates_noRibbon.png"), plot = p5, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)                        

    p5 = rbind(design_dt,dt_01,dt_02,dt_03,dt_04,dt_05,dt_06,dt_07,dt_08a,dt_08b,dt_08c) %>%
            .[Model=="design",Model:="00 design"] %>%
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
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("p5_catchability_covariates_relative_noRibbon.png"), plot = p5, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)   

        p5.list = list(dt_01,dt_02,dt_03,dt_04,dt_05,dt_06,dt_07,dt_08a,dt_08b,dt_08c)

        for(i in seq_along(p5.list))
        {
                    p5 = rbind(design_dt,rbindlist(p5.list[1:i])) %>%
                .[Model=="design",Model:="00 design"] %>%
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
            ggsave(filename=paste0("p5_catchability_covariates.",i,".png"), plot = p5, device = "png", path = dir_plot,
                    scale = 1.25, width = 16, height = 9, units = c("in"),
                    dpi = 300, limitsize = TRUE)                        

                    p5 = rbind(design_dt,rbindlist(p5.list[1:i])) %>%
                .[Model=="design",Model:="00 design"] %>%
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
            ggsave(filename=paste0("p5_catchability_covariates_relative.",i,".png"), plot = p5, device = "png", path = dir_plot,
                    scale = 1.25, width = 16, height = 9, units = c("in"),
                    dpi = 300, limitsize = TRUE)   
    
    
                    p5 = rbind(design_dt,rbindlist(p5.list[1:i])) %>%
                .[Model=="design",Model:="00 design"] %>%
                ggplot() +
                ylim(0,NA) +
                ylab("Predicted biomass (millions lbs)") +
                xlab("Year") +
                facet_wrap(~Category,scales="free") +
                geom_hline(yintercept=0) +
                geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
                viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
                viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
                theme_few(base_size=20)
            ggsave(filename=paste0("p5_catchability_covariates_noRibbon.",i,".png"), plot = p5, device = "png", path = dir_plot,
                    scale = 1.25, width = 16, height = 9, units = c("in"),
                    dpi = 300, limitsize = TRUE)                        

                    p5 = rbind(design_dt,rbindlist(p5.list[1:i])) %>%
                .[Model=="design",Model:="00 design"] %>%
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
                geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
                viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
                viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
                theme_few(base_size=20)
            ggsave(filename=paste0("p5_catchability_covariates_relative_noRibbon.",i,".png"), plot = p5, device = "png", path = dir_plot,
                    scale = 1.25, width = 16, height = 9, units = c("in"),
                    dpi = 300, limitsize = TRUE) 
        }


#_____________________________________________________________________________________________________________________________
# plot 6) compare 2022 data models
# data type: 05
    dt_01 = fread(paste0(proj.dir,"VAST/model_runs/2023-06-02/2022_pldg_mv_05_v_v_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="01 base"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    dt_02 = fread(paste0(proj.dir,"VAST/model_runs/2023-06-03/2022_pldg_mv_05_v_depth3_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="02 d3"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    dt_03 = fread(paste0(proj.dir,"VAST/model_runs/2023-06-03/2022_pldg_mv_05_v_depth3.hardness3_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="03 d3.h3"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    dt_04 = fread(paste0(proj.dir,"VAST/model_runs/2023-06-03/2022_pldg_mv_05_v_depth3.hardness3.complexity3_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="04 d3.h3.c3"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    dt_05 = fread(paste0(proj.dir,"VAST/model_runs/2023-06-03/2022_pldg_mv_05_v_depth3.hardness3.complexity3.slope3_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="05 d3.h3.c3.s3"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    dt_06a = fread(paste0(proj.dir,"VAST/model_runs/2023-06-04/2022_pldg_mv_05_v_depth3.hardness3.complexity3.slope3.island_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="06a d3.h3.c3.s3.island"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    dt_06b = fread(paste0(proj.dir,"VAST/model_runs/2023-06-04/2022_pldg_mv_05_v_depth3.hardness3.complexity3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="06b d3.h3.c3.s3.islandG"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
                                                                                            
    p6 = rbind(design_dt,dt_01,dt_02,dt_03,dt_04,dt_05,dt_06a,dt_06b) %>%
            .[Model=="design",Model:="00 design"] %>%
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
        ggsave(filename=paste0("p6_2022_index.png"), plot = p6, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)                        

    p6 = rbind(design_dt,dt_01,dt_02,dt_03,dt_04,dt_05,dt_06a,dt_06b) %>%
            .[Model=="design",Model:="00 design"] %>%
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
        ggsave(filename=paste0("p6_2022_index_relative.png"), plot = p6, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)   
   
   
    p6 = rbind(design_dt,dt_01,dt_02,dt_03,dt_04,dt_05,dt_06a,dt_06b) %>%
            .[Model=="design",Model:="00 design"] %>%
            ggplot() +
			ylim(0,NA) +
			ylab("Predicted biomass (millions lbs)") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("p6_2022_index_noRibbon.png"), plot = p6, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)                        

    p6 = rbind(design_dt,dt_01,dt_02,dt_03,dt_04,dt_05,dt_06a,dt_06b) %>%
            .[Model=="design",Model:="00 design"] %>%
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
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("p6_2022_index_relative_noRibbon.png"), plot = p6, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)   


#_____________________________________________________________________________________________________________________________
# plot 7) compare 2022 data models
# data type: 05
    dt_01 = fread(paste0(proj.dir,"VAST/model_runs/2023-06-04/2022_pldg_mv_05_v_depth3.hardness3.complexity3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="01 base"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    dt_02 = fread(paste0(proj.dir,"VAST/model_runs/2023-06-04/2022_pldg_mv_05_gear1_depth3.hardness3.complexity3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="02 gear1"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    dt_03 = fread(paste0(proj.dir,"VAST/model_runs/2023-06-05/2022_pldg_mv_05_gearSP1.5_depth3.hardness3.complexity3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="03 gearSP1.5"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
     dt_04 = fread(paste0(proj.dir,"VAST/model_runs/2023-06-04/2022_pldg_mv_05_gearSP1_depth3.hardness3.complexity3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="04 gearSP1"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    dt_05 = fread(paste0(proj.dir,"VAST/model_runs/2023-06-04/2022_pldg_mv_05_gear12_depth3.hardness3.complexity3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="05 gear12"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
    dt_06 = fread(paste0(proj.dir,"VAST/model_runs/2023-06-05/2022_pldg_mv_05_gearSP12.5_depth3.hardness3.complexity3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="06 gearSP12.5"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
     dt_07 = fread(paste0(proj.dir,"VAST/model_runs/2023-06-04/2022_pldg_mv_05_gearSP12_depth3.hardness3.complexity3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="07 gearSP12"] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
                                                                                                                                                                                           
    p7 = rbind(design_dt,dt_01,dt_02,dt_03,dt_04,dt_05,dt_06,dt_07) %>%
            .[Model=="design",Model:="00 design"] %>%
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
        ggsave(filename=paste0("p7_2022_index.png"), plot = p7, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)                        

    p7 = rbind(design_dt,dt_01,dt_02,dt_03,dt_04,dt_05,dt_06,dt_07) %>%
            .[Model=="design",Model:="00 design"] %>%
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
        ggsave(filename=paste0("p7_2022_index_relative.png"), plot = p7, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)   
   
   
    p7 = rbind(design_dt,dt_01,dt_02,dt_03,dt_04,dt_05,dt_06,dt_07) %>%
            .[Model=="design",Model:="00 design"] %>%
            ggplot() +
			ylim(0,NA) +
			ylab("Predicted biomass (millions lbs)") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("p7_2022_index_noRibbon.png"), plot = p7, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)                        

    p7 = rbind(design_dt,dt_01,dt_02,dt_03,dt_04,dt_05,dt_06,dt_07) %>%
            .[Model=="design",Model:="00 design"] %>%
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
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("p7_2022_index_relative_noRibbon.png"), plot = p7, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)   

        