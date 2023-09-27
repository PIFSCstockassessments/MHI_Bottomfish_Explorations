

# Nicholas Ducharme-Barth
# 2023/09/27
# plots for Tech Memo

# Copyright (c) 2023 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
library(data.table)
library(magrittr)
library(ggplot2)
library(viridis)
library(ggthemes)

#_____________________________________________________________________________________________________________________________
# define dirs
  proj_dir = "D:/HOME/SAP/2024_Deep7/"
  plot_dir = paste0(proj_dir,"Reports/tech_memo/")
  dir.create(plot_dir,recursive=TRUE)

#_____________________________________________________________________________________________________________________________
# load data
	index_dt = fread(file=paste0(proj_dir,"VAST/model_runs/comparison_plots/index_dt.csv")) %>%
            .[,category:=factor(category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
	index_summary_dt = fread(file=paste0(proj_dir,"VAST/model_runs/comparison_plots/index_summary_dt.csv")) 
	abundance_effect_dt = fread(file=paste0(proj_dir,"VAST/model_runs/comparison_plots/abundance_effect_dt.csv")) %>%
            .[,species_cd:=factor(species_cd,levels=c("'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
	q_effect_dt = fread(file=paste0(proj_dir,"VAST/model_runs/comparison_plots/q_effect_dt.csv")) %>%
            .[,species_cd:=factor(species_cd,levels=c("'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
		
  
  plot_empirical_dt = fread(file=paste0(proj_dir,"VAST/model_runs/comparison_plots/plot_empirical_ab_dt.csv"))
	plot_empirical_discrete_dt = fread(file=paste0(proj_dir,"VAST/model_runs/comparison_plots/plot_empirical_discrete_ab_dt.csv")) %>%
                                  .[,breaks:=factor(breaks,levels=c("Niihau","Kauai","Oahu","Maui Nui", "Big Island"))]
	load(file=paste0(proj_dir,"Data/",2022,"_05.bfish_combined_long_dt.RData"))
  load(paste0(proj_dir,"VAST/model_runs/2023-06-14/2022_lgdl_mv_05_gearSP12_depth3.hardness3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/fit.RData"))

  deep7_code_vec = tolower(c("ETCA","APRU","PRSI","HYQU","PRFI","PRZO","ETCO"))
	deep7_name_vec = c("Ehu (ETCA)", "Lehi (APRU)", "Kalekale (PRSI)", "Hapu'upu'u (HYQU)", "'Opakapaka (PRFI)", "Gindai (PRZO)", "Onaga (ETCO)")

#_____________________________________________________________________________________________________________________________
# plots
  # comparison
  p = index_dt %>% 
    #  .[category=="Total",category:="Deep7 (Total)"] %>%
     .[model_number %in% c(0,69)] %>%
     .[,model_name_plot := factor(as.character(model_number),levels=c("0","69"),labels=c("Design (Richards 2022)", "Model (this study)"))]
  p = p %>%   
     ggplot() +
		 ylim(0,NA) +
		 xlab("Year") +
		 facet_wrap(~category,scales="free") +
     ylab("Relative abundance") +
     geom_hline(yintercept=0) + 
     geom_hline(yintercept=1,linetype="dotted") +
     geom_ribbon(aes(x=time,ymin=l95_sc,ymax=u95_sc,group=model_name_plot,fill=model_name_plot),alpha=0.25) +
     geom_path(aes(x=time,y=estimate_sc,group=model_name_plot,color=model_name_plot),linewidth=1.5) +
     scale_colour_manual("Model",values=c(paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(uniqueN(p$model_name_plot)-1))) +
		 scale_fill_manual("Model",values=c(paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(uniqueN(p$model_name_plot)-1))) +
		 theme_few(base_size=20)
	ggsave(filename=paste0("index_comp_a.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1.25, width = 12, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

  # stepwise
  p = index_dt %>% 
      .[model_number %in% c(40,36:39,42,46,57,69)] %>%
      .[category == "Total" ] %>%
      .[,model_name_plot := factor(as.character(model_number),levels=as.character(c(40,36:39,42,46,57,69)),labels=c("01: Base (pldg)","02: + bs(depth,3)","03: + bs(hardness,3)","04: + bs(complexity,3)","05: + bs(slope,3)","06: + IslandG","07: + gearSP","08: - bs(complexity,3)","09: Final model (lgdl)"))] %>%
      .[,panel:=paste0("Step 0",as.numeric(model_name_plot))]
  p_lev = levels(p$model_name_plot)
  p_list2 = p_list = as.list(rep(NA,length(p_lev)))
  for(i in 1:length(p_lev))
  {
    if(i>1){
      p_list[[i]] = p[model_name_plot %in% p_lev[(i-1):i]] %>% .[,panel:=NULL] %>%.[,panel:=paste0("Step 0",i)]
    } else {
      p_list[[i]] = p[model_name_plot %in% p_lev[1:i]] %>% .[,panel:=NULL] %>%.[,panel:=paste0("Step 0",i)]
    }
    p_list2[[i]] = p[model_name_plot %in% p_lev[1:i]] %>% .[,panel:=NULL] %>%.[,panel:=paste0("Step 0",i)]
  }
  p_panel = rbindlist(p_list)
  p_panel2 = rbindlist(p_list2)
    p = p %>%   
      ggplot() +
      ylim(0,NA) +
      xlab("Year") +
      facet_wrap(~panel,scales="free") +
      ylab("Relative abundance") +
      geom_hline(yintercept=0) + 
      geom_hline(yintercept=1,linetype="dotted") +
      geom_path(data=p_panel2,aes(x=time,y=estimate_sc,group=model_name_plot),color="gray80",linewidth=1.5) +
      geom_path(data=p_panel,aes(x=time,y=estimate_sc,group=model_name_plot),color="gray60",linewidth=1.5) +
      geom_ribbon(aes(x=time,ymin=l95_sc,ymax=u95_sc,group=model_name_plot,fill=model_name_plot),alpha=0.25) +
      geom_path(aes(x=time,y=estimate_sc,group=model_name_plot,color=model_name_plot),linewidth=1.5) +
      scale_colour_manual("Model",values=c(paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(uniqueN(p$model_name_plot)-1))) +
      scale_fill_manual("Model",values=c(paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(uniqueN(p$model_name_plot)-1))) +
      theme_few(base_size=20)
    ggsave(filename=paste0("index_stepwise_a.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 12, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

# influ
p = q_effect_dt %>% .[model_number==69,.(component,species_cd,variable_d,value)] %>%
                setnames(.,"variable_d","gear_type") %>%
                .[,species_cd:=sapply(species_cd,function(x)tolower(gsub(")","",gsub("(","",strsplit(as.character(x)," ")[[1]][2],fixed=TRUE),fixed=TRUE)))] %>%
                merge(bfish_combined_long_dt[,.(year,species_cd,gear_type)],.,by=c("species_cd","gear_type"),allow.cartesian=TRUE) %>%
                .[,rho:=mean(value),by=.(component,species_cd)] %>%
                .[,diff:=value-rho] %>%
                .[,.(raw_influ=mean(diff)),by=.(year,species_cd,component)] %>%
                .[order(component,species_cd,year)] %>%
                .[,influ:=exp(raw_influ)] %>%
                .[,species_cd:=factor(species_cd,levels=deep7_code_vec,labels=deep7_name_vec)]
 p = p %>%   
     ggplot() +
		 xlab("Year") +
		 facet_wrap(~species_cd,scales="free") +
     ylab("Relative influence") +
     geom_hline(yintercept=1,linetype="dotted") +
     geom_path(aes(x=year,y=influ,group=component,color=component),linewidth=1.5) +
     scale_colour_manual("Model\ncomponent",values=c(paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(uniqueN(p$component)-1))) +
		 theme_few(base_size=20)
    ggsave(filename=paste0("index_influ_a.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 12, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)
