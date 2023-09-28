

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
library(scales)

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
    load(paste0(proj_dir,"VAST/model_runs/2023-06-14/2022_lgdl_mv_05_gearSP12_depth3.hardness3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/bfish_dt.RData"))
    load(paste0(proj_dir,"VAST/model_residual/2023-06-14/2022_lgdl_mv_05_gearSP12_depth3.hardness3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/resid_dt.RData"))
	load(paste0(proj_dir,"VAST/model_residual/2023-06-14/2022_lgdl_mv_05_gearSP12_depth3.hardness3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/residuals_dharma.RData"))

    # needed to define spatial domain and for predicting on to create index
	psu_table = fread(paste0(proj_dir,"Data/BFISH PSU lookup table.csv")) %>%
				.[,.(PSU,Island,lon_deg,lat_deg,STRATA,STRATA_2020,Depth_MEDIAN_m,med_slp,med_acr,BS_pct_over_136j,pctHB,pctHS)] %>%
				.[,substrate:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][1])] %>%
				.[,slope:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][2])]

  deep7_code_vec = tolower(c("ETCA","APRU","PRSI","HYQU","PRFI","PRZO","ETCO"))
	deep7_name_vec = c("Ehu (ETCA)", "Lehi (APRU)", "Kalekale (PRSI)", "Hapu'upu'u (HYQU)", "'Opakapaka (PRFI)", "Gindai (PRZO)", "Onaga (ETCO)")
    u_species = tolower(c("PRFI","ETCA","ETCO","PRSI","PRZO","HYQU","APRU"))
#_____________________________________________________________________________________________________________________________
# helper functions

  extract_discrete_abundance_effect = function(fit_internal,
                                                model_name,
                                                discrete_ab_variables = c("island","islandG"),
                                                target_species = c("prfi","etca","etco","prsi","przo","hyqu","apru"))
  {
    # checks
      discrete_ab_variables = intersect(discrete_ab_variables,colnames(fit_internal$covariate_data))
      keep_vec = rep(NA,length(discrete_ab_variables))
      for(i in 1:length(discrete_ab_variables))
      {
        keep_vec[i] =ifelse(sum(grep(paste0("\\b",discrete_ab_variables[i],"\\b"),fit_internal$X1_formula))>0,TRUE,FALSE)

      }
      discrete_ab_variables = discrete_ab_variables[keep_vec]

    if(fit_internal$ParHat[1] == "Model is not converged")
    {
      abundance_effect_dt = data.table(model_name=model_name,
                                       component=NA,
                                       species_cd=NA,
                                       category=NA,
                                       variable=NA,
                                       value=NA,
                                       value_link=NA,
                                       value_link_rescale=NA,
                                       value_rescale=NA)
    } else {
      ab_df_plot = fit_internal$covariate_data
   
      abundance_effect_dt.list = as.list(rep(NA,length(discrete_ab_variables)))
      for(i in 1:length(discrete_ab_variables))
      {
        abundance_effect_dt.list[[i]] = as.list(rep(NA,length(target_species)))
        tmp_cols = discrete_ab_variables[i]
        for(s in 1:length(target_species))
        {
          abundance_effect_dt.list[[i]][[s]] = as.data.table(ab_df_plot) %>%
            .[,..tmp_cols] %>%
            .[,value:=fit_internal$data_list$X1_gctp[,s,1,grep(discrete_ab_variables[i],dimnames(fit_internal$data_list$X1_gctp)[[4]])] %*% t(t(fit_internal$ParHat$gamma1_cp[s,grep(discrete_ab_variables[i],dimnames(fit_internal$data_list$X1_gctp)[[4]])]))] %>%
            .[,species_cd:=target_species[s]] %>%
            setnames(.,tmp_cols,"variable") %>%
            .[,value_link:=exp(value)] %>%
            .[,category:=tmp_cols] %>%
            .[,component:="1st"]
        }
        abundance_effect_dt.list[[i]] = rbindlist(abundance_effect_dt.list[[i]])
        rm(list="tmp_cols")
      }
      abundance_effect_dt = rbindlist(abundance_effect_dt.list) %>%
                            .[,model_name:=model_name] %>%
                            .[,value_link_rescale:=value_link/max(value_link),by=.(species_cd,category)] %>%
                            .[,value_rescale:=rescale(value),by=.(species_cd,category)] %>%
                            .[,.(model_name,component,species_cd,category,variable,value,value_link,value_link_rescale,value_rescale)] %>%
                            .[,category:=gsub("_sc","",category)]
     
    }
    return(unique(abundance_effect_dt))  
                    
  }


				# trim tails function comes from
				# https://stackoverflow.com/questions/44628130/ggplot2-dealing-with-extremes-values-by-setting-a-continuous-color-scale
				trim_tails <- function(range = c(-Inf, Inf)) scales::trans_new("trim_tails", 
	                transform = function(x) {
	                  force(range)
	                  desired_breaks <- scales::extended_breaks(n = 7)(x[x >= range[1] & x <= range[2]])
	                  break_increment <- diff(desired_breaks)[1]
	                  x[x < range[1]] <- range[1] - break_increment
	                  x[x > range[2]] <- range[2] + break_increment
	                  x
	                },
	                inverse = function(x) x,

	                breaks = function(x) {
	                  force(range)
	                  scales::extended_breaks(n = 7)(x)
	                },
	                format = function(x) {
	                  force(range)
	                  x[1] <- paste("<", range[1])
	                  x[length(x)] <- paste(">", range[2])
	                  x
	                })


#_____________________________________________________________________________________________________________________________
# plots
  # comparison
  nominal_agg_dt = bfish_dt %>%
               .[,.(year,species_cd,psu,weight_kg)] %>%
               .[,.(weight_kg=sum(weight_kg)),by=.(year,psu)] %>%
               .[,.(nominal=mean(weight_kg)),by=year] %>%
               .[,nominal_std:=nominal/mean(nominal)] %>%
               .[,year:=as.numeric(year)] %>%
               .[,category:="Total"] %>%
               .[,category:=factor(category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]


    nominal_sp_dt = bfish_dt %>%
               .[,.(year,species_cd,psu,weight_kg)] %>%
               .[,.(weight_kg=sum(weight_kg)),by=.(year,psu,species_cd)] %>%
               .[,.(nominal=mean(weight_kg)),by=.(year,species_cd)] %>%
               setnames(.,"species_cd","category") %>%
                .[,category:=deep7_name_vec[match(category,deep7_code_vec)]] %>%
                .[,category:=factor(category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))] %>%
               .[,nominal_std:=nominal/mean(nominal),by=category] %>%
               .[,year:=as.numeric(year)] 
    
    nominal_dt = rbind(nominal_agg_dt,nominal_sp_dt) %>% .[order(category)] %>%
               .[,panel:="Step 01"]


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
     geom_path(data=nominal_dt,aes(x=year,y=nominal_std),linetype="dashed",linewidth=1.5) +
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
      geom_path(data=nominal_dt[category=="Total"],aes(x=year,y=nominal_std),linetype="dashed",linewidth=1.5) +
      geom_path(data=p_panel2,aes(x=time,y=estimate_sc,group=model_name_plot),color="gray80",linewidth=1.5) +
      geom_path(data=p_panel,aes(x=time,y=estimate_sc,group=model_name_plot),color="gray40",linewidth=1.5) +
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
                .[,species_cd:=factor(species_cd,levels=deep7_code_vec,labels=deep7_name_vec)] %>%
                .[,species_cd:=factor(species_cd,levels=c("'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

 p = p %>%   
     ggplot() +
		 xlab("Year") +
		 facet_wrap(~species_cd,scales="free") +
     ylab("Relative influence") +
     geom_hline(yintercept=1,linetype="dotted") +
     geom_path(aes(x=as.numeric(year),y=influ,group=component,color=component),linewidth=1.5) +
     scale_colour_manual("Model\ncomponent",values=c(paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(uniqueN(p$component)-1))) +
		 theme_few(base_size=20)
    ggsave(filename=paste0("index_influ_a.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 12, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

t = bfish_combined_long_dt %>% .[,.(year,species_cd,gear_type)] %>% 
                               .[,.N,by=.(year,gear_type)] %>%
                               .[,N:=N/7] %>%
                               dcast(.,year~gear_type) %>%
                               .[,rf_prop:=research_fishing/(camera+research_fishing)] %>%
                               .[,std:=rf_prop/mean(rf_prop)]

# abundance covariates
 p = abundance_effect_dt %>%
              .[model_number==69] %>%
              setnames(.,c("variable_c","value"),c("variable","value_plot")) %>%
              .[,.(model_number,model_name_short,component,species_cd,category,variable,value_plot)] %>%
              na.omit(.) 
p = merge(p,index_summary_dt[,.(model_number,error_structure)],by="model_number") %>%
                .[,error_1st:=substr(error_structure,1,2)] %>%
                .[,value_link:=as.numeric(NA)] %>%
                .[error_1st=="pl",value_link:=1-exp(-exp(value_plot))] %>%
                .[error_1st=="lg",value_link:=boot::inv.logit(value_plot)] %>%
                .[,value_plot:=NULL] %>%
							  .[,value_link_rescale:=value_link/max(value_link),by=.(model_number,species_cd,category)] %>%
                setnames(.,c("value_link_rescale"),c("value_plot")) %>%
                .[,.(model_number,model_name_short,component,species_cd,category,variable,value_plot)] %>%
                na.omit(.)
 plot_empirical_dt_internal = copy(plot_empirical_dt) %>%
                          setnames(.,"empirical_sc_b","empirical_sc")


    p$model_name_short = factor(p$model_name_short,levels=unique(p$model_name_short),labels=c("09: Final model (lgdl)"))
    p = plot_empirical_dt_internal %>%
    .[gear_type=="camera"] %>%
    .[category %in% unique(p$category)] %>%
    .[,species_cd:=factor(species_cd,levels=c("'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))] %>%
    ggplot() +
    ylim(0,NA) +
    ylab("Relative encounter rate") +
	  xlab("Density covariate") +
	  facet_wrap(category~species_cd,scales="free",ncol=7) +
	  geom_rect(aes(xmin=breaks-x_width,xmax=breaks,ymin=0,ymax=empirical_sc,fill=gear_type),color="white") +
 	  geom_rect(data=plot_empirical_dt_internal[gear_type=="research_fishing" & category %in% unique(p$category)],aes(xmin=breaks,xmax=breaks+x_width,ymin=0,ymax=empirical_sc,fill=gear_type),color="white") +
    geom_line(data=p,aes(x=variable,y=value_plot,color=model_name_short,group=model_name_short)) +
    geom_hline(yintercept=0) +
    theme_few(base_size=20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.ticks.y=element_blank(),axis.text.y=element_blank(),strip.text.x = element_text(size = 10)) +
    scale_fill_manual("Gear type",values=c("gray60","gray80")) +
	  viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
    ggsave(filename=paste0("index_dens_covar_a.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 12, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

# discrete abundance effect

    # load(paste0(proj_dir,"VAST/model_runs/2023-06-14/2022_lgdl_mv_05_gearSP12_depth3.hardness3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/fit.RData"))
    discrete_ab_effect_dt = extract_discrete_abundance_effect(fit_internal = fit,model_name="09: Final model (lgdl)") %>%
                            unique(.) %>%
                            na.omit(.)

  discrete_ab_effect_dt$variable = factor(discrete_ab_effect_dt$variable,levels=c("Niihau","Kauai","Oahu","Maui Nui", "Big Island"))
  discrete_ab_effect_dt$species_cd=factor(discrete_ab_effect_dt$species_cd,levels=deep7_code_vec,labels=deep7_name_vec)
  

    p = plot_empirical_discrete_dt %>%
    .[category=="islandG"] %>%
    .[,species_cd:=factor(species_cd,levels=c("'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))] %>%
    ggplot() +
    ylab("Relative encounter rate") +
	  xlab("Island group") +
	  facet_wrap(category~species_cd,scales="free",ncol=7) +
	  geom_bar(aes(x=breaks,y=empirical_sc_b,fill=gear_type),stat="identity",color="white",position=position_dodge()) +
    geom_point(data=discrete_ab_effect_dt,aes(x=variable,y=value_link_rescale,color=model_name,group=model_name),position="jitter",size=2) +
    geom_hline(yintercept=0) +
    theme_few(base_size=20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.ticks.y=element_blank(),axis.text.y=element_blank(),strip.text.x = element_text(size = 10)) +
    scale_fill_manual("Gear type",values=c("gray60","gray80")) +
	  viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
    ggsave(filename=paste0("index_dens_covar_b.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 12, height = 6, units = c("in"),
            dpi = 300, limitsize = TRUE)
    
# catchability covariates

	    p = q_effect_dt %>% 
            .[model_number==69,.(component,species_cd,variable_d,value)] %>%
            .[component=="1st",value_link:=boot::inv.logit(value)] %>%
            .[component=="2nd",value_link:=exp(value)] %>%
            .[,species_cd:=sapply(species_cd,function(x)gsub(")","",gsub("(","",strsplit(as.character(x)," ")[[1]][2],fixed=TRUE),fixed=TRUE))] %>%
            .[,species_cd:=factor(species_cd,levels=c("PRFI","ETCA","ETCO","PRSI","PRZO","HYQU","APRU"))] %>%
            # change species cd to simple
			ggplot() +
			ylab("Relative Effect") +
			xlab("Gear type") +
			facet_grid(component~species_cd,scales="free_y") +
            geom_hline(yintercept=0) +
            geom_segment(aes(x=variable_d,xend=variable_d,y=0,yend=value_link)) +
			geom_point(aes(x=variable_d,y=value_link,fill=variable_d),shape=21,size=5) +
            theme_few(base_size=20) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            viridis::scale_fill_viridis("Gear",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
    ggsave(filename=paste0("index_q_covar_a.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 12, height = 6, units = c("in"),
            dpi = 300, limitsize = TRUE)

# spatial plots
# predict on annual psu table
			predict_knot = as.data.table(units::drop_units(fit$Report$D_gct)) %>%
				setnames(.,c("Site","Category","Time","value"),c("knot","species","year","density")) %>%
				.[,knot:=as.numeric(as.character(knot))] %>%
				.[,species:=factor(species,levels=c("prfi","etca","etco","prsi","przo","hyqu","apru"))] %>%
				.[,species_hw:=deep7_name_vec[match(species,deep7_code_vec)]] %>%
				.[,year:=as.numeric(as.character(year))] %>%
				.[,.(species,species_hw,year,knot,density)] %>%
				.[order(species,year,knot)]

                spatial_list = fit$spatial_list
                crs_ll = attr(spatial_list$loc_i,"origCRS")@projargs
                crs_eqd = attr(spatial_list$loc_i,"projCRS")@projargs
                	# needed for plotting and defining the barrier feature
                    hi_coast = rgdal::readOGR(dsn = paste0(proj_dir,"Data/GIS/Coastline"), layer = "Coastline")
                    hi_coast_sf = sf::st_as_sf(hi_coast)
                        hi_coast_eqd = sp::spTransform(hi_coast, crs_eqd)
                        hi_coast_eqd_sf = sf::st_as_sf(hi_coast_eqd) 

                    # convert to equi-distant coordinates for samples & psu
                        psu_coords = cbind(psu_table$lon_deg, psu_table$lat_deg)
                        psu_sp = sp::SpatialPoints(psu_coords)
                        sp::proj4string(psu_sp) = crs_ll
                        psu_sp_eqd = sp::spTransform(psu_sp, crs_eqd)
                        psu_table$lon_eqd = psu_sp_eqd@coords[,1]
                        psu_table$lat_eqd = psu_sp_eqd@coords[,2]
                        
				knot_coords = spatial_list$latlon_g[,2:1]
				knot_sp = sp::SpatialPoints(knot_coords)
				sp::proj4string(knot_sp) = crs_ll
				knot_sp_eqd = sp::spTransform(knot_sp, crs_eqd)
				knot_loc_dt = data.table(knot=1:nrow(knot_coords),lon=knot_sp@coords[,1],lat=knot_sp@coords[,2],lon_eqd=knot_sp_eqd@coords[,1],lat_eqd=knot_sp_eqd@coords[,2])
			predict_knot = merge(predict_knot,knot_loc_dt,by="knot")

			omega1_dt = as.data.table(fit$Report$Omega1_gc) %>%
				.[,knot:=1:nrow(fit$Report$Omega1_gc)] %>%
				melt(.,id.vars="knot") %>%
				setnames(.,c("variable","value"),c("species","omega1"))
			omega2_dt = as.data.table(fit$Report$Omega2_gc) %>%
				.[,knot:=1:nrow(fit$Report$Omega2_gc)] %>%
				melt(.,id.vars="knot") %>%
				setnames(.,c("variable","value"),c("species","omega2"))
			omega_dt = merge(omega1_dt,omega2_dt,by=c("knot","species"))

			epsilon1_dt = as.data.table(fit$Report$Epsilon1_gct) %>%
				.[,knot:=as.numeric(as.character(Site))] %>%
				.[,species:=factor(Category,levels=c("prfi","etca","etco","prsi","przo","hyqu","apru"))] %>%
				.[,year:=as.numeric(as.character(Time))] %>%
				setnames(.,"value","epsilon1") %>%
				.[,.(knot,species,year,epsilon1)]
			epsilon2_dt = as.data.table(fit$Report$Epsilon2_gct) %>%
				.[,knot:=as.numeric(as.character(Site))] %>%
				.[,species:=factor(Category,levels=c("prfi","etca","etco","prsi","przo","hyqu","apru"))] %>%
				.[,year:=as.numeric(as.character(Time))] %>%
				setnames(.,"value","epsilon2") %>%
				.[,.(knot,species,year,epsilon2)]
			epsilon_dt = merge(epsilon1_dt,epsilon2_dt,by=c("year","knot","species"))


			encounter_dt =	as.data.table(fit$Report$R1_gct) %>%
				.[,knot:=as.numeric(as.character(Site))] %>%
				.[,species:=factor(Category,levels=c("prfi","etca","etco","prsi","przo","hyqu","apru"))] %>%
				.[,year:=as.numeric(as.character(Time))] %>%
				setnames(.,"value","encounter_prob") %>%
				.[,.(knot,species,year,encounter_prob)]

			positive_dt =	as.data.table(fit$Report$R2_gct) %>%
				.[,knot:=as.numeric(as.character(Site))] %>%
				.[,species:=factor(Category,levels=c("prfi","etca","etco","prsi","przo","hyqu","apru"))] %>%
				.[,year:=as.numeric(as.character(Time))] %>%
				setnames(.,"value","positive_catch") %>%
				.[,.(knot,species,year,positive_catch)]
			component_dt = merge(encounter_dt,positive_dt,by=c("year","knot","species"))



			spatial_residual_dt = data.table(PSU=bfish_dt$psu,year=as.numeric(as.character(bfish_dt$year)),species=bfish_dt$species_cd,residual=residuals_dharma$scaledResiduals)

				psu_table$knot = 1:nrow(psu_table)
			predict_psu = merge(psu_table,predict_knot[,.(knot,species,species_hw,year,density)],by="knot",allow.cartesian=TRUE) %>%
						  .[,area_km2:=0.5^2] %>%
						  merge(.,omega_dt,by=c("knot","species")) %>%
						  merge(.,epsilon_dt,by=c("year","knot","species")) %>%
						  merge(.,component_dt,by=c("year","knot","species")) %>%
						 merge(.,spatial_residual_dt,by=c("PSU","year","species"),all.x=TRUE)

# filter for components that are effectively "turned-off"
			for(i in 1:length(u_species))
			{
				if(max(abs(predict_psu[species == u_species[i]]$omega1)) < 1e-8)
				{
					predict_psu = predict_psu %>%
								  .[species == u_species[i], omega1 := 0]
				}
				if(max(abs(predict_psu[species == u_species[i]]$epsilon1)) < 1e-8)
				{
					predict_psu = predict_psu %>%
								  .[species == u_species[i], epsilon1 := 0]
				}
				if(max(abs(predict_psu[species == u_species[i]]$omega2)) < 1e-8)
				{
					predict_psu = predict_psu %>%
								  .[species == u_species[i], omega2 := 0]
				}
				if(max(abs(predict_psu[species == u_species[i]]$epsilon2)) < 1e-8)
				{
					predict_psu = predict_psu %>%
								  .[species == u_species[i], epsilon2 := 0]
				}
			}

		p = copy(predict_psu) %>%
				.[,.(biomass=sum(density),lon_eqd=mean(lon_eqd),lat_eqd=mean(lat_eqd)),by=.(Island,PSU,knot,year)] %>%
				ggplot() + 
				facet_wrap(~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = biomass),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5,linewidth=0.5) +
				ggthemes::theme_few(base_size=20) + 
                theme(axis.text.x=element_blank(), #remove x axis labels
                    axis.ticks.x=element_blank(), #remove x axis ticks
                    axis.text.y=element_blank(),  #remove y axis labels
                    axis.ticks.y=element_blank()  #remove y axis ticks
                    ) +
				viridis::scale_color_viridis("Est. density",begin = 0.1,end = 0.8,direction = 1,option = "H",trans="log10")
	    ggsave(filename=paste0("sp_agg_dens.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 16, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)	

			# plot random effects distribution
				p = copy(predict_psu) %>%
					.[,.(knot,species,omega1,omega2,epsilon1,epsilon2)] %>%
					unique(.) %>%
					melt(.,id.vars=c("knot","species")) %>%
					.[,variable:=factor(variable,levels=c("omega1","omega2","epsilon1","epsilon2"))] %>%
                    .[,species:=factor(toupper(species),levels=c("PRFI","ETCA","ETCO","PRSI","PRZO","HYQU","APRU"))] %>%
                    ggplot() +
					facet_wrap(~species,scales="free_y",labeller=label_parsed) +
					xlab("Random effect") +
					ylab("Effect size") +
					geom_hline(yintercept=0) +
					geom_boxplot(aes(x=variable,y=value,fill=variable)) +
                    scale_x_discrete(labels = c(expression(omega[1]),expression(omega[2]),expression(epsilon[1]),expression(epsilon[2]))) +
					ggthemes::theme_few(base_size=20) + 
					viridis::scale_fill_viridis("Random\neffect",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE,labels = c(expression(omega[1]),expression(omega[2]),expression(epsilon[1]),expression(epsilon[2])))
			    ggsave(filename=paste0("sp_rf_dist.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 12, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

		p = copy(predict_psu) %>%
				.[,.(biomass=sum(density),lon_eqd=mean(lon_eqd),lat_eqd=mean(lat_eqd)),by=.(species,Island,PSU,knot,year)] %>%
				.[,species:=factor(toupper(species),levels=c("PRFI","ETCA","ETCO","PRSI","PRZO","HYQU","APRU"))] %>%
                ggplot() + 
				facet_grid(species~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = biomass),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5,linewidth=0.5) +
				ggthemes::theme_few(base_size=20) +
                theme(axis.text.x=element_blank(), #remove x axis labels
                    axis.ticks.x=element_blank(), #remove x axis ticks
                    axis.text.y=element_blank(),  #remove y axis labels
                    axis.ticks.y=element_blank()  #remove y axis ticks
                    ) + 
				viridis::scale_color_viridis("Est. density",begin = 0.1,end = 0.8,direction = 1,option = "H",trans="log10")
	    ggsave(filename=paste0("sp_sp_dens.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 16, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)


		p = copy(predict_psu) %>%
				.[,.(omega1=mean(omega1),omega2=mean(omega2),lon_eqd=mean(lon_eqd),lat_eqd=mean(lat_eqd)),by=.(species,Island,PSU,knot)] %>%
				.[,species:=factor(toupper(species),levels=c("PRFI","ETCA","ETCO","PRSI","PRZO","HYQU","APRU"))] %>%
                melt(.,id.vars=c("species","Island","PSU","knot","lon_eqd","lat_eqd")) %>%
                .[,variable:=factor(variable,levels=c("omega1","omega2"),labels=c(expression(omega[1]),expression(omega[2])))] %>%
                ggplot() + 
				facet_grid(species~variable,labeller=label_parsed) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = value),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5,linewidth=0.5) +
				ggthemes::theme_few(base_size=20) +
                theme(axis.text.x=element_blank(), #remove x axis labels
                    axis.ticks.x=element_blank(), #remove x axis ticks
                    axis.text.y=element_blank(),  #remove y axis labels
                    axis.ticks.y=element_blank()  #remove y axis ticks
                    ) + 
				scale_color_gradient2("Effect size",low = "blue",mid = "gray90",high ="red",trans = trim_tails(range = c(-1,1))) 
	    ggsave(filename=paste0("sp_sp_omega.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 9, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

		p = copy(predict_psu) %>%
				.[,.(eps=mean(epsilon1),lon_eqd=mean(lon_eqd),lat_eqd=mean(lat_eqd)),by=.(species,Island,PSU,knot,year)] %>%
				.[,species:=factor(toupper(species),levels=c("PRFI","ETCA","ETCO","PRSI","PRZO","HYQU","APRU"))] %>%
                ggplot() + 
				facet_grid(species~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = eps),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5,linewidth=0.5) +
				ggthemes::theme_few(base_size=20) +
                theme(axis.text.x=element_blank(), #remove x axis labels
                    axis.ticks.x=element_blank(), #remove x axis ticks
                    axis.text.y=element_blank(),  #remove y axis labels
                    axis.ticks.y=element_blank()  #remove y axis ticks
                    ) + 
				scale_color_gradient2(expression(epsilon[1]),low = "blue",mid = "gray90",high ="red",trans = trim_tails(range = c(-1,1))) 
	    ggsave(filename=paste0("sp_sp_eps1.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 16, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

		p = copy(predict_psu) %>%
				.[,.(eps=mean(epsilon2),lon_eqd=mean(lon_eqd),lat_eqd=mean(lat_eqd)),by=.(species,Island,PSU,knot,year)] %>%
				.[,species:=factor(toupper(species),levels=c("PRFI","ETCA","ETCO","PRSI","PRZO","HYQU","APRU"))] %>%
                ggplot() + 
				facet_grid(species~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = eps),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5,linewidth=0.5) +
				ggthemes::theme_few(base_size=20) +
                theme(axis.text.x=element_blank(), #remove x axis labels
                    axis.ticks.x=element_blank(), #remove x axis ticks
                    axis.text.y=element_blank(),  #remove y axis labels
                    axis.ticks.y=element_blank()  #remove y axis ticks
                    ) + 
				scale_color_gradient2(expression(epsilon[1]),low = "blue",mid = "gray90",high ="red",trans = trim_tails(range = c(-1,1))) 
	    ggsave(filename=paste0("sp_sp_eps2.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 16, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

		p = copy(predict_psu) %>%
				.[,.(val=mean(encounter_prob),lon_eqd=mean(lon_eqd),lat_eqd=mean(lat_eqd)),by=.(species,Island,PSU,knot,year)] %>%
				.[,species:=factor(toupper(species),levels=c("PRFI","ETCA","ETCO","PRSI","PRZO","HYQU","APRU"))] %>%
                ggplot() + 
				facet_grid(species~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = val),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5,linewidth=0.5) +
				ggthemes::theme_few(base_size=20) +
                theme(axis.text.x=element_blank(), #remove x axis labels
                    axis.ticks.x=element_blank(), #remove x axis ticks
                    axis.text.y=element_blank(),  #remove y axis labels
                    axis.ticks.y=element_blank()  #remove y axis ticks
                    ) + 
				viridis::scale_color_viridis("Encounter\nprobability",begin = 0.1,end = 0.8,direction = 1,option = "H")
	    ggsave(filename=paste0("sp_sp_encprob.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 16, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

		p = copy(predict_psu) %>%
				.[,.(val=mean(encounter_prob),lon_eqd=mean(lon_eqd),lat_eqd=mean(lat_eqd)),by=.(Island,PSU,knot,year)] %>%
                ggplot() + 
				facet_wrap(~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = val),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5,linewidth=0.5) +
				ggthemes::theme_few(base_size=20) +
                theme(axis.text.x=element_blank(), #remove x axis labels
                    axis.ticks.x=element_blank(), #remove x axis ticks
                    axis.text.y=element_blank(),  #remove y axis labels
                    axis.ticks.y=element_blank()  #remove y axis ticks
                    ) + 
				viridis::scale_color_viridis("Avg.\nEncounter\nprobability",begin = 0.1,end = 0.8,direction = 1,option = "H")
	    ggsave(filename=paste0("sp_agg_encprob.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 16, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

		p = copy(predict_psu) %>%
				.[,.(val=mean(positive_catch),lon_eqd=mean(lon_eqd),lat_eqd=mean(lat_eqd)),by=.(species,Island,PSU,knot,year)] %>%
				.[,species:=factor(toupper(species),levels=c("PRFI","ETCA","ETCO","PRSI","PRZO","HYQU","APRU"))] %>%
                ggplot() + 
				facet_grid(species~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = val),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5,linewidth=0.5) +
				ggthemes::theme_few(base_size=20) +
                theme(axis.text.x=element_blank(), #remove x axis labels
                    axis.ticks.x=element_blank(), #remove x axis ticks
                    axis.text.y=element_blank(),  #remove y axis labels
                    axis.ticks.y=element_blank()  #remove y axis ticks
                    ) + 
				viridis::scale_color_viridis("Positive\ncatch",begin = 0.1,end = 0.8,direction = 1,option = "H")
	    ggsave(filename=paste0("sp_sp_poscatch.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 16, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

		p = copy(predict_psu) %>%
				.[,.(val=mean(positive_catch),lon_eqd=mean(lon_eqd),lat_eqd=mean(lat_eqd)),by=.(Island,PSU,knot,year)] %>%
                ggplot() + 
				facet_wrap(~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = val),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5,linewidth=0.5) +
				ggthemes::theme_few(base_size=20) +
                theme(axis.text.x=element_blank(), #remove x axis labels
                    axis.ticks.x=element_blank(), #remove x axis ticks
                    axis.text.y=element_blank(),  #remove y axis labels
                    axis.ticks.y=element_blank()  #remove y axis ticks
                    ) + 
				viridis::scale_color_viridis("Avg.\npositive\ncatch",begin = 0.1,end = 0.8,direction = 1,option = "H")
	    ggsave(filename=paste0("sp_agg_poscatch.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 16, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

		p = copy(predict_psu) %>%
				.[,.(val=mean(residual,na.rm=TRUE),lon_eqd=mean(lon_eqd),lat_eqd=mean(lat_eqd)),by=.(species,Island,PSU,knot,year)] %>%
				.[,species:=factor(toupper(species),levels=c("PRFI","ETCA","ETCO","PRSI","PRZO","HYQU","APRU"))] %>%
                na.omit(.) %>%
                ggplot() + 
				facet_grid(species~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = val),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5,linewidth=0.5) +
				ggthemes::theme_few(base_size=20) +
                theme(axis.text.x=element_blank(), #remove x axis labels
                    axis.ticks.x=element_blank(), #remove x axis ticks
                    axis.text.y=element_blank(),  #remove y axis labels
                    axis.ticks.y=element_blank()  #remove y axis ticks
                    ) + 
				scale_color_gradient2("Avg.\nresidual",low = "blue",mid = "gray90",high ="red",midpoint = 0.5) 
	    ggsave(filename=paste0("sp_sp_resid.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 16, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)


		p = copy(predict_psu) %>%
				.[,.(val=mean(residual,na.rm=TRUE),lon_eqd=mean(lon_eqd),lat_eqd=mean(lat_eqd)),by=.(Island,PSU,knot,year)] %>%
                na.omit(.) %>%
                ggplot() + 
				facet_wrap(~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = val),size=1.5) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5,linewidth=0.5) +
				ggthemes::theme_few(base_size=20) +
                theme(axis.text.x=element_blank(), #remove x axis labels
                    axis.ticks.x=element_blank(), #remove x axis ticks
                    axis.text.y=element_blank(),  #remove y axis labels
                    axis.ticks.y=element_blank()  #remove y axis ticks
                    ) + 
				scale_color_gradient2("Avg.\nresidual",low = "blue",mid = "gray90",high ="red",midpoint = 0.5)
	    ggsave(filename=paste0("sp_agg_resid.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 16, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

 
  # stepwise for opakapaka
  p = index_dt %>% 
      .[model_number %in% c(69,72,73,74,79,77)] %>%
      .[category == "'Opakapaka (PRFI)" ] %>%
      .[,model_name_plot := factor(as.character(model_number),levels=as.character(c(69,72,73,74,79,77)),labels=c("01: Final MV","02: PRFI","03: PRFI (all lengths)","04a: Research fishing","04b: Data 01","05: Camera"))]
  
  p_list = as.list(rep(NA,3))
  p_list[[1]] = p[model_number %in% c(69,72,73,79)] %>% .[,panel:="Panel 01: Both gears"]
  p_list[[2]] = p[model_number %in% c(73,74)] %>% .[,panel:="Panel 02: Research fishing"]
  p_list[[3]] = p[model_number %in% c(73,79,77)] %>% .[,panel:="Panel 03: Camera"]
  p = rbindlist(p_list)
    p = p %>%   
      ggplot() +
      ylim(0,NA) +
      xlab("Year") +
      facet_wrap(~panel,scales="free") +
      ylab("Relative abundance") +
      geom_hline(yintercept=0) + 
      geom_hline(yintercept=1,linetype="dotted") +
      geom_ribbon(aes(x=time,ymin=l95_sc,ymax=u95_sc,group=model_name_plot,fill=model_name_plot),alpha=0.25) +
      geom_path(aes(x=time,y=estimate_sc,group=model_name_plot,color=model_name_plot),linewidth=1.5) +
      scale_colour_manual("Model",values=c(paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(uniqueN(p$model_name_plot)-1))) +
      scale_fill_manual("Model",values=c(paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(uniqueN(p$model_name_plot)-1))) +
      theme_few(base_size=20)
    ggsave(filename=paste0("paka_index_stepwise_a.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 12, height = 6, units = c("in"),
            dpi = 300, limitsize = TRUE)

  # sensitivity for opakapaka
  p = index_dt %>% 
      .[model_number %in% c(74,83,89,77,82,88)] %>%
      .[category == "'Opakapaka (PRFI)" ] %>%
      .[,model_name_plot := factor(as.character(model_number),levels=as.character(c(74,83,89,77,82,88)),labels=c("PRFI","+ moonphase","+ month","PRFI","+ moonphase","+ month"))]
  p_list = as.list(rep(NA,2))
  p_list[[1]] = p[model_number %in% c(74,83,89)] %>% .[,panel:="Panel 01: Research fishing"]
  p_list[[2]] = p[model_number %in% c(77,82,88)] %>% .[,panel:="Panel 02: Camera"]
   p = rbindlist(p_list)
    p = p %>%   
      ggplot() +
      ylim(0,NA) +
      xlab("Year") +
      facet_wrap(~panel,scales="free") +
      ylab("Relative abundance") +
      geom_hline(yintercept=0) + 
      geom_hline(yintercept=1,linetype="dotted") +
      geom_ribbon(aes(x=time,ymin=l95_sc,ymax=u95_sc,group=model_name_plot,fill=model_name_plot),alpha=0.25) +
      geom_path(aes(x=time,y=estimate_sc,group=model_name_plot,color=model_name_plot),linewidth=1.5) +
      scale_colour_manual("Model",values=c(paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(uniqueN(p$model_name_plot)-1))) +
      scale_fill_manual("Model",values=c(paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(uniqueN(p$model_name_plot)-1))) +
      theme_few(base_size=20)
    ggsave(filename=paste0("paka_index_sens_a.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 8, height = 6, units = c("in"),
            dpi = 300, limitsize = TRUE) 


  # sensitivity for full model
  p1= index_dt %>% 
      .[model_number %in% c(46,51,49)] %>%
      .[category == "Total" ] %>%
      .[,model_name_plot := factor(as.character(model_number),levels=as.character(c(46,51,49)),labels=c("01: Base (Step 07)","01: Comp. 1","01: Comp. 1&2"))] %>%
      .[,panel:="Panel 01: Vessel RE"]
  p2= index_dt %>% 
      .[model_number %in% c(57,69,62,66)] %>%
      .[category == "Total" ] %>%
      .[,model_name_plot := factor(as.character(model_number),levels=as.character(c(57,69,62,66)),labels=c("02: Base (Step 08)","02: lgdl","02: pldl","02: lgdg"))] %>%
      .[,panel:="Panel 02: Error structure"]
  p3= index_dt %>% 
      .[model_number %in% c(57,64,68,67,71)] %>%
      .[category == "Total" ] %>%
      .[,model_name_plot := factor(as.character(model_number),levels=as.character(c(57,64,68,67,71)),labels=c("03: Base (Step 08)","03: Treat. 01","03: Treat. 02","03: Treat. 03","03: Treat. 04"))] %>%
      .[,panel:="Panel 03: Data treatment"]
  p4= index_dt %>% 
      .[model_number %in% c(57,65)] %>%
      .[category == "Total" ] %>%
      .[,model_name_plot := factor(as.character(model_number),levels=as.character(c(57,65)),labels=c("04: Base (Step 08)","04: 60m radius"))] %>%
      .[,panel:="Panel 04: Area-swept"]
  p5= index_dt %>% 
      .[model_number %in% c(57,61)] %>%
      .[category == "Total" ] %>%
      .[,model_name_plot := factor(as.character(model_number),levels=as.character(c(57,61)),labels=c("05: Base (Step 08)","05: Outlier"))] %>%
      .[,panel:="Panel 04: Lehi outlier"]
 
  p_list = list(p1,p2,p3,p4,p5)
  p = rbindlist(p_list)
    p = p %>%   
      ggplot() +
      ylim(0,NA) +
      xlab("Year") +
      facet_wrap(~panel,scales="free") +
      ylab("Relative abundance") +
      geom_hline(yintercept=0) + 
      geom_hline(yintercept=1,linetype="dotted") +
      geom_ribbon(aes(x=time,ymin=l95_sc,ymax=u95_sc,group=model_name_plot,fill=model_name_plot),alpha=0.25) +
      geom_path(aes(x=time,y=estimate_sc,group=model_name_plot,color=model_name_plot),linewidth=1.5) +
      scale_colour_manual("Model",values=c(paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(2),paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(3),paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(4),paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(1),paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(1))) +
      scale_fill_manual("Model",values=c(paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(2),paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(3),paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(4),paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(1),paste0("gray",30),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(1))) +
      theme_few(base_size=20)
	ggsave(filename=paste0("index_sens_a.png"), plot = p, device = "png", path = plot_dir,
	  			scale = 1.25, width = 12, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

