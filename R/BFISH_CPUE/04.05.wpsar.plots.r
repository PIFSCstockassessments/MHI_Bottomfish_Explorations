

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
  plot_dir = paste0(proj_dir,"Reports/wpsar/")
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
 
# abundance covariates
 p = abundance_effect_dt %>%
              .[model_number %in% c(74,77)] %>%
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
                na.omit(.) %>%
                .[species_cd == "'Opakapaka (PRFI)"] %>%
                .[,species_cd:=factor(species_cd,levels=c("'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

plot_empirical_dt_internal = copy(plot_empirical_dt) %>%
                          setnames(.,"empirical_sc_b","empirical_sc") %>%
                          .[,species_cd:=factor(species_cd,levels=c("'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]


    p$model_name_short = factor(p$model_name_short,levels=unique(p$model_name_short),labels=c("Research fishing","Camera"))
    p = plot_empirical_dt_internal %>%
    .[species_cd == "'Opakapaka (PRFI)"] %>%
    .[gear_type=="camera"] %>%
    .[category %in% unique(p$category)]  %>%
    ggplot() +
    ylim(0,NA) +
    ylab("Relative encounter rate") +
	  xlab("Density covariate") +
	  facet_wrap(category~species_cd,scales="free",ncol=7) +
	  geom_rect(aes(xmin=breaks-x_width,xmax=breaks,ymin=0,ymax=empirical_sc,fill=gear_type),color="white") +
 	  geom_rect(data=plot_empirical_dt_internal[species_cd == "'Opakapaka (PRFI)"&gear_type=="research_fishing" & category %in% unique(p$category)],aes(xmin=breaks,xmax=breaks+x_width,ymin=0,ymax=empirical_sc,fill=gear_type),color="white") +
    geom_line(data=p,aes(x=variable,y=value_plot,color=model_name_short,group=model_name_short)) +
    geom_hline(yintercept=0) +
    theme_few(base_size=20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.ticks.y=element_blank(),axis.text.y=element_blank(),strip.text.x = element_text(size = 10)) +
    scale_fill_manual("Gear type",values=c("gray60","gray80")) +
	  viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
    ggsave(filename=paste0("index_dens_covar_a_opakapaka.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 12, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

 # stepwise
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
      .[model_number %in% c(92,93,91,94,95,96,97)] %>%
      .[category == "Total" ] %>%
      .[,model_name_plot := factor(as.character(model_number),levels=as.character(c(92,93,91,94,95,96,97)),labels=c("01: Base (lgdl)","02: + bs(depth,3)","03: + bs(hardness,3)","04: + bs(slope,3)","05: + IslandG","06: + gearSP","07: -spatial RE"))] %>%
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
    ggsave(filename=paste0("index_stepwise_lgdl.png"), plot = p, device = "png", path = plot_dir,
            scale = 1.25, width = 12, height = 9, units = c("in"),
            dpi = 300, limitsize = TRUE)

            