

# Nicholas Ducharme-Barth
# 2023/06/04
# effects plots

# Copyright (c) 2023 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(ggplot2)
	library(ggthemes)
  library(scales)

#_____________________________________________________________________________________________________________________________
# set working directory
	proj_dir = "D:/HOME/SAP/2024_Deep7/"
	dir_vec = dir(paste0(proj_dir,"VAST/model_runs/"),recursive=TRUE)
	ests_vec = dir_vec[grep("/fit.RData",dir_vec,fixed=TRUE)]
  ests_vec = ests_vec[grep("/2022_",ests_vec,fixed=TRUE)]

  deep7_code_vec = tolower(c("ETCA","APRU","PRSI","HYQU","PRFI","PRZO","ETCO"))
	deep7_name_vec = c("Ehu (ETCA)", "Lehi (APRU)", "Kalekale (PRSI)", "Hapu'upu'u (HYQU)", "'Opakapaka (PRFI)", "Gindai (PRZO)", "Onaga (ETCO)")


#_____________________________________________________________________________________________________________________________
# load data
    data_flag = 2022
    data_treatment = "05"
	  load(file=paste0(proj_dir,"Data/",data_flag,"_",data_treatment,".bfish_combined_long_dt.RData"))

    psu_table = fread(paste0(proj_dir,"Data/BFISH PSU lookup table.csv")) %>%
        .[STRATA=="SB_H_S",STRATA:="SB_A_S"] %>%
				.[,.(PSU,STRATA,Island,Depth_MEDIAN_m,med_acr,BS_pct_over_136j,pctHS)] %>%
				setnames(.,c("PSU","STRATA","Island","Depth_MEDIAN_m","med_acr","BS_pct_over_136j","pctHS"),c("psu","strata","island","depth","complexity","hardness","slope")) %>%
				.[complexity>300,complexity:=300] %>%
        .[is.na(complexity),complexity:=12.20] %>%
        .[is.na(hardness),hardness:=0.094888] %>%
				.[,islandG:=ifelse(island=="Niihau","Kauai",island)] %>%
				.[,islandG:=factor(islandG,levels=c("Big Island","Maui Nui","Oahu","Kauai"))]
    
    plot_empirical_dt = bfish_combined_long_dt %>%
                        .[,.(model_sampling_unit,psu,gear_type,species_cd,weight_kg)] %>%
                        merge(.,psu_table,by="psu") %>%
                        setnames(.,"weight_kg","binary") %>%
                        .[,binary:=ifelse(binary>0,1,0)] %>%
                        melt(.,id.vars=c("psu","model_sampling_unit","gear_type","strata","island","islandG","species_cd","binary")) %>%
                        setnames(.,"variable","category") %>%
                        .[,breaks:=cut(value,breaks=10,right=FALSE,include.lowest=TRUE),by=category] %>%
                        .[,breaks:=sapply(breaks,function(x)as.numeric(gsub("[","",strsplit(as.character(x),",")[[1]][1],fixed=TRUE)))] %>%
                        .[breaks<0,breaks:=0] %>%
                        .[,empirical:=mean(binary),by=.(gear_type,species_cd,category,breaks)] %>%
                        .[,.(gear_type,species_cd,category,breaks,empirical)] %>%
                        unique(.) 
  
    plot_empirical_dt2 = plot_empirical_dt %>%                    
                        .[,.(category,breaks)] %>%
                        unique(.) %>%
                        .[order(category,breaks)] %>%
                        .[,x_width:=0.5*abs(mean(diff(breaks))),by=category]
    plot_empirical_dt = merge(plot_empirical_dt,plot_empirical_dt2) %>%
                        .[,y_max_a:=max(empirical),by=.(category,species_cd,gear_type)] %>%
                        .[,empirical_sc_a:=empirical/y_max_a] %>%
                        .[,y_max_b:=max(empirical),by=.(category,species_cd)] %>%
                        .[,empirical_sc_b:=empirical/y_max_b] %>%
                        .[,species_cd:=factor(species_cd,levels=deep7_code_vec,labels=deep7_name_vec)]

    plot_empirical_discrete_dt = bfish_combined_long_dt %>%
                        .[,.(model_sampling_unit,psu,gear_type,species_cd,weight_kg)] %>%
                        merge(.,psu_table[,.(psu,island,islandG)],by="psu") %>%
                        setnames(.,"weight_kg","binary") %>%
                        .[,binary:=ifelse(binary>0,1,0)] %>%
                        melt(.,id.vars=c("psu","model_sampling_unit","gear_type","species_cd","binary")) %>%
                        setnames(.,"variable","category") %>%
                        .[,breaks:=value] %>%
                        .[,empirical:=mean(binary),by=.(gear_type,species_cd,category,breaks)] %>%
                        .[,.(gear_type,species_cd,category,breaks,empirical)] %>%
                        unique(.) %>%
                        .[,y_max_a:=max(empirical),by=.(category,species_cd,gear_type)] %>%
                        .[,empirical_sc_a:=empirical/y_max_a] %>%
                        .[,y_max_b:=max(empirical),by=.(category,species_cd)] %>%
                        .[,empirical_sc_b:=empirical/y_max_b] %>%
                        .[,breaks:=factor(breaks,levels=c("Niihau","Kauai","Oahu","Maui Nui", "Big Island"))] %>%
                        .[,species_cd:=factor(species_cd,levels=deep7_code_vec,labels=deep7_name_vec)]

#_____________________________________________________________________________________________________________________________
# define helper functions

  extract_continuous_abundance_effect = function(fit_internal,
                                                model_name,
                                                continuous_ab_variables = c("depth_sc","complexity_sc","hardness_sc","slope_sc"),
                                                target_species = c("prfi","etca","etco","prsi","przo","hyqu","apru"))
  {
    # checks
      continuous_ab_variables = intersect(continuous_ab_variables,colnames(fit_internal$covariate_data))
      keep_vec = rep(NA,length(continuous_ab_variables))
      for(i in 1:length(continuous_ab_variables))
      {
        keep_vec[i] =ifelse(sum(grep(continuous_ab_variables[i],dimnames(fit_internal$data_list$X1_gctp)[[4]]))>0,TRUE,FALSE)
      }
      continuous_ab_variables = continuous_ab_variables[keep_vec]

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
      ab_df_plot = ab_df_plot[,-which(colnames(ab_df_plot)%in%continuous_ab_variables)]
      colnames(ab_df_plot)[which(colnames(ab_df_plot)%in%gsub("_sc","",continuous_ab_variables))] = paste0(colnames(ab_df_plot)[which(colnames(ab_df_plot)%in%gsub("_sc","",continuous_ab_variables))],"_sc")
            
      abundance_effect_dt.list = as.list(rep(NA,length(continuous_ab_variables)))
      for(i in 1:length(continuous_ab_variables))
      {
        abundance_effect_dt.list[[i]] = as.list(rep(NA,length(target_species)))
        tmp_cols = continuous_ab_variables[i]
        for(s in 1:length(target_species))
        {
          abundance_effect_dt.list[[i]][[s]] = as.data.table(ab_df_plot) %>%
            .[,..tmp_cols] %>%
            .[,value:=fit_internal$data_list$X1_gctp[,s,1,grep(continuous_ab_variables[i],dimnames(fit_internal$data_list$X1_gctp)[[4]])] %*% t(t(fit_internal$ParHat$gamma1_cp[s,grep(continuous_ab_variables[i],dimnames(fit_internal$data_list$X1_gctp)[[4]])]))] %>%
            .[,species_cd:=target_species[s]] %>%
            setnames(.,tmp_cols,"variable") %>%
            .[,value_link:=exp(value)] %>%
            .[,category:=tmp_cols] %>%
            .[,component:="1st"] %>%
            .[,value_link_rescale:=value_link/max(value_link)] %>%
            .[,value_rescale:=rescale(value)]
        }
        abundance_effect_dt.list[[i]] = rbindlist(abundance_effect_dt.list[[i]])
        rm(list="tmp_cols")
      }
      abundance_effect_dt = rbindlist(abundance_effect_dt.list) %>%
                            .[,model_name:=model_name] %>%
                            .[,.(model_name,component,species_cd,category,variable,value,value_link,value_link_rescale,value_rescale)] %>%
                            .[,category:=gsub("_sc","",category)]
     
    }
    return(abundance_effect_dt)  
                    
  }


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


#_____________________________________________________________________________________________________________________________
# extract & plot fitted continuous abundance effects
  # only extract those that have continuous abundance covariates
  continuous_abundance_models = ests_vec[-grep("_island",ests_vec,fixed=TRUE)]
  continuous_abundance_models = continuous_abundance_models[-grep("_v_TRUE",continuous_abundance_models,fixed=TRUE)]
  # order by last modified
  continuous_abundance_models = continuous_abundance_models[order(sapply(paste0(proj_dir,"VAST/model_runs/",continuous_abundance_models),file.mtime))]

  continuous_ab_effect_dt.list = as.list(rep(NA,length(continuous_abundance_models)))
  
  for(i in seq_along(continuous_ab_effect_dt.list))
  {
    load(paste0(proj_dir,"VAST/model_runs/",continuous_abundance_models[i]))
    tmp_model = strsplit(strsplit(continuous_abundance_models[i],"/")[[1]][2],"_")[[1]]
    tmp_model_name = paste0(tmp_model[1],".",tmp_model[4],".",
                            tmp_model[5],".",
                            paste0(sapply(strsplit(tmp_model[6],"[.]"),function(x)paste0(substr(x,1,1),substr(x,nchar(x),nchar(x)))),collapse=""))
    continuous_ab_effect_dt.list[[i]] = extract_continuous_abundance_effect(fit_internal = fit,model_name=tmp_model_name) %>%
      .[,model_name_order:=i]
    rm(list=c("fit","tmp_model","tmp_model_name"))
  }

  continuous_ab_effect_dt = rbindlist(continuous_ab_effect_dt.list) %>%
                            unique(.) %>%
                            na.omit(.)
  continuous_ab_effect_dt$model_name = factor(continuous_ab_effect_dt$model_name,levels=unique(continuous_ab_effect_dt$ model_name))
  continuous_ab_effect_dt$species_cd=factor(continuous_ab_effect_dt$species_cd,levels=deep7_code_vec,labels=deep7_name_vec)

    p = plot_empirical_dt %>%
    .[gear_type=="camera"] %>%
    ggplot() +
    ylab("Relative encounter rate") +
	  xlab("Habitat covariate") +
	  facet_wrap(category~species_cd,scales="free",ncol=7) +
	  geom_rect(aes(xmin=breaks-x_width,xmax=breaks,ymin=0,ymax=empirical_sc_a,fill=gear_type),color="white") +
 	  geom_rect(data=plot_empirical_dt[gear_type=="research_fishing"],aes(xmin=breaks,xmax=breaks+x_width,ymin=0,ymax=empirical_sc_a,fill=gear_type),color="white") +
    geom_line(data=continuous_ab_effect_dt,aes(x=variable,y=value_link_rescale,color=model_name,group=model_name)) +
    geom_hline(yintercept=0) +
    theme_few(base_size=20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.ticks.y=element_blank(),axis.text.y=element_blank(),strip.text.x = element_text(size = 10)) +
    scale_fill_manual("Gear type",values=c("gray60","gray80")) +
	  viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
    ggsave(filename=paste0("continuous_abundance_effects_a.png"), plot = p, device = "png", path = paste0(proj_dir,"VAST/model_runs/comparison_plots/"),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

    p = plot_empirical_dt %>%
    .[gear_type=="camera"] %>%
    ggplot() +
    ylab("Relative encounter rate") +
	  xlab("Habitat covariate") +
	  facet_wrap(category~species_cd,scales="free",ncol=7) +
	  geom_rect(aes(xmin=breaks-x_width,xmax=breaks,ymin=0,ymax=empirical_sc_b,fill=gear_type),color="white") +
 	  geom_rect(data=plot_empirical_dt[gear_type=="research_fishing"],aes(xmin=breaks,xmax=breaks+x_width,ymin=0,ymax=empirical_sc_b,fill=gear_type),color="white") +
    geom_line(data=continuous_ab_effect_dt,aes(x=variable,y=value_link_rescale,color=model_name,group=model_name)) +
    geom_hline(yintercept=0) +
    theme_few(base_size=20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.ticks.y=element_blank(),axis.text.y=element_blank(),strip.text.x = element_text(size = 10)) +
    scale_fill_manual("Gear type",values=c("gray60","gray80")) +
	  viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
    ggsave(filename=paste0("continuous_abundance_effects_b.png"), plot = p, device = "png", path = paste0(proj_dir,"VAST/model_runs/comparison_plots/"),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
#_____________________________________________________________________________________________________________________________
# extract & plot fitted discrete abundance effects
  # only extract those that have continuous abundance covariates
  discrete_abundance_models = ests_vec[grep("island",ests_vec,fixed=TRUE)]
  # order by last modified
  discrete_abundance_models = discrete_abundance_models[order(sapply(paste0(proj_dir,"VAST/model_runs/",discrete_abundance_models),file.mtime))]

  discrete_ab_effect_dt.list = as.list(rep(NA,length(discrete_abundance_models)))
  
  for(i in seq_along(discrete_ab_effect_dt.list))
  {
    load(paste0(proj_dir,"VAST/model_runs/",discrete_abundance_models[i]))
    tmp_model = strsplit(strsplit(discrete_abundance_models[i],"/")[[1]][2],"_")[[1]]
    tmp_model_name = paste0(tmp_model[1],".",tmp_model[4],".",
                            tmp_model[5],".",
                            paste0(sapply(strsplit(tmp_model[6],"[.]"),function(x)paste0(substr(x,1,1),substr(x,nchar(x),nchar(x)))),collapse=""))
    discrete_ab_effect_dt.list[[i]] = extract_discrete_abundance_effect(fit_internal = fit,model_name=tmp_model_name) %>%
      .[,model_name_order:=i]
    rm(list=c("fit","tmp_model","tmp_model_name"))
  }

  discrete_ab_effect_dt = rbindlist(discrete_ab_effect_dt.list) %>%
                            unique(.) %>%
                            na.omit(.)
  discrete_ab_effect_dt$model_name = factor(discrete_ab_effect_dt$model_name,levels=unique(discrete_ab_effect_dt$ model_name))
  discrete_ab_effect_dt$variable = factor(discrete_ab_effect_dt$variable,levels=c("Niihau","Kauai","Oahu","Maui Nui", "Big Island"))
  discrete_ab_effect_dt$species_cd=factor(discrete_ab_effect_dt$species_cd,levels=deep7_code_vec,labels=deep7_name_vec)
  
    
    p = plot_empirical_discrete_dt %>%
    ggplot() +
    ylab("Relative encounter rate") +
	  xlab("Island group") +
	  facet_wrap(category~species_cd,scales="free",ncol=7) +
	  geom_bar(aes(x=breaks,y=empirical_sc_a,fill=gear_type),stat="identity",color="white",position=position_dodge()) +
    geom_point(data=discrete_ab_effect_dt,aes(x=variable,y=value_link_rescale,color=model_name,group=model_name),position="jitter",size=2) +
    geom_hline(yintercept=0) +
    theme_few(base_size=20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.ticks.y=element_blank(),axis.text.y=element_blank(),strip.text.x = element_text(size = 10)) +
    scale_fill_manual("Gear type",values=c("gray60","gray80")) +
	  viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
    ggsave(filename=paste0("discrete_abundance_effects_a.png"), plot = p, device = "png", path = paste0(proj_dir,"VAST/model_runs/comparison_plots/"),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

    p = plot_empirical_discrete_dt %>%
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
    ggsave(filename=paste0("discrete_abundance_effects_b.png"), plot = p, device = "png", path = paste0(proj_dir,"VAST/model_runs/comparison_plots/"),
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
          