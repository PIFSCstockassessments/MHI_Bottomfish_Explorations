

# Nicholas Ducharme-Barth
# 2023/06/04
# gather shiny data

# Copyright (c) 2023 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(ggplot2)
	library(ggthemes)


#_____________________________________________________________________________________________________________________________
# define helper function	
	read_vast_parameter_estimates = function(tmp_par)
	{	
		target_var = c("time_for_run","objective","number_of_coefficients","AIC","max_gradient","diagnostics")
		start_index = grep('^\\$',tmp_par)+1
		end_index = c(start_index[-1]-3,length(tmp_par)-1)
		variable_name = gsub('\\$','',tmp_par[start_index-1])
		
		# subset to important quantities
		start_index = start_index[which(variable_name %in% target_var)]
		end_index = end_index[which(variable_name %in% target_var)]
		variable_name = variable_name[which(variable_name %in% target_var)]
		
		# define storage
		out_list = as.list(rep(NA,length(variable_name)))
		names(out_list) = variable_name
		
		out_list$time_for_run = as.numeric(gsub("Time difference of ","",gsub(" secs","",tmp_par[start_index[which(variable_name=="time_for_run")]],fixed=TRUE),fixed=TRUE))
		out_list$objective = as.numeric(gsub("[1] ","",tmp_par[start_index[which(variable_name=="objective")]],fixed=TRUE))
		out_list$number_of_coefficients = as.numeric(strsplit(trimws(tmp_par[start_index[which(variable_name=="number_of_coefficients")]+1]),"\\s+")[[1]])
		names(out_list$number_of_coefficients) = strsplit(trimws(tmp_par[start_index[which(variable_name=="number_of_coefficients")]]),"\\s+")[[1]]
		out_list$AIC = as.numeric(gsub("[1] ","",tmp_par[start_index[which(variable_name=="AIC")]],fixed=TRUE))
		out_list$max_gradient = as.numeric(gsub("[1] ","",tmp_par[start_index[which(variable_name=="max_gradient")]],fixed=TRUE))

		out_list$diagnostics = as.data.frame(t(unname(sapply(tmp_par[(start_index[which(variable_name=="diagnostics")]+1):end_index[which(variable_name=="diagnostics")]],function(x)strsplit(trimws(x),"\\s+")[[1]][-1]))))
		colnames(out_list$diagnostics) = strsplit(trimws(tmp_par[start_index[which(variable_name=="diagnostics")]]),"\\s+")[[1]]
		out_list$diagnostics$starting_value = as.numeric(out_list$diagnostics$starting_value)
		out_list$diagnostics$Lower = as.numeric(out_list$diagnostics$Lower)
		out_list$diagnostics$MLE = as.numeric(out_list$diagnostics$MLE)
		out_list$diagnostics$Upper = as.numeric(out_list$diagnostics$Upper)
		out_list$diagnostics$final_gradient = as.numeric(out_list$diagnostics$final_gradient)
		return(out_list)
	}

  extract_abundance_effects = function(fit_internal,
                                                model_name,
                                                continuous_variables = c("depth","complexity","hardness","slope"),
												discrete_variables = c("island"),
                                                target_species = c("prfi","etca","etco","prsi","przo","hyqu","apru"))
  {
	# checks
      variables = colnames(fit_internal$covariate_data)[as.vector(unlist(sapply(c(discrete_variables,continuous_variables),function(x)grep(x,colnames(fit_internal$covariate_data),fixed=TRUE))))]
      keep_vec_1 = rep(NA,length(variables))
	  process_formula = strsplit(trimws(gsub("\\d"," ",paste0(gsub("[\\(|\\~|\\=|\\)|\\+|\\,]", " ", fit_internal$X1_formula),collapse=""))),"\\s+")[[1]]
      for(i in 1:length(variables))
      {
        keep_vec_1[i] =ifelse(sum(grep(paste0("^",variables[i],"$"),process_formula))>0,TRUE,FALSE)
      }

      variables = variables[keep_vec_1]
	  keep_column_names = variables
	  keep_column_names[grep("_sc",variables,fixed=TRUE)] = gsub("_sc","",keep_column_names[grep("_sc",variables,fixed=TRUE)],fixed=TRUE)

    if(fit_internal$ParHat[1] == "Model is not converged")
    {
      abundance_effect_dt = data.table(model_name=model_name,
                                       component=NA,
                                       species_cd=NA,
                                       category=NA,
                                       variable_c=NA,
									   variable_d=NA,
                                       value=NA,
                                       value_link=NA,
                                       value_link_rescale=NA,
                                       value_rescale=NA)
    } else {
      ab_df_plot = fit_internal$covariate_data
	  if(length(keep_column_names)==1)
	  {
		ab_df_plot = data.frame(t=ab_df_plot[,keep_column_names])
	  } else {
		ab_df_plot = ab_df_plot[,keep_column_names]
	  }
      colnames(ab_df_plot) = variables

            
      abundance_effect_dt.list = as.list(rep(NA,length(variables)))
      for(i in 1:length(variables))
      {
        abundance_effect_dt.list[[i]] = as.list(rep(NA,length(target_species)))
        tmp_cols = variables[i]
        for(s in 1:length(target_species))
        {
          abundance_effect_dt.list[[i]][[s]] = as.data.table(ab_df_plot) %>%
            .[,..tmp_cols] %>%
            .[,value:=fit_internal$data_list$X1_gctp[,s,1,grep(variables[i],dimnames(fit_internal$data_list$X1_gctp)[[4]])] %*% t(t(fit_internal$ParHat$gamma1_cp[s,grep(variables[i],dimnames(fit_internal$data_list$X1_gctp)[[4]])]))] %>%
            .[,species_cd:=target_species[s]]
		  if(keep_column_names[i] %in% continuous_variables)
		  {
			abundance_effect_dt.list[[i]][[s]] = abundance_effect_dt.list[[i]][[s]] %>%
				setnames(.,tmp_cols,"variable_c") %>%
				.[,variable_d:=NA]
		  } else {
						abundance_effect_dt.list[[i]][[s]] = abundance_effect_dt.list[[i]][[s]] %>%
				setnames(.,tmp_cols,"variable_d") %>%
				.[,variable_c:=NA]
		  }
            abundance_effect_dt.list[[i]][[s]] = abundance_effect_dt.list[[i]][[s]] %>%
            .[,value_link:=exp(value)] %>%
            .[,category:=tmp_cols] %>%
            .[,component:="1st"] %>%
            .[,.(component,species_cd,category,variable_c,variable_d,value,value_link)]

        }
        abundance_effect_dt.list[[i]] = rbindlist(abundance_effect_dt.list[[i]])
        rm(list="tmp_cols")
      }
      abundance_effect_dt = rbindlist(abundance_effect_dt.list) %>%
                            .[,model_name:=model_name] %>%
							.[,value_link_rescale:=value_link/max(value_link),by=.(species_cd,category)] %>%
            				.[,value_rescale:=scales::rescale(value),by=.(species_cd,category)] %>%
                            .[,.(model_name,component,species_cd,category,variable_c,variable_d,value,value_link,value_link_rescale,value_rescale)] %>%
                            .[,category:=gsub("_sc","",category)] %>%
							unique(.)
     
    }
    return(abundance_effect_dt)  
                    
  }

  
  extract_catchability_effect = function(fit_internal,
                                                model_name,
                                                discrete_variables = c("gear_type","platform"),
												continuous_variables =c("time","month","lunarphase"))
  {
    # checks
      variables = intersect(c(discrete_variables,continuous_variables),colnames(fit_internal$catchability_data))
      keep_vec_1 = rep(NA,length(variables))
	  keep_vec_2 = rep(NA,length(variables))
      for(i in 1:length(variables))
      {
        keep_vec_1[i] =ifelse(sum(grep(variables[i],fit_internal$Q1_formula))>0,TRUE,FALSE)
        keep_vec_2[i] =ifelse(sum(grep(variables[i],fit_internal$Q2_formula))>0,TRUE,FALSE)
      }

      variables_1 = variables[keep_vec_1]
	  variables_2 = variables[keep_vec_2]

    if(fit_internal$ParHat[1] == "Model is not converged")
    {
      q_effect_dt = data.table(model_name=model_name,
                                       component=NA,
                                       species_cd=NA,
                                       category=NA,
                                       variable_c=NA,
									   variable_d=NA,
                                       value=NA,
                                       value_link=NA,
                                       value_link_rescale=NA,
                                       value_rescale=NA)
    } else {
	  # fit_internal$catchability_data$a_i = fit_internal$data_list$a_i
      q_df_1 = fit_internal$catchability_data[,c("year","category",variables_1)]
      q_df_2 = fit_internal$catchability_data[,c("year","category",variables_2)]

      q_effect_dt_1.list = as.list(rep(NA,length(variables_1)))
      for(i in 1:length(variables_1))
      {
		  tmp_cols = c("year","category",variables_1[i])
          q_effect_dt_1.list[[i]] = as.data.table(q_df_1) %>%
            .[,..tmp_cols] %>%
            .[,value:=fit_internal$data_list$Q1_ik[,grep(variables_1[i],dimnames(fit_internal$data_list$Q1_ik)[[2]])] %*% t(t(fit_internal$ParHat$lambda1_k[grep(variables_1[i],dimnames(fit_internal$data_list$Q1_ik)[[2]])]))] %>%
            setnames(.,c("category",variables_1[i]),c("species_cd","variable"))
		  if(variables_1[i]%in%continuous_variables)
		  {
			q_effect_dt_1.list[[i]] = q_effect_dt_1.list[[i]] %>%
				setnames(.,"variable","variable_c") %>%
				.[,variable_d:=NA]
		  } else {
			q_effect_dt_1.list[[i]] = q_effect_dt_1.list[[i]] %>%
				setnames(.,"variable","variable_d") %>%
				.[,variable_c:=NA]
		  }
		  q_effect_dt_1.list[[i]] = q_effect_dt_1.list[[i]] %>%
            .[,value_link:=exp(value)] %>%
            .[,category:=variables_1[i]] %>%
            .[,component:="1st"] %>%
			.[,.(model_name,component,species_cd,category,variable_c,variable_d,value,value_link)]
        
        rm(list="tmp_cols")
      }
      q_effect_dt_1 = rbindlist(q_effect_dt_1.list) %>%
                            .[,model_name:=model_name] %>%
                            .[,value_link_rescale:=value_link/max(value_link),by=.(species_cd,category)] %>%
                            .[,value_rescale:=scales::rescale(value),by=.(species_cd,category)] %>%
                            .[,.(model_name,component,species_cd,category,variable_c,variable_d,value,value_link,value_link_rescale,value_rescale)]     
	  if(length(variables_2)>0)
	  {
		q_effect_dt_2.list = as.list(rep(NA,length(variables_2)))
		for(i in 1:length(variables_2))
		{
			tmp_cols = c("year","category",variables_2[i])
			q_effect_dt_2.list[[i]] = as.data.table(q_df_2) %>%
				.[,..tmp_cols] %>%
				.[,value:=fit_internal$data_list$Q2_ik[,grep(variables_2[i],dimnames(fit_internal$data_list$Q2_ik)[[2]])] %*% t(t(fit_internal$ParHat$lambda2_k[grep(variables_2[i],dimnames(fit_internal$data_list$Q2_ik)[[2]])]))] %>%
				setnames(.,c("category",variables_2[i]),c("species_cd","variable"))
		  if(variables_2[i]%in%continuous_variables)
		  {
			q_effect_dt_2.list[[i]] = q_effect_dt_2.list[[i]] %>%
				setnames(.,"variable","variable_c") %>%
				.[,variable_d:=NA]
		  } else {
			q_effect_dt_2.list[[i]] = q_effect_dt_2.list[[i]] %>%
				setnames(.,"variable","variable_d") %>%
				.[,variable_c:=NA]
		  }
		  q_effect_dt_2.list[[i]] = q_effect_dt_2.list[[i]] %>%
				.[,value_link:=exp(value)] %>%
				.[,category:=variables_2[i]] %>%
				.[,component:="2nd"] %>%
				.[,.(model_name,component,species_cd,category,variable_c,variable_d,value,value_link)]
			
			rm(list="tmp_cols")
		}
		q_effect_dt_2 = rbindlist(q_effect_dt_2.list) %>%
								.[,model_name:=model_name] %>%
								.[,value_link_rescale:=value_link/max(value_link),by=.(species_cd,category)] %>%
								.[,value_rescale:=scales::rescale(value),by=.(species_cd,category)] %>%
								.[,.(model_name,component,species_cd,category,variable_c,variable_d,value,value_link,value_link_rescale,value_rescale)]     
		
	  } else {
		      q_effect_dt_2 = data.table(model_name=model_name,
                                       component="2nd",
                                       species_cd=NA,
                                       category=NA,
                                       variable_c=NA,
									   variable_d=NA,
                                       value=NA,
                                       value_link=NA,
                                       value_link_rescale=NA,
                                       value_rescale=NA)
	  }

	  q_effect_dt = rbind(q_effect_dt_1,q_effect_dt_2) %>% unique(.)
    }
    return(unique(q_effect_dt))               
  }


#_____________________________________________________________________________________________________________________________
# set working directory
	proj_dir = "D:/HOME/SAP/2024_Deep7/"

    deep7_code_vec = tolower(c("ETCA","APRU","PRSI","HYQU","PRFI","PRZO","ETCO"))
	deep7_name_vec = c("Ehu (ETCA)", "Lehi (APRU)", "Kalekale (PRSI)", "Hapu'upu'u (HYQU)", "'Opakapaka (PRFI)", "Gindai (PRZO)", "Onaga (ETCO)")

	dir_vec = dir(paste0(proj_dir,"VAST/model_runs/"),recursive=TRUE)
	ests_vec = dir_vec[grep("/index_dt.csv",dir_vec,fixed=TRUE)]
   	complete_vec = dir_vec[grep("/setup/parameter_estimates.txt",dir_vec,fixed=TRUE)]

	ests_vec_stem = unname(sapply(ests_vec,function(x)paste0(strsplit(x,"[/]")[[1]][1:2],collapse="/")))
	complete_vec_stem = unname(sapply(complete_vec,function(x)paste0(strsplit(x,"[/]")[[1]][1:2],collapse="/")))

	# ignore write-out file
    if("comparison_plots/index_dt.csv" %in% ests_vec)
    {
        ests_vec = ests_vec[-which(ests_vec == "comparison_plots/index_dt.csv")]
		ests_vec_stem = ests_vec_stem[-which(ests_vec_stem == "comparison_plots/index_dt.csv")]
    }

	ests_vec = ests_vec[which(ests_vec_stem %in% complete_vec_stem)]


#_____________________________________________________________________________________________________________________________
# load psu data
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
# bring in indices
    # design based bfish index data is in kg so convert to lbs (millions)
	design_dt = fread(paste0(proj_dir,"Data/2022_design_based_estimates.csv")) %>%
				.[,.(Model,Category,Time,Estimate,CV)] %>%
				 setnames(.,c("Model","Category","Time","Estimate","CV"),c("model_name_short","category","time","estimate","cv")) %>%
				.[model_name_short=="design",model_name_short:="00 design"] %>%
				.[,name:="design"] %>%
				.[,date:=as.POSIXct(NA)] %>%
				.[,data_year:=2022] %>%		  
				.[,error_structure:="design"] %>%
				.[,species := "mv"] %>%
				.[,data_treatment := "design"] %>%
				.[,q_config := "design"] %>%
				.[,ab_config := "design"] %>%
				.[,lehi_filter := FALSE] %>%
				.[,knot_dist:=NA] %>%
				.[,fine_scale := NA] %>%
				.[,all_lengths:="FALSE"] %>%
				.[,gears_used := "both"] %>%
                .[,estimate:=estimate*2.20462262185] %>%
				.[,estimate:=estimate/1000000] %>%
				.[,l95:=exp(log(estimate)-2*sqrt(log(cv^2+1)))] %>%
				.[,u95:=exp(log(estimate)+2*sqrt(log(cv^2+1)))] %>%
                .[,category:=deep7_name_vec[match(category,deep7_code_vec)]] %>%
                .[is.na(category),category:="Total"] %>%
                .[,category:=factor(category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))] %>%
			    .[,estimate_sc:=estimate/mean(estimate),by=.(category)] %>%
			    .[,l95_sc:=exp(log(estimate_sc)-2*sqrt(log(cv^2+1)))] %>%
			    .[,u95_sc:=exp(log(estimate_sc)+2*sqrt(log(cv^2+1)))] %>%
				.[,.(name,date,data_year,error_structure,species,data_treatment,q_config,ab_config,lehi_filter,knot_dist,fine_scale,all_lengths,gears_used,model_name_short,category,time,estimate,cv,l95,u95,estimate_sc,l95_sc,u95_sc)]


#_____________________________________________________________________________________________________________________________
# extract other indices
  # order by last modified
  ests_vec = ests_vec[order(sapply(paste0(proj_dir,"VAST/model_runs/",ests_vec),file.mtime))]

  q_effect_dt.list = abundance_effect_dt.list = index_summary_dt.list = index_dt.list = as.list(rep(NA,length(ests_vec)))
  
  for(i in seq_along(index_dt.list))
  {
	tmp_date = strsplit(ests_vec[i],"/")[[1]][1]
	tmp_name = strsplit(ests_vec[i],"/")[[1]][2]
    tmp_model = strsplit(strsplit(ests_vec[i],"/")[[1]][2],"_")[[1]]
    tmp_model_name = paste0(tmp_model[1],".",tmp_model[4],".",
                            tmp_model[5],".",
                            paste0(sapply(strsplit(tmp_model[6],"[.]"),function(x)paste0(substr(x,1,1),substr(x,nchar(x),nchar(x)))),collapse=""))
 	index_dt.list[[i]] = fread(paste0(proj_dir,"VAST/model_runs/",ests_vec[i])) %>%
					  	.[,name:=tmp_name] %>%
						.[,date:=tmp_date] %>%
						.[,data_year:=strsplit(name,"_")[[1]][1]] %>%		  
						.[,error_structure:=strsplit(name,"_")[[1]][2]] %>%
						.[,species := strsplit(name,"_")[[1]][3]] %>%
						.[,data_treatment := strsplit(name,"_")[[1]][4]] %>%
						.[,q_config := strsplit(name,"_")[[1]][5]] %>%
						.[,ab_config := strsplit(name,"_")[[1]][6]] %>%
						.[,lehi_filter := strsplit(name,"_")[[1]][7]] %>%
						.[,knot_dist := as.numeric(strsplit(name,"_")[[1]][8])] %>%
						.[,fine_scale := strsplit(name,"_")[[1]][9]] %>%
						.[,all_lengths:="FALSE"] %>%
						.[,gears_used := "both"] %>%
						.[,date:=as.POSIXct(date)] %>%
						 setnames(.,c("Model","Category","Time","Estimate","CV"),c("model_name_short","category","time","estimate","cv")) %>%
                     .[,model_name_short:=paste0(ifelse(i<10,"0",""),i," ",tmp_model_name)] %>%
                     .[,category:=deep7_name_vec[match(category,deep7_code_vec)]] %>%
                     .[is.na(category),category:="Total"] %>%
                     .[,category:=factor(category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))] %>%
					 .[,estimate_sc:=estimate/mean(estimate),by=.(category)] %>%
			         .[,l95_sc:=exp(log(estimate_sc)-2*sqrt(log(cv^2+1)))] %>%
			         .[,u95_sc:=exp(log(estimate_sc)+2*sqrt(log(cv^2+1)))] %>%
					 .[,.(name,date,data_year,error_structure,species,data_treatment,q_config,ab_config,lehi_filter,knot_dist,fine_scale,all_lengths,gears_used,model_name_short,category,time,estimate,cv,l95,u95,estimate_sc,l95_sc,u95_sc)]
		
		if(length(strsplit(tmp_name,"_")[[1]])>=14)
		{
			index_dt.list[[i]]$all_lengths = strsplit(tmp_name,"_")[[1]][13]
			index_dt.list[[i]]$gears_used = strsplit(tmp_name,"_")[[1]][14]
		} 
		parameter_estimates = read_vast_parameter_estimates(readLines(paste0(proj_dir,"VAST/model_runs/",tmp_date,"/",tmp_name,"/setup/parameter_estimates.txt")))
		parameter_estimates_final = read_vast_parameter_estimates(readLines(paste0(proj_dir,"VAST/model_runs/",tmp_date,"/",tmp_name,"/parameter_estimates.txt")))
		var_par = grep("\\bL_",parameter_estimates_final$diagnostics$Param)
		problem_par = sum(which(abs(parameter_estimates_final$diagnostics[var_par,"MLE"])<1e-3|abs(parameter_estimates_final$diagnostics[var_par,"MLE"])>1.5e1))

		# Add check to see if residuals are calculated, then make a summary data.table
				index_summary_dt.list[[i]] = data.table(date=tmp_date,
										  name=tmp_name) %>%
  		                        .[,model_name_short:=paste0(ifelse(i<10,"0",""),i," ",tmp_model_name)] %>%
								.[,data_year:=strsplit(name,"_")[[1]][1]] %>%		  
								.[,error_structure:=strsplit(name,"_")[[1]][2]] %>%
								.[,species := strsplit(name,"_")[[1]][3]] %>%
								.[,data_treatment := strsplit(name,"_")[[1]][4]] %>%
								.[,q_config := strsplit(name,"_")[[1]][5]] %>%
								.[,ab_config := strsplit(name,"_")[[1]][6]] %>%
								.[,lehi_filter := strsplit(name,"_")[[1]][7]] %>%
								.[,knot_dist := as.numeric(strsplit(name,"_")[[1]][8])] %>%
								.[,fine_scale := strsplit(name,"_")[[1]][9]] %>%
								.[,all_lengths:="FALSE"] %>%
								.[,gears_used := "both"] %>%
								.[,date:=as.POSIXct(date)] %>%
								.[,runtime:=parameter_estimates$time_for_run] %>%
								.[,nll:=parameter_estimates$objective] %>%
								.[,n_par:=parameter_estimates$number_of_coefficients[1]] %>%
								.[,n_fixed:=parameter_estimates$number_of_coefficients[2]] %>%
								.[,n_random:=parameter_estimates$number_of_coefficients[3]] %>%
								.[,aic:=parameter_estimates$AIC] %>%
								.[,aic_all:=2*n_par+2*nll] %>%
								.[,mgc:=parameter_estimates$max_gradient] %>%
								# .[,pdh:=parameter_estimates$SD$pdHess] %>%
								.[,bad_param:=problem_par] %>%
								cbind(.,t(t(table(parameter_estimates$diagnostics$Param)))) %>%
								.[,V2:=NULL] %>%
								setnames(.,"V1","parameter")
		if(length(strsplit(tmp_name,"_")[[1]])>=14)
		{
			index_summary_dt.list[[i]]$all_lengths = strsplit(tmp_name,"_")[[1]][13]
			index_summary_dt.list[[i]]$gears_used = strsplit(tmp_name,"_")[[1]][14]
		} 
		if(tmp_model[6]!="v"|tmp_model[5]!="v")
		{
			load(paste0(proj_dir,"VAST/model_runs/",tmp_date,"/",tmp_name,"/fit.RData"))
		}						
		if(tmp_model[6]!="v")
		{
			if(strsplit(tmp_name,"_")[[1]][3]=="mv")
			{
				abundance_effect_dt.list[[i]] = extract_abundance_effects(fit_internal=fit,model_name=paste0(ifelse(i<10,"0",""),i," ",tmp_model_name))
			} else {
				abundance_effect_dt.list[[i]] = extract_abundance_effects(fit_internal=fit,model_name=paste0(ifelse(i<10,"0",""),i," ",tmp_model_name),target_species=strsplit(tmp_name,"_")[[1]][3])
			}
			
		} else {
			abundance_effect_dt.list[[i]] = data.table(model_name=paste0(ifelse(i<10,"0",""),i," ",tmp_model_name),
                                       component=NA,
                                       species_cd=NA,
                                       category=NA,
                                       variable_c=NA,
									   variable_d=NA,
                                       value=NA,
                                       value_link=NA,
                                       value_link_rescale=NA,
                                       value_rescale=NA)
		}
		if(tmp_model[5]!="v")
		{
			q_effect_dt.list[[i]] = extract_catchability_effect(fit_internal=fit,model_name=paste0(ifelse(i<10,"0",""),i," ",tmp_model_name))
		} else {
			q_effect_dt.list[[i]] = data.table(model_name=paste0(ifelse(i<10,"0",""),i," ",tmp_model_name),
                                       component=NA,
                                       species_cd=NA,
                                       category=NA,
                                       variable_c=NA,
									   variable_d=NA,
                                       value=NA,
                                       value_link=NA,
                                       value_link_rescale=NA,
                                       value_rescale=NA)
		}

    rm(list=c("fit","tmp_model","tmp_model_name","tmp_name","tmp_date","parameter_estimates_final","var_par","problem_par","parameter_estimates"))
  }

  index_dt = rbind(design_dt,rbindlist(index_dt.list)) %>%
			.[,model_number:=as.numeric(sapply(model_name_short,function(x)strsplit(x,"\\s+")[[1]][1]))] %>%
			.[order(model_number,time,category)]

  index_summary_dt = rbindlist(index_summary_dt.list) %>%
				 dcast(.,date+name+model_name_short+data_year+error_structure+species+data_treatment+q_config+ab_config+lehi_filter+knot_dist+fine_scale+all_lengths+gears_used+runtime+nll+n_par+n_fixed+n_random+aic+aic_all+mgc+bad_param~parameter) %>%
					merge(.,data.table(model_name_short="00 design",q_config="design",ab_config="design"),by=c("model_name_short","q_config","ab_config"),all=TRUE) %>%
  					.[,model_number:=as.numeric(sapply(model_name_short,function(x)strsplit(x,"\\s+")[[1]][1]))] %>%
					.[order(model_number)] %>%
					.[,.(model_number,model_name_short,q_config,ab_config,data_year,error_structure,species,data_treatment,lehi_filter,knot_dist,fine_scale,all_lengths,gears_used,runtime,nll,n_par,n_fixed,n_random,aic,aic_all,mgc,bad_param,L_epsilon1_z,L_epsilon2_z,L_eta1_z,L_eta2_z,L_omega1_z,L_omega2_z,beta1_ft,beta2_ft,gamma1_cp,lambda1_k,lambda2_k,logSigmaM,logkappa1,logkappa2)]

  abundance_effect_dt = rbindlist(abundance_effect_dt.list) %>%
					setnames(.,"model_name","model_name_short") %>%
					merge(.,data.table(model_name_short="00 design"),by=c("model_name_short"),all=TRUE) %>%
  					.[,model_number:=as.numeric(sapply(model_name_short,function(x)strsplit(x,"\\s+")[[1]][1]))] %>%
					.[order(model_number)] %>%
					.[,.(model_number,model_name_short,component,species_cd,category,variable_c,variable_d,value,value_link,value_link_rescale,value_rescale)] %>%
					.[,species_cd:=deep7_name_vec[match(species_cd,deep7_code_vec)]] %>%
                    .[,species_cd:=factor(species_cd,levels=c("'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

  q_effect_dt = rbindlist(q_effect_dt.list) %>%
					setnames(.,"model_name","model_name_short") %>%
					merge(.,data.table(model_name_short="00 design"),by=c("model_name_short"),all=TRUE) %>%
  					.[,model_number:=as.numeric(sapply(model_name_short,function(x)strsplit(x,"\\s+")[[1]][1]))] %>%
					.[order(model_number)] %>%
					.[,.(model_number,model_name_short,component,species_cd,category,variable_c,variable_d,value,value_link,value_link_rescale,value_rescale)] %>%
					.[,species_cd:=deep7_name_vec[match(species_cd,deep7_code_vec)]] %>%
                    .[,species_cd:=factor(species_cd,levels=c("'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]


#_____________________________________________________________________________________________________________________________
# save
	fwrite(index_dt,file=paste0(proj_dir,"VAST/model_runs/comparison_plots/index_dt.csv"))
	fwrite(index_summary_dt,file=paste0(proj_dir,"VAST/model_runs/comparison_plots/index_summary_dt.csv"))
	fwrite(abundance_effect_dt,file=paste0(proj_dir,"VAST/model_runs/comparison_plots/abundance_effect_dt.csv"))
	fwrite(q_effect_dt,file=paste0(proj_dir,"VAST/model_runs/comparison_plots/q_effect_dt.csv"))

	fwrite(plot_empirical_dt,file=paste0(proj_dir,"VAST/model_runs/comparison_plots/plot_empirical_ab_dt.csv"))
	fwrite(plot_empirical_discrete_dt,file=paste0(proj_dir,"VAST/model_runs/comparison_plots/plot_empirical_discrete_ab_dt.csv"))

