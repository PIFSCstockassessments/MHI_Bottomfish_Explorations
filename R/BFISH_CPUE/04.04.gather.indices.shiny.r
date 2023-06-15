

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
				.[,.(name,date,data_year,error_structure,species,data_treatment,q_config,ab_config,lehi_filter,knot_dist,fine_scale,model_name_short,category,time,estimate,cv,l95,u95,estimate_sc,l95_sc,u95_sc)]


#_____________________________________________________________________________________________________________________________
# extract other indices
  # order by last modified
  ests_vec = ests_vec[order(sapply(paste0(proj_dir,"VAST/model_runs/",ests_vec),file.mtime))]

  index_summary_dt.list = index_dt.list = as.list(rep(NA,length(ests_vec)))
  
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
						.[,date:=as.POSIXct(date)] %>%
						 setnames(.,c("Model","Category","Time","Estimate","CV"),c("model_name_short","category","time","estimate","cv")) %>%
                     .[,model_name_short:=paste0(ifelse(i<10,"0",""),i," ",tmp_model_name)] %>%
                     .[,category:=deep7_name_vec[match(category,deep7_code_vec)]] %>%
                     .[is.na(category),category:="Total"] %>%
                     .[,category:=factor(category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))] %>%
					 .[,estimate_sc:=estimate/mean(estimate),by=.(category)] %>%
			         .[,l95_sc:=exp(log(estimate_sc)-2*sqrt(log(cv^2+1)))] %>%
			         .[,u95_sc:=exp(log(estimate_sc)+2*sqrt(log(cv^2+1)))] %>%
					 .[,.(name,date,data_year,error_structure,species,data_treatment,q_config,ab_config,lehi_filter,knot_dist,fine_scale,model_name_short,category,time,estimate,cv,l95,u95,estimate_sc,l95_sc,u95_sc)]
	
		parameter_estimates = read_vast_parameter_estimates(readLines(paste0(proj_dir,"VAST/model_runs/",tmp_date,"/",tmp_name,"/setup/parameter_estimates.txt")))
		parameter_estimates_final = read_vast_parameter_estimates(readLines(paste0(proj_dir,"VAST/model_runs/",tmp_date,"/",tmp_name,"/parameter_estimates.txt")))
		var_par = grep("\\bL_",parameter_estimates_final$diagnostics$Param)
		problem_par = sum(which(abs(parameter_estimates_final$diagnostics[var_par,"MLE"])<1e-3|abs(parameter_estimates_final$diagnostics[var_par,"MLE"])>1.5e1))

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

    rm(list=c("tmp_model","tmp_model_name","tmp_name","tmp_date","parameter_estimates_final","var_par","problem_par","parameter_estimates"))
  }

  index_dt = rbind(design_dt,rbindlist(index_dt.list)) %>%
			.[,model_number:=as.numeric(sapply(model_name_short,function(x)strsplit(x,"\\s+")[[1]][1]))] %>%
			.[order(model_number,time,category)]

  index_summary_dt = rbindlist(index_summary_dt.list) %>%
				 dcast(.,date+name+model_name_short+data_year+error_structure+species+data_treatment+q_config+ab_config+lehi_filter+knot_dist+fine_scale+runtime+nll+n_par+n_fixed+n_random+aic+aic_all+mgc+bad_param~parameter) %>%
					merge(.,data.table(model_name_short="00 design",q_config="design",ab_config="design"),by=c("model_name_short","q_config","ab_config"),all=TRUE) %>%
  					.[,model_number:=as.numeric(sapply(model_name_short,function(x)strsplit(x,"\\s+")[[1]][1]))] %>%
					.[order(model_number)] %>%
					.[,.(model_number,model_name_short,q_config,ab_config,data_year,error_structure,species,data_treatment,lehi_filter,knot_dist,fine_scale,runtime,nll,n_par,n_fixed,n_random,aic,aic_all,mgc,bad_param,L_epsilon1_z,L_epsilon2_z,L_eta1_z,L_eta2_z,L_omega1_z,L_omega2_z,beta1_ft,beta2_ft,gamma1_cp,lambda1_k,lambda2_k,logSigmaM,logkappa1,logkappa2)]

#_____________________________________________________________________________________________________________________________
# save
	fwrite(index_dt,file=paste0(proj_dir,"VAST/model_runs/comparison_plots/index_dt.csv"))
	fwrite(index_summary_dt,file=paste0(proj_dir,"VAST/model_runs/comparison_plots/index_summary_dt.csv"))


