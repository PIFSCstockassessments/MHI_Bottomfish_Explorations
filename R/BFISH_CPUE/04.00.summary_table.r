

# Nicholas Ducharme-Barth
# 01/06/2023
# Extract model fit characteristics from setup/ folder

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
	dir_vec = dir(paste0(proj.dir,"VAST/model_runs/"),recursive=TRUE)
	ests_vec = dir_vec[grep("setup/parameter_estimates.txt",dir_vec,fixed=TRUE)]
	complete_vec = dir_vec[grep("index_dt.csv",dir_vec,fixed=TRUE)]

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
# extract summaries		
	summary_dt.list = as.list(rep(NA,length(ests_vec)))
	for(i in seq_along(ests_vec))
	{
		tmp_date = strsplit(ests_vec[i],"/")[[1]][1]
		tmp_name = strsplit(ests_vec[i],"/")[[1]][2]
		parameter_estimates = read_vast_parameter_estimates(readLines(paste0(proj.dir,"VAST/model_runs/",tmp_date,"/",tmp_name,"/setup/parameter_estimates.txt")))

		tmp_complete = grep(paste0(tmp_date,"/",tmp_name),complete_vec,fixed=TRUE)
		if(sum(tmp_complete)>0)
		{
			parameter_estimates_final = read_vast_parameter_estimates(readLines(paste0(proj.dir,"VAST/model_runs/",tmp_date,"/",tmp_name,"/parameter_estimates.txt")))
			var_par = grep("\\bL_",parameter_estimates_final$diagnostics$Param)
			problem_par = sum(which(abs(parameter_estimates_final$diagnostics[var_par,"MLE"])<1e-3|abs(parameter_estimates_final$diagnostics[var_par,"MLE"])>1.5e1))
			rm(list=c("parameter_estimates_final","var_par"))
		} else {
			problem_par = NA
		}
		summary_dt.list[[i]] = data.table(date=tmp_date,
										  name=tmp_name) %>%
								.[,data_year:=strsplit(name,"_")[[1]][1]] %>%		  
								.[,error_structure:=strsplit(name,"_")[[1]][2]] %>%
								.[,species := strsplit(name,"_")[[1]][3]] %>%
								.[,data_treatment := strsplit(name,"_")[[1]][4]] %>%
								.[,q_config := strsplit(name,"_")[[1]][5]] %>%
								.[,ab_config := strsplit(name,"_")[[1]][6]] %>%
								.[,lehi_filter := strsplit(name,"_")[[1]][7]] %>%
								.[,fine_scale := strsplit(name,"_")[[1]][9]] %>%
								.[,all_lengths:="FALSE"] %>%
								.[,gears_used:="both"] %>%
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
								.[,complete:=ifelse(sum(tmp_complete)>0,TRUE,FALSE)] %>%
								.[,bad_param:=problem_par] %>%
								cbind(.,t(t(table(parameter_estimates$diagnostics$Param)))) %>%
								.[,V2:=NULL] %>%
								setnames(.,"V1","parameter")
		if(length(strsplit(tmp_name,"_")[[1]])>=14)
		{
			summary_dt.list[[i]]$all_lengths = strsplit(tmp_name,"_")[[1]][13]
			summary_dt.list[[i]]$gears_used = strsplit(tmp_name,"_")[[1]][14]
		} 
		
		rm(list=c("parameter_estimates","tmp_date","tmp_name","tmp_complete","problem_par"))
	}

	summary_dt = rbindlist(summary_dt.list) %>%
				 dcast(.,date+name+data_year+error_structure+species+data_treatment+q_config+ab_config+lehi_filter+fine_scale+all_lengths+gears_used+runtime+nll+n_par+n_fixed+n_random+aic+aic_all+mgc+complete+bad_param~parameter)

	fwrite(summary_dt,file=paste0(proj.dir,"VAST/model_runs/comparison_plots/summary_dt.",as.character(format(Sys.time(),format="%Y-%m-%d")),".csv"))

	