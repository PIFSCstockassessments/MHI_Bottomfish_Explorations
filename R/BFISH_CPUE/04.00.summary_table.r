

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
	ests_vec = dir_vec[grep("setup/parameter_estimates.RData",dir_vec,fixed=TRUE)]

	summary_dt.list = as.list(rep(NA,length(ests_vec)))
	for(i in seq_along(ests_vec))
	{
		load(paste0(proj.dir,"VAST/model_runs/",ests_vec[i]))

		summary_dt.list[[i]] = data.table(date=strsplit(ests_vec[i],"/")[[1]][1],
										  name=strsplit(ests_vec[i],"/")[[1]][2]) %>%
								.[,data_year:=strsplit(name,"_")[[1]][1]] %>%		  
								.[,error_structure:=strsplit(name,"_")[[1]][2]] %>%
								.[,species := strsplit(name,"_")[[1]][3]] %>%
								.[,data_treatment := strsplit(name,"_")[[1]][4]] %>%
								.[,q_config := strsplit(name,"_")[[1]][5]] %>%
								.[,ab_config := strsplit(name,"_")[[1]][6]] %>%
								.[,lehi_filter := strsplit(name,"_")[[1]][7]] %>%
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
								.[,pdh:=parameter_estimates$SD$pdHess] %>%
								cbind(.,t(t(table(parameter_estimates$diagnostics$Param)))) %>%
								.[,V2:=NULL] %>%
								setnames(.,"V1","parameter") 
		
		rm(list=c("parameter_estimates"))
	}

	summary_dt = rbindlist(summary_dt.list) %>%
				 dcast(.,date+name+data_year+error_structure+species+data_treatment+q_config+ab_config+lehi_filter+fine_scale+runtime+nll+n_par+n_fixed+n_random+aic+aic_all+mgc+pdh~parameter)

	fwrite(summary_dt,file=paste0(proj.dir,"VAST/model_runs/comparison_plots/summary_dt.",as.character(format(Sys.time(),format="%Y-%m-%d")),".csv"))

	