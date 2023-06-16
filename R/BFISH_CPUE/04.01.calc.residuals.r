

# Nicholas Ducharme-Barth
# 2023/06/04
# gather and calc residuals

# Copyright (c) 2023 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(VAST)
	library(FishStatsUtils)

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
    psu_table = fread(paste0(proj_dir,"Data/BFISH PSU lookup table.csv")) %>%
        .[STRATA=="SB_H_S",STRATA:="SB_A_S"] %>%
				.[,.(PSU,STRATA,Island,Depth_MEDIAN_m,med_acr,BS_pct_over_136j,pctHS)] %>%
				setnames(.,c("PSU","STRATA","Island","Depth_MEDIAN_m","med_acr","BS_pct_over_136j","pctHS"),c("psu","strata","island","depth","complexity","hardness","slope")) %>%
				.[complexity>300,complexity:=300] %>%
        .[is.na(complexity),complexity:=12.20] %>%
        .[is.na(hardness),hardness:=0.094888] %>%
				.[,islandG:=ifelse(island=="Niihau","Kauai",island)] %>%
				.[,islandG:=factor(islandG,levels=c("Big Island","Maui Nui","Oahu","Kauai"))]

#_____________________________________________________________________________________________________________________________
# extract other indices
  # order by last modified
  ests_vec = ests_vec[order(sapply(paste0(proj_dir,"VAST/model_runs/",ests_vec),file.mtime))]

  n_sim = 250
  rs = 123

  i = 57
  for(i in seq_along(residual_dt.list))
  {
	tmp_date = strsplit(ests_vec[i],"/")[[1]][1]
	tmp_name = strsplit(ests_vec[i],"/")[[1]][2]
    tmp_model = strsplit(strsplit(ests_vec[i],"/")[[1]][2],"_")[[1]]
    tmp_model_name = paste0(tmp_model[1],".",tmp_model[4],".",
                            tmp_model[5],".",
                            paste0(sapply(strsplit(tmp_model[6],"[.]"),function(x)paste0(substr(x,1,1),substr(x,nchar(x),nchar(x)))),collapse=""))
	tmp_model_name_short = paste0(ifelse(i<10,"0",""),i," ",tmp_model_name)
	tmp_data_flag = strsplit(tmp_name,"_")[[1]][1]
	tmp_data_treatment = strsplit(tmp_name,"_")[[1]][4]
	tmp_species = strsplit(tmp_name,"_")[[1]][3]
	tmp_lehi_filter = strsplit(tmp_name,"_")[[1]][7]

	# bring in bfish data
		load(file=paste0(proj_dir,"Data/",tmp_data_flag,"_",tmp_data_treatment,".bfish_combined_long_dt.RData"))
		# subset to species
		if(tmp_species == "mv")
		{ 
			tmp_target_species = c("prfi","etca","etco","prsi","przo","hyqu","apru")
		} else{
			tmp_target_species = tmp_species
		}
		tmp_bfish_df = bfish_combined_long_dt %>% .[species_cd %in% tmp_target_species] %>% as.data.frame(.)
		
		# remove sample with large lehi observation
		if(tmp_lehi_filter)
		{
			tmp_bfish_df =  subset(tmp_bfish_df,design_sampling_unit!="2021_Fall_32293")
		}
		rm(list=c("tmp_data_flag","tmp_data_treatment","tmp_species","tmp_lehi_filter","tmp_target_species"))
		
	
	# bring in model
	# ~2 minutes
	load(paste0(proj_dir,"VAST/model_runs/",tmp_date,"/",tmp_name,"/fit.RData"))
	fit = reload_model(x = fit,check_gradient=FALSE)

	sim_fit = matrix(NA,nrow=length(fit$data_list$b_i),ncol=n_sim)
	
	# calc residuals
	# ~3 hours
	for(j in 1:n_sim)
	{
		# ~20 seconds
		sim_dat = simulate_data( fit,type = 1,random_seed = list(rs+j,NULL)[[1+is.null(rs)]])
		sim_fit[,j] = sim_dat$b_i
		rm(list="sim_dat")
	}

	# save
	save_dir = paste0(proj_dir,"VAST/model_residual/",tmp_date,"/",tmp_name,"/")
	dir.create(save_dir,recursive=TRUE)
	save(sim_fit,file=paste0(save_dir,"sim_fit.RData"))


	residuals_dharma = summary(fit,what="residuals",type=1,n_samples=n_sim,random_seed=rs)
	dev.off()
	residuals_dharma$simulatedResponse = sim_fit
	residuals_dharma$observedResponse = as.vector(fit$data_list$b_i)

	residuals_dharma2 = DHARMa::createDHARMa(simulatedResponse=sim_fit,
				observedResponse=as.vector(fit$data_list$b_i),
				fittedPredictedResponse=as.vector(fit$Report$D_i),
				integer=FALSE)
	
	residuals_dharma3 = DHARMa::createDHARMa(simulatedResponse=sim_fit,
				observedResponse=as.vector(fit$data_list$b_i),
				fittedPredictedResponse=as.vector(fit$Report$D_i),
				integer=TRUE)

	# package
		resid_dt = as.data.table(tmp_bfish_df) %>%
					.[,residual:=residuals_dharma$scaledResiduals] %>%
					.[,residual2:=residuals_dharma2$scaledResiduals] %>%
					.[,residual3:=residuals_dharma3$scaledResiduals] %>%
					.[,model_name_short:=tmp_model_name_short]

	# save
	save(residuals_dharma,file=paste0(save_dir,"residuals_dharma.RData"))
	save(residuals_dharma2,file=paste0(save_dir,"residuals_dharma2.RData"))
	save(residuals_dharma3,file=paste0(save_dir,"residuals_dharma3.RData"))
	save(resid_dt,file=paste0(save_dir,"resid_dt.RData"))

    rm(list=c("tmp_model","tmp_model_name","tmp_name","tmp_date","tmp_model_name_short","fit","sim_fit","tmp_bfish_df","residuals_dharma","residuals_dharma2","residuals_dharma3","resid_dt","save_dir"))
  }    