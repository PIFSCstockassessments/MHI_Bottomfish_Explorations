

# Nicholas Ducharme-Barth
# 2024/10/23
# calc residuals for tech memo v3

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
	proj_dir = "C:/Users/Nicholas.Ducharme-Ba/Temp_Storage/HOME/SAP/2024_Deep7/"

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
  residual_flag = 1

  i = 69 # "2023-06-14/2022_lgdl_mv_05_gearSP12_depth3.hardness3.slope3.islandG_TRUE_7.5_TRUE_TRUE_pit_noxval/index_dt.csv"

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
	save_dir = paste0(proj_dir,"VAST/model_residual_v3/",tmp_date,"/",tmp_name,"/")
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
	
	# residuals_dharma3 = DHARMa::createDHARMa(simulatedResponse=sim_fit,
	# 			observedResponse=as.vector(fit$data_list$b_i),
	# 			fittedPredictedResponse=as.vector(fit$Report$D_i),
	# 			integer=TRUE)

	# package
		resid_dt = as.data.table(tmp_bfish_df) %>%
					.[,residual:=residuals_dharma$scaledResiduals] %>%
					.[,residual2:=residuals_dharma2$scaledResiduals] %>%
					# .[,residual3:=residuals_dharma3$scaledResiduals] %>%
					.[,model_name_short:=tmp_model_name_short]
    

    # plot & residual tests
        if(residual_flag==1)
        {
            residuals_object = get("residuals_dharma")
        } else{
            residuals_object = get("residuals_dharma2")
        }
			# basic QQ & residual v. predicted plot
				 	png(filename = paste0(save_dir,"dharma_agg_qq.png"),width = 16, height = 9, units = "in",res=300,bg=NA)
						plot(residuals_object)
					dev.off()
				# quantile test: residual v. predicted by model covariate
					png(filename = paste0(save_dir,"dharma_agg_quant_resid.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_quant = DHARMa::testQuantiles(residuals_object)
					dev.off()
					png(filename = paste0(save_dir,"dharma_agg_quant_resid_year.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_quant_year = DHARMa::testQuantiles(residuals_object,predictor=as.factor(resid_dt$year))
						DHARMa::plotResiduals(residuals_object, form = as.factor(resid_dt$year))
					dev.off()
					png(filename = paste0(save_dir,"dharma_agg_quant_resid_lon.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_quant_lon = DHARMa::testQuantiles(residuals_object,predictor=resid_dt$lon)
					dev.off()
					png(filename = paste0(save_dir,"dharma_agg_quant_resid_lat.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_quant_lat = DHARMa::testQuantiles(residuals_object,predictor=resid_dt$lat)
					dev.off()
					png(filename = paste0(save_dir,"dharma_agg_quant_resid_gear.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_quant_gear = DHARMa::testQuantiles(residuals_object,predictor=as.factor(resid_dt$gear_type))
					dev.off()
				# make residual plots by factor variable
					residual_factor_vec = c("island","strata","strata_2020",
											"substrate","slope","depth_strata","depth_strata_2020",
											"complexity","hardness","season","month","platform")
					for(v in 1:length(residual_factor_vec))
					{
						png(filename = paste0(save_dir,"dharma_agg_quant_resid_",residual_factor_vec[v],".png"),width = 9, height = 9, units = "in",res=300,bg=NA)
							tmp_quant_factor = DHARMa::testQuantiles(residuals_object,predictor=as.factor(as.character(as.data.frame(resid_dt)[,residual_factor_vec[v]])))
						dev.off()
					}
					residual_continuous_vec = c("depth","time","jd","lunar_phase")
					residual_continuous_vec_pvalue = rep(NA,length(residual_continuous_vec))
					for(v in 1:length(residual_continuous_vec))
					{
						png(filename = paste0(save_dir,"dharma_agg_quant_resid_",residual_continuous_vec[v],".png"),width = 9, height = 9, units = "in",res=300,bg=NA)
							tmp_quant_continuous = DHARMa::testQuantiles(residuals_object,predictor=as.data.frame(resid_dt)[,residual_continuous_vec[v]])
						dev.off()
						residual_continuous_vec_pvalue[v] = tmp_quant_continuous$p.value
						rm(list=c("tmp_quant_continuous"))
					}
				# test for uniformity in residuals, overdispersion, outliers
				 	png(filename = paste0(save_dir,"dharma_agg_residual_tests.png"),width = 16, height = 9, units = "in",res=300,bg=NA)
						tmp_test = DHARMa::testResiduals(residuals_object)
					dev.off()
				# test for zero inflation
				 	png(filename = paste0(save_dir,"dharma_agg_test_zi.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_zinf = DHARMa::testZeroInflation(residuals_object)
					dev.off()
				# test for over-dispersion
				 	png(filename = paste0(save_dir,"dharma_agg_test_od.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						DHARMa::testDispersion(residuals_object,alternative="two.sided")
					dev.off()
				# test for spatial autocorrelation
				# DHARMa::testSpatialAutocorrelation needs unique locations
					resid_dt$spatial_group = as.numeric(as.factor(paste0(resid_dt$lon,"_",resid_dt$lat)))
					tmp_spatial_group_dt = data.table(spatial_group=resid_dt$spatial_group,x=resid_dt$lon,y=resid_dt$lat) %>%
								   .[,.(x=mean(x),y=mean(y)),by=spatial_group]
					residuals_spatial_group = DHARMa::recalculateResiduals(residuals_object, group = resid_dt$spatial_group)	
				 	png(filename = paste0(save_dir,"dharma_agg_test_spcorr.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_dharma_sp = DHARMa::testSpatialAutocorrelation(residuals_spatial_group, x = tmp_spatial_group_dt$x, y = tmp_spatial_group_dt$y)
						text(min(tmp_spatial_group_dt$x),min(tmp_spatial_group_dt$y),paste0("Moran's I p-value: ",round(tmp_dharma_sp$p.value,digits=2)),adj=c(0,0))
					dev.off()
				# test for temporal autocorrelation
					resid_dt$temporal_group = factor(resid_dt$year,levels=as.character(sort(unique(resid_dt$year))))
					residuals_temporal_group = DHARMa::recalculateResiduals(residuals_object, group = resid_dt$temporal_group)	
				 	png(filename = paste0(save_dir,"dharma_agg_test_tcorr.png"),width = 9, height = 9, units = "in",res=300,bg=NA)		
						tmp_ar = DHARMa::testTemporalAutocorrelation(residuals_temporal_group, time=levels(resid_dt$temporal_group))
					dev.off()

				resid_test_agg_dt = data.table(model_name_short=tmp_model_name_short,
                            type="agg",
							KS_stat=tmp_test$uniformity$statistic,
							KS_pvalue=tmp_test$uniformity$p.value,
							dispersion_stat=tmp_test$dispersion$statistic,
							dispersion_pvalue=tmp_test$dispersion$p.value,
							outlier_stat=tmp_test$outliers$statistic,
							outlier_pvalue=tmp_test$outliers$p.value,
							zinf_stat=tmp_zinf$statistic,
							zinf_pvalue=tmp_zinf$p.value,							
							moran_pvalue=tmp_dharma_sp$p.value,
							DW_stat=tmp_ar$statistic,
							DW_pvalue=tmp_ar$p.value,
							quant_pvalue=tmp_quant$p.value,
							quant_lon_pvalue=tmp_quant_lon$p.value,
							quant_lat_pvalue=tmp_quant_lat$p.value,
							quant_depth_pvalue=residual_continuous_vec_pvalue[1],
							quant_time_pvalue=residual_continuous_vec_pvalue[2],
							quant_jd_pvalue=residual_continuous_vec_pvalue[3],
							quant_lunar_pvalue=residual_continuous_vec_pvalue[4])
				# clean-up
					rm(list=c("residual_continuous_vec_pvalue","tmp_spatial_group_dt","tmp_zinf","tmp_test","tmp_ar","tmp_dharma_sp","tmp_quant","tmp_quant_lon","tmp_quant_lat"))
		
# calc residuals by species
	    u_species = unique(resid_dt$species_cd)

		resid_test_dt.list = as.list(rep(NA,length(u_species)))

		for(j in 1:length(u_species))
		{
			resid_test_dt.list[[j]] = data.table(model_name_short=tmp_model_name_short,
                            type=u_species[j],
							KS_stat=NA,
							KS_pvalue=NA,
							dispersion_stat=NA,
							dispersion_pvalue=NA,
							outlier_stat=NA,
							outlier_pvalue=NA,
							zinf_stat=NA,
							zinf_pvalue=NA,							
							moran_pvalue=NA,
							DW_stat=NA,
							DW_pvalue=NA)

			tmp_idx = which(resid_dt$species_cd == u_species[j])
			tmp_bfish = as.data.frame(resid_dt)[tmp_idx,]
			tmp_residuals = residuals_dharma
			tmp_residuals$simulatedResponse = tmp_residuals$simulatedResponse[tmp_idx,]
			tmp_residuals$observedResponse = tmp_residuals$observedResponse[tmp_idx]
			tmp_residuals$nObs = length(tmp_idx)
			tmp_residuals$scaledResiduals = tmp_residuals$scaledResiduals[tmp_idx]
			tmp_residuals$fittedPredictedResponse = tmp_residuals$fittedPredictedResponse[tmp_idx]

			# basic QQ & residual v. predicted plot
				 	png(filename = paste0(save_dir,"dharma_",u_species[j],"_qq.png"),width = 16, height = 9, units = "in",res=300,bg=NA)
						plot(tmp_residuals)
					dev.off()
				# quantile test: residual v. predicted by model covariate
					png(filename = paste0(save_dir,"dharma_",u_species[j],"_quant_resid.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_quant = DHARMa::testQuantiles(tmp_residuals)
					dev.off()
					png(filename = paste0(save_dir,"dharma_",u_species[j],"_quant_resid_year.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_quant_year = DHARMa::testQuantiles(tmp_residuals,predictor=as.factor(tmp_bfish$year))
						DHARMa::plotResiduals(tmp_residuals, form = as.factor(tmp_bfish$year))
					dev.off()
					png(filename = paste0(save_dir,"dharma_",u_species[j],"_quant_resid_lon.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_quant_lon = DHARMa::testQuantiles(tmp_residuals,predictor=tmp_bfish$lon)
					dev.off()
					png(filename = paste0(save_dir,"dharma_",u_species[j],"_quant_resid_lat.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_quant_lat = DHARMa::testQuantiles(tmp_residuals,predictor=tmp_bfish$lat)
					dev.off()
                    png(filename = paste0(save_dir,"dharma_",u_species[j],"_quant_resid_gear.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_quant_gear = DHARMa::testQuantiles(tmp_residuals,predictor=as.factor(tmp_bfish$gear_type))
					dev.off()
					for(v in 1:length(residual_factor_vec))
					{
                    	png(filename = paste0(save_dir,"dharma_",u_species[j],"_quant_resid_",residual_factor_vec[v],".png"),width = 9, height = 9, units = "in",res=300,bg=NA)
							tmp_quant_factor = DHARMa::testQuantiles(tmp_residuals,predictor=as.factor(as.character(tmp_bfish[,residual_factor_vec[v]])))
						dev.off()
					}
					tmp_continuous_vec_pvalue = rep(NA,length(residual_continuous_vec))
					for(v in 1:length(residual_continuous_vec))
					{
                    	png(filename = paste0(save_dir,"dharma_",u_species[j],"_quant_resid_",residual_continuous_vec[v],".png"),width = 9, height = 9, units = "in",res=300,bg=NA)
							tmp_quant_continuous = DHARMa::testQuantiles(tmp_residuals,predictor=tmp_bfish[,residual_continuous_vec[v]])
						dev.off()
						tmp_continuous_vec_pvalue[v] = tmp_quant_continuous$p.value
						rm(list=c("tmp_quant_continuous"))
					}
				# test for uniformity in residuals, overdispersion, outliers
				 	png(filename = paste0(save_dir,"dharma_",u_species[j],"_residual_tests.png"),width = 16, height = 9, units = "in",res=300,bg=NA)
						tmp_test = DHARMa::testResiduals(tmp_residuals)
					dev.off()
				# test for zero inflation
				 	png(filename = paste0(save_dir,"dharma_",u_species[j],"_test_zi.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_zinf = DHARMa::testZeroInflation(tmp_residuals)
					dev.off()
				# test for over-dispersion
				 	png(filename = paste0(save_dir,"dharma_",u_species[j],"_test_od.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						DHARMa::testDispersion(tmp_residuals,alternative="two.sided")
					dev.off()
				# test for spatial autocorrelation
				# DHARMa::testSpatialAutocorrelation needs unique locations
					tmp_bfish$spatial_group = as.numeric(as.factor(paste0(tmp_bfish$lon,"_",tmp_bfish$lat)))
					tmp_spatial_group_dt = data.table(spatial_group=tmp_bfish$spatial_group,x=tmp_bfish$lon,y=tmp_bfish$lat) %>%
								   .[,.(x=mean(x),y=mean(y)),by=spatial_group]
					tmp_residuals_spatial_group = DHARMa::recalculateResiduals(tmp_residuals, group = tmp_bfish$spatial_group)	
				 	png(filename = paste0(save_dir,"dharma_",u_species[j],"_test_spcorr.png"),width = 9, height = 9, units = "in",res=300,bg=NA)
						tmp_dharma_sp = DHARMa::testSpatialAutocorrelation(tmp_residuals_spatial_group, x = tmp_spatial_group_dt$x, y = tmp_spatial_group_dt$y)
						text(min(tmp_spatial_group_dt$x),min(tmp_spatial_group_dt$y),paste0("Moran's I p-value: ",round(tmp_dharma_sp$p.value,digits=2)),adj=c(0,0))
					dev.off()
				# test for temporal autocorrelation
					tmp_bfish$temporal_group = factor(tmp_bfish$year,levels=as.character(sort(unique(tmp_bfish$year))))
					tmp_residuals_temporal_group = DHARMa::recalculateResiduals(tmp_residuals, group = tmp_bfish$temporal_group)	
				 	png(filename = paste0(save_dir,"dharma_",u_species[j],"_test_tcorr.png"),width = 9, height = 9, units = "in",res=300,bg=NA)		
						tmp_ar = DHARMa::testTemporalAutocorrelation(tmp_residuals_temporal_group, time=levels(tmp_bfish$temporal_group))
					dev.off()


				resid_test_dt.list[[j]]$KS_stat=tmp_test$uniformity$statistic
				resid_test_dt.list[[j]]$KS_pvalue=tmp_test$uniformity$p.value
				resid_test_dt.list[[j]]$dispersion_stat=tmp_test$dispersion$statistic
				resid_test_dt.list[[j]]$dispersion_pvalue=tmp_test$dispersion$p.value
				resid_test_dt.list[[j]]$outlier_stat=tmp_test$outliers$statistic
				resid_test_dt.list[[j]]$outlier_pvalue=tmp_test$outliers$p.value
				resid_test_dt.list[[j]]$zinf_stat=tmp_zinf$statistic
				resid_test_dt.list[[j]]$zinf_pvalue=tmp_zinf$p.value				
				resid_test_dt.list[[j]]$moran_pvalue=tmp_dharma_sp$p.value
				resid_test_dt.list[[j]]$DW_stat=tmp_ar$statistic
				resid_test_dt.list[[j]]$DW_pvalue=tmp_ar$p.value
				resid_test_dt.list[[j]]$quant_pvalue=tmp_quant$p.value
				resid_test_dt.list[[j]]$quant_lon_pvalue=tmp_quant_lon$p.value
				resid_test_dt.list[[j]]$quant_lat_pvalue=tmp_quant_lat$p.value
				resid_test_dt.list[[j]]$quant_depth_pvalue=tmp_continuous_vec_pvalue[1]
				resid_test_dt.list[[j]]$quant_time_pvalue=tmp_continuous_vec_pvalue[2]
				resid_test_dt.list[[j]]$quant_jd_pvalue=tmp_continuous_vec_pvalue[3]
				resid_test_dt.list[[j]]$quant_lunar_pvalue=tmp_continuous_vec_pvalue[4]

				# clean-up
					rm(list=c("tmp_zinf","tmp_test","tmp_ar","tmp_idx","tmp_bfish","tmp_residuals","tmp_spatial_group_dt","tmp_dharma_sp","tmp_residuals_spatial_group","tmp_residuals_temporal_group","tmp_quant","tmp_quant_lon","tmp_quant_lat","tmp_continuous_vec_pvalue"))
		}
		
		resid_test_dt = rbind(resid_test_agg_dt,rbindlist(resid_test_dt.list)) %>%
	                .[,type:=deep7_name_vec[match(type,deep7_code_vec)]] %>%
                    .[is.na(type),type:="Total"] %>%
                    .[,type:=factor(type,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]


	# save
	save(residuals_dharma,file=paste0(save_dir,"residuals_dharma.RData"))
	save(residuals_dharma2,file=paste0(save_dir,"residuals_dharma2.RData"))
	# save(residuals_dharma3,file=paste0(save_dir,"residuals_dharma3.RData"))
	save(resid_dt,file=paste0(save_dir,"resid_dt.RData"))
	save(resid_test_dt,file=paste0(save_dir,"resid_test_dt.RData"))

    rm(list=c("resid_test_dt.list","resid_test_agg_dt","resid_test_dt","u_species","tmp_model","tmp_model_name","tmp_name","tmp_date","tmp_model_name_short","fit","sim_fit","tmp_bfish_df","residuals_dharma","residuals_dharma2","resid_dt","save_dir","residuals_object"))
    gc()
#   }    
