

# Nicholas Ducharme-Barth
# 09/08/2022
# Set up simple spatiotemporal model example using VAST
# 1) Bring data (including spatial data)
# 2) Set-up barrier mesh (including conversion to equal distant projection)
# 3) Fit multivariate model using Poisson-link delta-Gamma model, add functionality to sum biomass across categories
# 4) Examine diagnostics
# 5) Calculate index

# This combines topics covered in VAST vignettes (https://github.com/James-Thorson-NOAA/VAST/wiki)

# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(VAST)
	library(ggplot2)
	library(ggthemes)
	# rgdal functions sourced directly
	# sp functions sourced directly
	# DHARMa functions sourced directly
	# sf functions sourced directly
	# ggthemes functions sourced directly
#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"
	working_dir = paste0(proj.dir,"VAST/model_runs/pl-dgamma_mv7_vanilla/")
	dir.create(working_dir,recursive=TRUE)

#_____________________________________________________________________________________________________________________________
# define data_flag
	# data_flag = "" # only loads data up through 2020
	data_flag = "2021_" # includes data through 2021
#_____________________________________________________________________________________________________________________________
# 1) bring in data
	load(file=paste0(proj.dir,"Data/",data_flag,"bfish_combined_long_dt.RData"))
	# collapse weight_kg across bait_type
		sample_weight_dt = bfish_combined_long_dt[,.(weight_kg=sum(weight_kg)),by=.(sample_id,species_cd)]
		unique_dt = copy(bfish_combined_long_dt) %>% .[,bait_type:=NULL] %>% .[,weight_kg:=NULL] %>% unique(.)
		bfish_combined_long_dt = merge(unique_dt,sample_weight_dt,by=c("sample_id","species_cd"))
	
	# subset to species
	bfish_df = bfish_combined_long_dt %>% .[species_cd %in% c("prfi","etco","etca","prsi","przo","hyqu","apru")] %>% as.data.frame(.)

	# remove sample with large lehi observation
	bfish_df =  subset(bfish_df,sample_id!="20210831_183445")

	# needed to define spatial domain and for predicting on to create index
	psu_table = fread(paste0(proj.dir,"Data/BFISH PSU lookup table.csv")) %>%
				.[,.(PSU,Island,lon_deg,lat_deg,STRATA,STRATA_2020,Depth_MEDIAN_m,med_slp,med_acr,BS_pct_over_136j,pctHB,pctHS)] %>%
				.[,substrate:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][1])] %>%
				.[,slope:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][2])]
	# needed for plotting and defining the barrier feature
	hi_coast = rgdal::readOGR(dsn = paste0(proj.dir,"Data/GIS/Coastline"), layer = "Coastline")
	hi_coast_sf = sf::st_as_sf(hi_coast)

#_____________________________________________________________________________________________________________________________
# 2) check for outliers & remove

#_____________________________________________________________________________________________________________________________
# 3) set-up spatial domain
	# convert to equal distant projection
		# get original lat-lon crs
		crs_ll = sp::CRS(sp::proj4string(hi_coast))
		# use two-point equi-distant projection
		crs_eqd = sp::CRS("+proj=tpeqd +lat_1=20.45 +lon_1=-158.9 +lat_2=20.45 +lon_2=-156.2 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
	  	hi_coast_eqd = sp::spTransform(hi_coast, crs_eqd)
	  	hi_coast_eqd_sf = sf::st_as_sf(hi_coast_eqd) 

	# convert to equi-distant coordinates for samples & psu
		bfish_coords = cbind(bfish_df$lon, bfish_df$lat)
		bfish_sp = sp::SpatialPoints(bfish_coords)
		sp::proj4string(bfish_sp) = crs_ll
		bfish_sp_eqd = sp::spTransform(bfish_sp, crs_eqd)
		bfish_df$lon_eqd = bfish_sp_eqd@coords[,1]
		bfish_df$lat_eqd = bfish_sp_eqd@coords[,2]

		psu_coords = cbind(psu_table$lon_deg, psu_table$lat_deg)
		psu_sp = sp::SpatialPoints(psu_coords)
		sp::proj4string(psu_sp) = crs_ll
		psu_sp_eqd = sp::spTransform(psu_sp, crs_eqd)
		psu_table$lon_eqd = psu_sp_eqd@coords[,1]
		psu_table$lat_eqd = psu_sp_eqd@coords[,2]

	# define extrapolation grid based on PSUs
		input_grid = cbind(psu_coords[,2],psu_coords[,1],0.5^2)
		colnames(input_grid) = c("Lat","Lon","Area_km2")
		Extrapolation_List = make_extrapolation_info(Region="user",projargs=slot(crs_eqd,"projargs"),input_grid=input_grid)
	# add strata definitions for each Island group, STRATA and STRATA_2020
		a_el_orig = Extrapolation_List$a_el[,1]
		a_el_tmp = matrix(0,nrow=length(a_el_orig),ncol=1+length(unique(psu_table$Island))+length(unique(psu_table$STRATA))+length(unique(psu_table$STRATA_2020)))
		colnames(a_el_tmp) = c("All_areas","Niihau","Kauai","Oahu","Maui Nui","Big Island",sort(unique(psu_table$STRATA)),sort(unique(psu_table$STRATA_2020)))
		a_el_tmp[,1] = a_el_orig
		for(i in 1:5)
		{
			a_el_tmp[which(psu_table$Island == colnames(a_el_tmp)[1+i]),1+i] = a_el_orig[which(psu_table$Island == colnames(a_el_tmp)[1+i])]
		}
		for(i in 1:length(unique(psu_table$STRATA)))
		{
			a_el_tmp[which(psu_table$STRATA == colnames(a_el_tmp)[6+i]),6+i] = a_el_orig[which(psu_table$STRATA == colnames(a_el_tmp)[6+i])]
		}
		for(i in 1:length(unique(psu_table$STRATA_2020)))
		{
			a_el_tmp[which(psu_table$STRATA_2020 == colnames(a_el_tmp)[16+i]),16+i] = a_el_orig[which(psu_table$STRATA_2020 == colnames(a_el_tmp)[16+i])]
		}

		Extrapolation_List$a_el = a_el_tmp
		units(Extrapolation_List$a_el) = "km^2"

	# define input mesh
		# setting the boundary helps restrict the placement of knots to where the data is
		# the closer to zero convex is, the more intricate the boundary can become and the larger resolution is needed
		# also a more intricate boundary will result in longer mesh generation time and a greater number of knots
		# this takes a few seconds
		A = proc.time()
		mesh_boundary = INLA::inla.nonconvex.hull(cbind(psu_table$lon_eqd, psu_table$lat_eqd), convex = -0.005,resolution=205)
		B = proc.time()
		print(paste0("mesh boundary took ",round((B-A)[3],digits=2)," seconds to generate"))
		# in inla.mesh.2d the max.n.strict argument restricts the number of knots more than setting max.n. However, the cutoff distance between knots is what really controls the number of knots.
		# for this example leave cutoff == 7.5 (in km) to keep the number of knots manageable, but for full analysis can decrease to << 4.
		mesh_inla = INLA::inla.mesh.2d(loc = psu_sp_eqd,boundary = mesh_boundary,max.n.strict=c(500,16),cutoff=7.5)

		# define knot locations for spatial_list
			mesh_coords_eqd = sp::SpatialPoints(mesh_inla$loc[,1:2])
			sp::proj4string(mesh_coords_eqd) = crs_eqd
			mesh_coords_ll = sp::spTransform(mesh_coords_eqd, crs_ll)
			intensity_loc = mesh_coords_ll@coords 
		
		fine_scale = FALSE
		spatial_list = make_spatial_info(n_x=nrow(mesh_coords_ll@coords),
    					   Lon_i=bfish_df$lon,
    					   Lat_i=bfish_df$lat,
    					   Extrapolation_List=Extrapolation_List,
    					   Method = "Barrier",
    					   anisotropic_mesh = mesh_inla,
    					   grid_size_km = 0.5,
    					   grid_size_LL = 0.5/110,
    					   fine_scale = fine_scale,
    					   Save_Results = FALSE,
    					   LON_intensity=intensity_loc[,1],
    					   LAT_intensity=intensity_loc[,2],
    					   map_data=hi_coast)

		# add plot to show the mesh structure
		# pull-out needed quantities for plotting
			t.sub=1:nrow(spatial_list$MeshList$anisotropic_mesh$graph$tv)
		    idx = cbind(spatial_list$MeshList$anisotropic_mesh$graph$tv[t.sub,c(1:3,1), drop=FALSE], NA)
		    x = spatial_list$MeshList$anisotropic_mesh$loc[t(idx), 1]
	        y = spatial_list$MeshList$anisotropic_mesh$loc[t(idx), 2]
		    idx_normal = idx[-spatial_list$MeshList$anisotropic_mesh_triangles_over_land,]
	        x_normal = spatial_list$MeshList$anisotropic_mesh$loc[t(idx_normal), 1]
	        y_normal =spatial_list$MeshList$anisotropic_mesh$loc[t(idx_normal), 2]
	        idx_boundary = mesh_boundary$idx
	        # add NA at each non-consecutive index
	        non_consec_idx = c(1,which(diff(idx_boundary[,2])<0)+1)
	        new_idx_boundary = matrix(c(NA,1),nrow=1,ncol=2)
	        for(i in 1:(length(non_consec_idx)-1))
	        {
	        	new_idx_boundary = rbind(new_idx_boundary,idx_boundary[non_consec_idx[i]:non_consec_idx[i+1],],c(idx_boundary[non_consec_idx[i+1],1],NA))
	        }
	        # remove duplicate indices for plotting
	        new_idx_boundary = unique(new_idx_boundary)
	        x_boundary = mesh_boundary$loc[t(new_idx_boundary), 1]
	        y_boundary = mesh_boundary$loc[t(new_idx_boundary), 2] 
 			
 			normal_col = "black"
 			barrier_col = "gray70"
 			boundary_col = "hotpink"
 			land_col = "gray90"
	 		
	 		png(filename = paste0(working_dir,"spatial_mesh.png"),width = 16, height = 9, units = "in",res=300)
	 		par(mar=c(5,5,1,1))
			sp::plot(hi_coast_eqd,axes=TRUE,col=land_col,ylim=c(-300,300),xlab="Eastings (km)",ylab="Northings (km)",cex=1.5,cex.axis=1.5,cex.lab=1.5,las=1)
			lines(x,y,col=barrier_col)
			points(x,y,pch=16,cex=0.5,col=barrier_col)
			lines(x_normal,y_normal,col=normal_col)
			lines(x_boundary,y_boundary,col=boundary_col,lwd=3)
			points(x_normal,y_normal,pch=16,cex=0.5,col=normal_col)
			legend("bottomleft",legend=c("normal knot","barrier knot","normal edge","barrier edge","knot boundary","land"),lwd=c(NA,NA,1,1,3,1),pch=c(16,16,NA,NA,NA,NA),col=c(normal_col,barrier_col,normal_col,barrier_col,boundary_col,NA),fill=c(NA,NA,NA,NA,NA,land_col),border=c(NA,NA,NA,NA,NA,"black"),bty="n",cex=1.5)
			dev.off()
#_____________________________________________________________________________________________________________________________
# 4) fit a basic VAST model
# Error structure: delta-gamma
# Fixed effects: Year
# Random effects: spatial + spatiotemporal
	
	# make mean-variance plots
		meanvar_dt = as.data.table(bfish_df) %>%
					 .[weight_kg>0] %>%
					 .[,.(mean=log(mean(weight_kg)),var=log(var(weight_kg)),.N),by=.(year,species_cd)] %>%
					 .[!is.na(var)&is.finite(var)]
		meanvar_lm = lm(formula = meanvar_dt$var ~ meanvar_dt$mean)
		text_dt = data.table(x=1,y=0,label=paste0("slope = ",round(meanvar_lm$coefficients[2],digits=2)))


		p = copy(meanvar_dt) %>%
			.[,species_cd:=factor(species_cd,levels=c("prfi","etca","etco","prsi","przo","hyqu","apru"),labels=c("Opakapaka", "Ehu", "Onaga", "Kalekale", "Gindai", "Hapuupuu", "Lehi"))] %>%
			ggplot() +
			xlab("log(Mean)") +
			ylab("log(Variance)") +
			geom_abline(slope=1,intercept=meanvar_lm$coefficients[1],linetype="dashed",color="gray70") +
			geom_abline(slope=2,intercept=meanvar_lm$coefficients[1],linetype="dashed",color="gray70") +
			geom_abline(slope=3,intercept=meanvar_lm$coefficients[1],linetype="dashed",color="gray70") +
			geom_abline(slope=meanvar_lm$coefficients[2],intercept=meanvar_lm$coefficients[1],linetype="dashed",color="black",size=2) +			
			geom_point(aes(x=mean,y=var,fill=species_cd,size=N),shape=21) +
			theme_few(base_size=20) +
			geom_text(data=text_dt,aes(x=x,y=y,label=label),size=7) +
     		viridis::scale_fill_viridis("Species\ncode",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
			ggsave(filename=paste0("mean_var.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	# make settings
		bias.correct = FALSE
		settings = make_settings(n_x=spatial_list$n_x,
								 Region="user",
								 fine_scale=fine_scale,
								 purpose="index2",
								 use_anisotropy = FALSE,
								 FieldConfig=matrix( c("IID","IID","IID","IID","IID","IID"), ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) ),
								 Options=c("treat_nonencounter_as_zero"=TRUE ),
								 bias.correct=bias.correct,
								 max_cells=Inf,
								 ObsModel=c(2,4))
		settings$grid_size_km = 0.5
		settings$Options = c( settings$Options, "range_fraction"=0.01 )

		# settings$Options = c( settings$Options, "report_additional_variables"=TRUE )

	# set-up model
		fit_setup = fit_model( settings=settings,
 				   				Lon_i=bfish_df$lon,
    			   				Lat_i=bfish_df$lat,
          						t_i=as.integer(bfish_df$year),
          						b_i=bfish_df$weight_kg,
          						a_i=rep(pi * (0.02760333457^2),nrow(bfish_df)), # assumed area swept from the MOUSS camera converted to km2; Ault et al 2018
	    	  					c_i = as.numeric(factor(bfish_df[,'species_cd'],levels=c("prfi","etca","etco","prsi","przo","hyqu","apru")))-1,
	    	  					category_names=c("prfi","etca","etco","prsi","przo","hyqu","apru"),
	    	  					covariate_data = NULL,
	    	  					catchability_data = NULL,
          						working_dir = working_dir,
          						newtonsteps = 0,
          						# extrapolation list args
          						projargs=slot(crs_eqd,"projargs"),input_grid=input_grid,
          						extrapolation_list = Extrapolation_List,
          						# spatial list args
    	  						Method = "Barrier",anisotropic_mesh = mesh_inla,grid_size_LL = 0.5/110,Save_Results = FALSE,LON_intensity=intensity_loc[,1],LAT_intensity=intensity_loc[,2],
    	  						spatial_list = spatial_list,build_model=TRUE,test_fit=FALSE)
		# fit_setup$parameter_estimates$SD
		
		# turn off estimation of variance for spatiotemporal random effect for 3rd category (species)
		# and fix to low value since goes to zero when freely estimated
		modified_map = fit_setup$tmb_list$Map
		modified_map$L_omega1_z = factor(c(1,2,3,4,5,NA,NA),levels=c(1,2,3,4,5))
		modified_map$L_epsilon1_z = factor(c(1,2,NA,3,4,5,6),levels=c(1,2,3,4,5,6))
		modified_map$L_omega2_z = factor(c(1,2,3,NA,4,NA,NA),levels=c(1,2,3,4))
		modified_map$L_epsilon2_z = factor(c(1,2,NA,3,4,5,NA),levels=c(1,2,3,4,5))

		modified_parameters = fit_setup$tmb_list$Parameters
		modified_parameters$L_omega1_z = c(1,1,1,1,1,0.00001,0.00001)
		modified_parameters$L_epsilon1_z = c(1,1,0.00001,1,1,1,1)
		modified_parameters$L_omega2_z = c(1,1,1,0.00001,1,0.00001,0.00001)
		modified_parameters$L_epsilon2_z = c(1,1,0.00001,1,1,1,0.00001)

		fit = fit_model( settings=settings,
 				   				Lon_i=bfish_df$lon,
    			   				Lat_i=bfish_df$lat,
          						t_i=as.integer(bfish_df$year),
          						b_i=bfish_df$weight_kg,
          						a_i=rep(pi * (0.02760333457^2),nrow(bfish_df)), # assumed area swept from the MOUSS camera converted to km2; Ault et al 2018
	    	  					c_i = as.numeric(factor(bfish_df[,'species_cd'],levels=c("prfi","etca","etco","prsi","przo","hyqu","apru")))-1,
	    	  					category_names=c("prfi","etca","etco","prsi","przo","hyqu","apru"),
	    	  					covariate_data = NULL,
	    	  					catchability_data = NULL,
          						working_dir = working_dir,
          						newtonsteps = 1,
          						# extrapolation list args
          						projargs=slot(crs_eqd,"projargs"),input_grid=input_grid,
          						extrapolation_list = Extrapolation_List,
          						# spatial list args
    	  						Method = "Barrier",anisotropic_mesh = mesh_inla,grid_size_LL = 0.5/110,Save_Results = FALSE,LON_intensity=intensity_loc[,1],LAT_intensity=intensity_loc[,2],
    	  						spatial_list = spatial_list,
    	  						Map = modified_map,
    	  						Parameters = modified_parameters,
    	  						test_fit=TRUE)

		# re-fit model using Map and Parameters from fit in order to get aggregate index
		dir.create(paste0(working_dir,"agg/"),recursive=TRUE)
		fit_agg = fit_model( settings=settings,
 				   				Lon_i=bfish_df$lon,
    			   				Lat_i=bfish_df$lat,
          						t_i=as.integer(bfish_df$year),
          						b_i=bfish_df$weight_kg,
          						a_i=rep(pi * (0.02760333457^2),nrow(bfish_df)), # assumed area swept from the MOUSS camera converted to km2; Ault et al 2018
	    	  					c_i = as.numeric(factor(bfish_df[,'species_cd'],levels=c("prfi","etca","etco","prsi","przo","hyqu","apru")))-1,
	    	  					category_names=c("prfi","etca","etco","prsi","przo","hyqu","apru"),
   	  						  	Expansion_cz = cbind( c(0,3,3,3,3,3,3), c(0,0,1,2,3,4,5)),
	    	  					covariate_data = NULL,
	    	  					catchability_data = NULL,
          						working_dir = paste0(working_dir,"agg/"),
          						newtonsteps = 1,
          						# extrapolation list args
          						projargs=slot(crs_eqd,"projargs"),input_grid=input_grid,
          						extrapolation_list = Extrapolation_List,
          						# spatial list args
    	  						Method = "Barrier",anisotropic_mesh = mesh_inla,grid_size_LL = 0.5/110,Save_Results = FALSE,LON_intensity=intensity_loc[,1],LAT_intensity=intensity_loc[,2],
    	  						spatial_list = spatial_list,
    	  						Map = fit$tmb_list$Map,
    	  						Parameters = fit$tmb_list$Parameters)

		# Sdreport = fit$parameter_estimates$SD
		# par_biascorrect = TMB:::as.list.sdreport( Sdreport, what="Est. (bias.correct)", report=TRUE )

		index = plot_biomass_index( fit,DirName=working_dir)
		index_agg = plot_biomass_index( fit_agg,DirName=paste0(working_dir,"agg/"))

	save(fit,file=paste0(working_dir,"fit.RData"))
		save(fit_agg,file=paste0(working_dir,"agg/fit_agg.RData"))
		bfish_df$pred_weight_kg = fit$Report$D_i
		bfish_dt = as.data.table(bfish_df)		
		save(bfish_dt,file=paste0(working_dir,"bfish_dt.RData"))
		save(spatial_list,file=paste0(working_dir,"spatial_list.RData"))
		save(Extrapolation_List,file=paste0(working_dir,"Extrapolation_List.RData"))

		# calculate RMSE
		rmse_agg_dt = copy(bfish_dt) %>%
				  .[,.(rmse=sqrt(mean((weight_kg-pred_weight_kg)^2)))] %>%
				  .[,type:="agg"]
		rmse_spec_dt = copy(bfish_dt) %>%
				  .[,.(rmse=sqrt(mean((weight_kg-pred_weight_kg)^2))),by=species_cd] %>%
				  setnames(.,"species_cd","type")
		rmse_year_dt = copy(bfish_dt) %>%
				  .[,.(rmse=sqrt(mean((weight_kg-pred_weight_kg)^2))),by=year] %>%
				  setnames(.,"year","type")
		rmse_island_dt = copy(bfish_dt) %>%
				  .[,.(rmse=sqrt(mean((weight_kg-pred_weight_kg)^2))),by=island] %>%
				  setnames(.,"island","type")
		rmse_gear_dt = copy(bfish_dt) %>%
				  .[,.(rmse=sqrt(mean((weight_kg-pred_weight_kg)^2))),by=gear_type] %>%
				  setnames(.,"gear_type","type")
		rmse_strata_dt = copy(bfish_dt) %>%
				  .[,.(rmse=sqrt(mean((weight_kg-pred_weight_kg)^2))),by=strata] %>%
				  setnames(.,"strata","type")				  
		rmse_dt = rbind(rmse_agg_dt,rmse_spec_dt,rmse_year_dt,rmse_island_dt,rmse_gear_dt,rmse_strata_dt)
		fwrite(rmse_dt,file=paste0(working_dir,"rmse_dt.csv"))

#_____________________________________________________________________________________________________________________________
# 5) Model diagnostics & plot output
	# need to simulate responses to calculate DHARMa style residuals
	# can use built in functions from sdmTMB package but DHARMa package has some built in tests/plots
	n_sim = 500
	sim_fit = matrix(NA,nrow=length(fit$data_list$b_i),ncol=n_sim)
	rs = 123
	for(i in 1:n_sim)
	{
		sim_dat = simulate_data( fit,type = 1,random_seed = list(rs+i,NULL)[[1+is.null(rs)]])
		sim_fit[,i] = sim_dat$b_i
	}

	save(sim_fit,file=paste0(working_dir,"sim_fit.RData"))


	pred = fit$Report$D_i
	residuals_dharma =  DHARMa::createDHARMa(simulatedResponse = sim_fit,observedResponse = bfish_df$weight_kg,fittedPredictedResponse = pred)	

	residual_dt = data.table(type="agg",
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

	# basic QQ & residual v. predicted plot
	 	png(filename = paste0(working_dir,"dharma_agg_qq.png"),width = 16, height = 9, units = "in",res=300)
			plot(residuals_dharma)
		dev.off()
	# residual v. predicted by model covariate
		png(filename = paste0(working_dir,"dharma_agg_resid_year.png"),width = 9, height = 9, units = "in",res=300)
			DHARMa::plotResiduals(residuals_dharma, form = as.factor(bfish_df$year))
		dev.off()
		png(filename = paste0(working_dir,"dharma_agg_resid_lon.png"),width = 9, height = 9, units = "in",res=300)
			DHARMa::plotResiduals(residuals_dharma, form = bfish_df$lon)
		dev.off()
		png(filename = paste0(working_dir,"dharma_agg_resid_lat.png"),width = 9, height = 9, units = "in",res=300)
			DHARMa::plotResiduals(residuals_dharma, form = bfish_df$lat)
		dev.off()
	# test for uniformity in residuals, overdispersion, outliers
	 	png(filename = paste0(working_dir,"dharma_agg_residual_tests.png"),width = 16, height = 9, units = "in",res=300)
			resid_test = DHARMa::testResiduals(residuals_dharma)
		dev.off()
	# test for zero inflation
	 	png(filename = paste0(working_dir,"dharma_agg_test_zi.png"),width = 9, height = 9, units = "in",res=300)
			zinf_test = DHARMa::testZeroInflation(residuals_dharma)
		dev.off()
	# test for over-dispersion
	 	png(filename = paste0(working_dir,"dharma_agg_test_od.png"),width = 9, height = 9, units = "in",res=300)
			DHARMa::testDispersion(residuals_dharma,alternative="two.sided")
		dev.off()
	# test for spatial autocorrelation
	# DHARMa::testSpatialAutocorrelation needs unique locations
		bfish_df$spatial_group = as.numeric(as.factor(paste0(bfish_df$lon,"_",bfish_df$lat)))
		spatial_group_dt = data.table(spatial_group=bfish_df$spatial_group,x=bfish_df$lon,y=bfish_df$lat) %>%
					   .[,.(x=mean(x),y=mean(y)),by=spatial_group]
		residuals_spatial_group = DHARMa::recalculateResiduals(residuals_dharma, group = bfish_df$spatial_group)	
	 	png(filename = paste0(working_dir,"dharma_agg_test_spcorr.png"),width = 9, height = 9, units = "in",res=300)
			dharma_sp = DHARMa::testSpatialAutocorrelation(residuals_spatial_group, x = spatial_group_dt$x, y = spatial_group_dt$y)
			text(min(spatial_group_dt$x),min(spatial_group_dt$y),paste0("Moran's I p-value: ",round(dharma_sp$p.value,digits=2)),adj=c(0,0))
		dev.off()
	# test for temporal autocorrelation
		bfish_df$temporal_group = factor(bfish_df$year,levels=as.character(sort(unique(bfish_df$year))))
		residuals_temporal_group = DHARMa::recalculateResiduals(residuals_dharma, group = bfish_df$temporal_group)	
	 	png(filename = paste0(working_dir,"dharma_agg_test_tcorr.png"),width = 9, height = 9, units = "in",res=300)		
			resid_autocorr_test = DHARMa::testTemporalAutocorrelation(residuals_temporal_group, time=levels(bfish_df$temporal_group))
		dev.off()


		residual_dt$KS_stat=resid_test$uniformity$statistic
		residual_dt$KS_pvalue=resid_test$uniformity$p.value
		residual_dt$dispersion_stat=resid_test$dispersion$statistic
		residual_dt$dispersion_pvalue=resid_test$dispersion$p.value
		residual_dt$outlier_stat=resid_test$outliers$statistic
		residual_dt$outlier_pvalue=resid_test$outliers$p.value
		residual_dt$zinf_stat=zinf_test$statistic
		residual_dt$zinf_pvalue=zinf_test$p.value
		residual_dt$moran_pvalue=dharma_sp$p.value
		residual_dt$DW_stat=resid_autocorr_test$statistic
		residual_dt$DW_pvalue=resid_autocorr_test$p.value

	# calc residuals by species
		u_species = unique(bfish_df$species_cd)

		residual_dt.list = as.list(rep(NA,length(u_species)))
		for(i in 1:length(u_species))
		{
			residual_dt.list[[i]] = data.table(type=u_species[i],
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

			tmp_idx = which(bfish_df$species_cd == u_species[i])
			tmp_bfish = bfish_df[tmp_idx,]
			tmp_residuals = residuals_dharma
			tmp_residuals$simulatedResponse = tmp_residuals$simulatedResponse[tmp_idx,]
			tmp_residuals$observedResponse = tmp_residuals$observedResponse[tmp_idx]
			tmp_residuals$nObs = length(tmp_idx)
			tmp_residuals$scaledResiduals = tmp_residuals$scaledResiduals[tmp_idx]
			tmp_residuals$fittedPredictedResponse = tmp_residuals$fittedPredictedResponse[tmp_idx]

			# basic QQ & residual v. predicted plot
				 	png(filename = paste0(working_dir,"dharma_",u_species[i],"_qq.png"),width = 16, height = 9, units = "in",res=300)
						plot(tmp_residuals)
					dev.off()
				# residual v. predicted by model covariate
					png(filename = paste0(working_dir,"dharma_",u_species[i],"_resid_year.png"),width = 9, height = 9, units = "in",res=300)
						DHARMa::plotResiduals(tmp_residuals, form = as.factor(tmp_bfish$year))
					dev.off()
					png(filename = paste0(working_dir,"dharma_",u_species[i],"_resid_lon.png"),width = 9, height = 9, units = "in",res=300)
						DHARMa::plotResiduals(tmp_residuals, form = tmp_bfish$lon)
					dev.off()
					png(filename = paste0(working_dir,"dharma_",u_species[i],"_resid_lat.png"),width = 9, height = 9, units = "in",res=300)
						DHARMa::plotResiduals(tmp_residuals, form = tmp_bfish$lat)
					dev.off()
				# test for uniformity in residuals, overdispersion, outliers
				 	png(filename = paste0(working_dir,"dharma_",u_species[i],"_residual_tests.png"),width = 16, height = 9, units = "in",res=300)
						tmp_test = DHARMa::testResiduals(tmp_residuals)
					dev.off()
				# test for zero inflation
				 	png(filename = paste0(working_dir,"dharma_",u_species[i],"_test_zi.png"),width = 9, height = 9, units = "in",res=300)
						tmp_zinf = DHARMa::testZeroInflation(tmp_residuals)
					dev.off()
				# test for over-dispersion
				 	png(filename = paste0(working_dir,"dharma_",u_species[i],"_test_od.png"),width = 9, height = 9, units = "in",res=300)
						DHARMa::testDispersion(tmp_residuals,alternative="two.sided")
					dev.off()
				# test for spatial autocorrelation
				# DHARMa::testSpatialAutocorrelation needs unique locations
					tmp_bfish$spatial_group = as.numeric(as.factor(paste0(tmp_bfish$lon,"_",tmp_bfish$lat)))
					tmp_spatial_group_dt = data.table(spatial_group=tmp_bfish$spatial_group,x=tmp_bfish$lon,y=tmp_bfish$lat) %>%
								   .[,.(x=mean(x),y=mean(y)),by=spatial_group]
					tmp_residuals_spatial_group = DHARMa::recalculateResiduals(tmp_residuals, group = tmp_bfish$spatial_group)	
				 	png(filename = paste0(working_dir,"dharma_",u_species[i],"_test_spcorr.png"),width = 9, height = 9, units = "in",res=300)
						tmp_dharma_sp = DHARMa::testSpatialAutocorrelation(residuals_spatial_group, x = tmp_spatial_group_dt$x, y = tmp_spatial_group_dt$y)
						text(min(tmp_spatial_group_dt$x),min(tmp_spatial_group_dt$y),paste0("Moran's I p-value: ",round(tmp_dharma_sp$p.value,digits=2)),adj=c(0,0))
					dev.off()
				# test for temporal autocorrelation
					tmp_bfish$temporal_group = factor(tmp_bfish$year,levels=as.character(sort(unique(tmp_bfish$year))))
					tmp_residuals_temporal_group = DHARMa::recalculateResiduals(tmp_residuals, group = tmp_bfish$temporal_group)	
				 	png(filename = paste0(working_dir,"dharma_",u_species[i],"_test_tcorr.png"),width = 9, height = 9, units = "in",res=300)		
						tmp_ar = DHARMa::testTemporalAutocorrelation(tmp_residuals_temporal_group, time=levels(tmp_bfish$temporal_group))
					dev.off()


				residual_dt.list[[i]]$KS_stat=tmp_test$uniformity$statistic
				residual_dt.list[[i]]$KS_pvalue=tmp_test$uniformity$p.value
				residual_dt.list[[i]]$dispersion_stat=tmp_test$dispersion$statistic
				residual_dt.list[[i]]$dispersion_pvalue=tmp_test$dispersion$p.value
				residual_dt.list[[i]]$outlier_stat=tmp_test$outliers$statistic
				residual_dt.list[[i]]$outlier_pvalue=tmp_test$outliers$p.value
				residual_dt.list[[i]]$zinf_stat=tmp_zinf$statistic
				residual_dt.list[[i]]$zinf_pvalue=tmp_zinf$p.value				
				residual_dt.list[[i]]$moran_pvalue=tmp_dharma_sp$p.value
				residual_dt.list[[i]]$DW_stat=tmp_ar$statistic
				residual_dt.list[[i]]$DW_pvalue=tmp_ar$p.value

				# clean-up
					rm(list=c("tmp_zinf","tmp_test","tmp_ar","tmp_idx","tmp_bfish","tmp_residuals","tmp_spatial_group_dt","tmp_dharma_sp","tmp_residuals_spatial_group","tmp_residuals_temporal_group"))
		}
		

		residual_dt = rbind(residual_dt,rbindlist(residual_dt.list))
		fwrite(residual_dt,file=paste0(working_dir,"residual_dt.csv"))

#_____________________________________________________________________________________________________________________________
# 6) Make plots
# predicted density
		deep7_code_vec = tolower(c("ETCA","APRU","PRSI","HYQU","PRFI","PRZO","ETCO"))
		deep7_name_vec = c("Ehu", "Lehi", "Kalekale", "Hapu'upu'u", "'Opakapaka", "Gindai", "Onaga")
		

		# predict on annual psu table
			predict_knot = as.data.table(drop_units(fit$Report$D_gct)) %>%
				setnames(.,c("Site","Category","Time","value"),c("knot","species","year","density")) %>%
				.[,knot:=as.numeric(as.character(knot))] %>%
				.[,species:=factor(species,levels=c("prfi","etca","etco","prsi","przo","hyqu","apru"))] %>%
				.[,species_hw:=deep7_name_vec[match(species,deep7_code_vec)]] %>%
				.[,year:=as.numeric(as.character(year))] %>%
				.[,.(species,species_hw,year,knot,density)] %>%
				.[order(species,year,knot)]

				knot_coords = spatial_list$latlon_x[,2:1]
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

			spatial_residual_dt = data.table(PSU=bfish_df$psu,year=as.numeric(as.character(bfish_df$year)),species=bfish_df$species_cd,residual=residuals_dharma$scaledResiduals)

				psu_table$knot = fit$spatial_list$NN_Extrap$nn.idx
			predict_psu = merge(psu_table,predict_knot[,.(knot,species,species_hw,year,density)],by="knot",allow.cartesian=TRUE) %>%
						  .[,area_km2:=0.5^2] %>%
						  merge(.,omega_dt,by=c("knot","species")) %>%
						  merge(.,epsilon_dt,by=c("year","knot","species")) %>%
						  merge(.,component_dt,by=c("year","knot","species")) %>%
						  merge(.,spatial_residual_dt,by=c("PSU","year","species"),all.x=TRUE)

		save(predict_psu,file=paste0(working_dir,"predict_psu.RData"))

			p = copy(predict_psu) %>%
				.[,.(biomass=sum(density),lon_eqd=mean(lon_eqd),lat_eqd=mean(lat_eqd)),by=.(Island,PSU,knot,year)] %>%
				ggplot() + 
				facet_wrap(~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = biomass),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5) +
				ggthemes::theme_few(base_size=20) + 
				viridis::scale_color_viridis("Est. density",begin = 0.1,end = 0.8,direction = 1,option = "H",trans="log10")
			ggsave(filename=paste0("pred_dens_agg.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

			# plot random effects
				p = copy(predict_psu) %>%
					.[,.(knot,species,omega1,omega2)] %>%
					unique(.) %>%
					melt(.,id.vars=c("knot","species")) %>%
					.[,variable:=factor(variable,levels=c("omega1","omega2"),labels=c("Omega 1","Omega 2"))] %>%
					ggplot() +
					ggtitle("Spatial random effect") +
					facet_wrap(~variable) +
					xlab("Species") +
					ylab("Random effect") +
					geom_hline(yintercept=0) +
					geom_boxplot(aes(x=species,y=value,fill=species)) +
					ggthemes::theme_few(base_size=20) + 
					viridis::scale_fill_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
			ggsave(filename=paste0("pred_omega_agg.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)	
				
				p = copy(predict_psu) %>%
					.[,.(knot,year,species,epsilon1,epsilon2)] %>%
					unique(.) %>%
					melt(.,id.vars=c("knot","species","year")) %>%
					.[,variable:=factor(variable,levels=c("epsilon1","epsilon2"),labels=c("Epsilon 1","Epsilon 2"))] %>%
					ggplot() +
					ggtitle("Spatiotemporal random effect") +
					facet_wrap(~variable) +
					xlab("Species") +
					ylab("Random effect") +
					geom_hline(yintercept=0) +
					geom_boxplot(aes(x=species,y=value,fill=species)) +
					ggthemes::theme_few(base_size=20) + 
					viridis::scale_fill_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
			ggsave(filename=paste0("pred_epsilon_agg.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)				

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

			for(i in 1:length(u_species))
			{

				p = copy(predict_psu) %>%
				.[species==u_species[i]] %>%
				.[,.(biomass=sum(density),lon_eqd=mean(lon_eqd),lat_eqd=mean(lat_eqd)),by=.(Island,PSU,knot,year)] %>%
				ggplot() + 
				ggtitle(paste0(unique(predict_psu[species==u_species[i]]$species_hw)," (",toupper(u_species[i]),")")) +
				facet_wrap(~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = biomass),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5) +
				ggthemes::theme_few(base_size=20) + 
				viridis::scale_color_viridis("Est. density",begin = 0.1,end = 0.8,direction = 1,option = "H",trans="log10")
				ggsave(filename=paste0("pred_dens_",u_species[i],".png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	  			# omega
	  			p = copy(predict_psu) %>%
				.[species==u_species[i]] %>%
				.[,.(PSU,knot,lon_eqd,lat_eqd,omega1,omega2)] %>%
				unique(.) %>%
				melt(.,id.vars=c("PSU","knot","lon_eqd","lat_eqd")) %>%
				.[,variable:=factor(variable,levels=c("omega1","omega2"),labels=c("Omega 1","Omega 2"))] %>%
				ggplot() + 
				ggtitle(paste0(unique(predict_psu[species==u_species[i]]$species_hw)," (",toupper(u_species[i]),")")) +
				facet_wrap(~variable) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = value),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5) +
				ggthemes::theme_few(base_size=20) + 
				scale_color_gradient2("Spatial\nrandom\neffect",low = "blue",mid = "gray90",high ="red",trans = trim_tails(range = c(-3,3))) 
				ggsave(filename=paste0("pred_omega_",u_species[i],".png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	  			# epsilon
	  			p = copy(predict_psu) %>%
				.[species==u_species[i]] %>%
				.[,.(PSU,year,knot,lon_eqd,lat_eqd,epsilon1,epsilon2)] %>%
				unique(.) %>%
				melt(.,id.vars=c("PSU","year","knot","lon_eqd","lat_eqd")) %>%
				.[,variable:=factor(variable,levels=c("epsilon1","epsilon2"),labels=c("Epsilon 1","Epsilon 2"))] %>%
				ggplot() + 
				ggtitle(paste0(unique(predict_psu[species==u_species[i]]$species_hw)," (",toupper(u_species[i]),")")) +
				facet_grid(variable~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = value),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5) +
				ggthemes::theme_few(base_size=20) + 
				scale_color_gradient2("Spatiotemporal\nrandom\neffect",low = "blue",mid = "gray90",high ="red",trans = trim_tails(range = c(-3,3))) 
				ggsave(filename=paste0("pred_epsilon_",u_species[i],".png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	  			# encounter
	  			p = copy(predict_psu) %>%
				.[species==u_species[i]] %>%
				.[,.(PSU,year,knot,lon_eqd,lat_eqd,encounter_prob)] %>%
				unique(.) %>%
				ggplot() + 
				ggtitle(paste0(unique(predict_psu[species==u_species[i]]$species_hw)," (",toupper(u_species[i]),")")) +
				facet_wrap(~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = encounter_prob),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5) +
				ggthemes::theme_few(base_size=20) + 
				viridis::scale_color_viridis("Encounter\nprobability",begin = 0.1,end = 0.8,direction = 1,option = "H")
				ggsave(filename=paste0("pred_enc_",u_species[i],".png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
	  			
	  			# positive
	  			p = copy(predict_psu) %>%
				.[species==u_species[i]] %>%
				.[,.(PSU,year,knot,lon_eqd,lat_eqd,positive_catch)] %>%
				unique(.) %>%
				ggplot() + 
				ggtitle(paste0(unique(predict_psu[species==u_species[i]]$species_hw)," (",toupper(u_species[i]),")")) +
				facet_wrap(~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, color = positive_catch),size=0.05) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5) +
				ggthemes::theme_few(base_size=20) + 
				viridis::scale_color_viridis("Positive\ncatch",begin = 0.1,end = 0.8,direction = 1,option = "H",trans="log10")
				ggsave(filename=paste0("pred_pos_",u_species[i],".png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	  			p = copy(predict_psu) %>%
				.[species==u_species[i]] %>%
				.[,.(Island,PSU,year,knot,lon_eqd,lat_eqd,residual)] %>%
				.[!is.na(residual)] %>%
				unique(.) %>%
				ggplot() + 
				ggtitle(paste0(unique(predict_psu[species==u_species[i]]$species_hw)," (",toupper(u_species[i]),")")) +
				facet_wrap(Island~year,nrow=5,scales="free") +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, fill = residual),shape=21,color="white",size=0.75) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5) +
				ggthemes::theme_few() + 
				scale_fill_gradient2("PIT\nresidual",low = "blue",mid = "gray90",high ="red",midpoint=0.5) 
				ggsave(filename=paste0("pred_residual_psu_",u_species[i],".png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

	  			p = copy(predict_psu) %>%
				.[species==u_species[i]] %>%
				.[,.(PSU,year,knot,lon_eqd,lat_eqd,residual)] %>%
				.[!is.na(residual)] %>%
				.[,.(residual=mean(residual,na.rm=TRUE)),by=.(year,knot)] %>%
				merge(.,unique(predict_knot[,.(knot,lon_eqd,lat_eqd)]),by="knot") %>%
				unique(.) %>%
				ggplot() + 
				ggtitle(paste0(unique(predict_psu[species==u_species[i]]$species_hw)," (",toupper(u_species[i]),")")) +
				facet_wrap(~year) +
				xlab("Eastings") +
				ylab("Northings") +
				geom_point(aes(x = lon_eqd, y = lat_eqd, fill = residual),color="white",size=2.25,shape=21) +
				# geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5) +
				ggthemes::theme_few(base_size=20) + 
				scale_fill_gradient2("PIT\nresidual",low = "blue",mid = "gray90",high ="red",midpoint=0.5) 
				ggsave(filename=paste0("pred_residual_knot_",u_species[i],".png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
			}

	mv_dt = as.data.table(index$Table) %>%
			.[,Category:=factor(Category,levels=paste0("Category_",1:7),labels=c("prfi","etca","etco","prsi","przo","hyqu","apru"))] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="mv"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)] %>%
			.[,Estimate:=Estimate*2.20462262185] %>%
			  .[,Estimate:=Estimate/1000000] %>%
			  .[,l95:=exp(log(Estimate)-2*CV)] %>%
			  .[,u95:=exp(log(Estimate)+2*CV)]
	mv_agg_dt = as.data.table(index_agg$Table) %>%
			.[Category=="Category_7"] %>%
			.[,Category:="agg"] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="mv"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)] %>%
			.[,Estimate:=Estimate*2.20462262185] %>%
			  .[,Estimate:=Estimate/1000000] %>%
			  .[,l95:=exp(log(Estimate)-2*CV)] %>%
			  .[,u95:=exp(log(Estimate)+2*CV)]

index_dt = rbind(mv_dt,mv_agg_dt)

p = copy(index_dt) %>%
	  .[Category!="agg"] %>%
	ggplot() +
	ylim(0,NA) +
	ylab("Predicted biomass (millions lbs)") +
	xlab("Year") +
	facet_wrap(~Category,scales="free") +
	geom_hline(yintercept=0) +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Category),alpha=0.25) +
	geom_path(aes(x=Time,y=Estimate,group=Model,color=Category),size=1.5) +
	viridis::scale_color_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
	viridis::scale_fill_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
	theme_few(base_size=20)
	ggsave(filename=paste0("index_by_species.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)


p = copy(index_dt) %>%
	  .[Category!="agg"] %>%
	  .[,Estimate:=Estimate/mean(Estimate),by=Category] %>%
			  .[,l95:=exp(log(Estimate)-2*CV)] %>%
			  .[,u95:=exp(log(Estimate)+2*CV)] %>%	  
	ggplot() +
	ylim(0,NA) +
	ylab("Relative abundance") +
	xlab("Year") +
	facet_wrap(~Category,scales="free") +
	geom_hline(yintercept=0) +
	geom_hline(yintercept=1,linetype="dashed") +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Category),alpha=0.25) +
	geom_path(aes(x=Time,y=Estimate,group=Model,color=Category),size=1.5) +
	viridis::scale_color_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
	viridis::scale_fill_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
	theme_few(base_size=20)
	ggsave(filename=paste0("index_by_species_relative.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

p = copy(index_dt) %>%
	  .[Category!="agg"] %>%
	ggplot() +
	ylim(0,NA) +
	ylab("Predicted biomass (millions lbs)") +
	xlab("Year") +
	geom_hline(yintercept=0) +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Category,fill=Category),alpha=0.25) +
	geom_area(aes(x=Time,y=Estimate,group=Category,fill=Category)) +
	geom_path(data=index_dt[Category=="agg"],aes(x=Time,y=Estimate),color="black",size=1.25) +
	geom_point(data=index_dt[Category=="agg"],aes(x=Time,y=Estimate),color="black",size=4) +
	geom_segment(data=index_dt[Category=="agg"],aes(x=Time,xend=Time,y=l95,yend=u95),color="black",size=1.25) +
	# geom_ribbon(data=index_dt[Category=="agg"],aes(x=Time,ymin=l95,ymax=u95),color="black",linetype="dashed",fill=NA,size=1.25) +
	viridis::scale_color_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
	viridis::scale_fill_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
	theme_few(base_size=20)
	ggsave(filename=paste0("index_by_species_stack.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
