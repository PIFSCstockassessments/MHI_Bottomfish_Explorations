

# Nicholas Ducharme-Barth
# 05/02/2022
# Set up simple spatiotemporal model example using VAST
# 1) Bring data (including spatial data)
# 2) Set-up barrier mesh (including conversion to equal distant projection)
# 3) Fit two species (e.g. prfi,etca) model using delta-Gamma distribution, add functionality to sum biomass across categories
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
	working_dir = paste0(proj.dir,"VAST/model_runs/tutorial/dgamma_mv_sum/")
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
	bfish_df = bfish_combined_long_dt %>% .[species_cd %in% c("prfi","etco","etca")] %>% as.data.frame(.)

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
								 ObsModel=c(2,3))
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
	    	  					c_i = as.numeric(factor(bfish_df[,'species_cd'],levels=c("prfi","etca","etco")))-1,
	    	  					category_names=c("prfi","etca","etco"),
	    	  					Expansion_cz = cbind( c(0,3,3), c(0,0,1)),
	    	  					covariate_data = NULL,
	    	  					catchability_data = NULL,
          						working_dir = working_dir,
          						newtonsteps = 1,
          						# extrapolation list args
          						projargs=slot(crs_eqd,"projargs"),input_grid=input_grid,
          						extrapolation_list = Extrapolation_List,
          						# spatial list args
    	  						Method = "Barrier",anisotropic_mesh = mesh_inla,grid_size_LL = 0.5/110,Save_Results = FALSE,LON_intensity=intensity_loc[,1],LAT_intensity=intensity_loc[,2],
    	  						spatial_list = spatial_list,build_model=FALSE)
		
		# turn off estimation of variance for spatiotemporal random effect for 3rd category (species)
		# and fix to low value since goes to zero when freely estimated
		modified_map = fit_setup$tmb_list$Map
		# modified_map$L_omega2_z = factor(c(1,2,NA),levels=c(1,2))
		modified_map$L_epsilon2_z = factor(c(1,2,NA),levels=c(1,2))

		modified_parameters = fit_setup$tmb_list$Parameters
		# modified_parameters$L_omega2_z
		modified_parameters$L_epsilon2_z = c(1,1,0.00001)

		fit = fit_model( settings=settings,
 				   				Lon_i=bfish_df$lon,
    			   				Lat_i=bfish_df$lat,
          						t_i=as.integer(bfish_df$year),
          						b_i=bfish_df$weight_kg,
          						a_i=rep(pi * (0.02760333457^2),nrow(bfish_df)), # assumed area swept from the MOUSS camera converted to km2; Ault et al 2018
	    	  					c_i = as.numeric(factor(bfish_df[,'species_cd'],levels=c("prfi","etca","etco")))-1,
	    	  					category_names=c("prfi","etca","etco"),
	  						    Expansion_cz = cbind( c(0,3,3), c(0,0,1)),
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
    	  						Parameters = modified_parameters)

		# Sdreport = fit$parameter_estimates$SD
		# par_biascorrect = TMB:::as.list.sdreport( Sdreport, what="Est. (bias.correct)", report=TRUE )

		index = plot_biomass_index( fit,DirName=working_dir)

#_____________________________________________________________________________________________________________________________
# 5) Model diagnostics & plot output
	# need to simulate responses to calculate DHARMa style residuals
	# can use built in functions from sdmTMB package but DHARMa package has some built in tests/plots
	n_sim = 500
	sim_fit = matrix(NA,nrow=length(fit$data_list$b_i),ncol=n_sim)
	rs = 123
	for(i in 1:n_sim)
	{
		sim_fit[,i] = simulate_data( fit,type = 1,random_seed = list(rs+i,NULL)[[1+is.null(rs)]])$b_i
	}
	pred = fit$Report$D_i
	residuals_dharma =  DHARMa::createDHARMa(simulatedResponse = sim_fit,observedResponse = bfish_df$weight_kg,fittedPredictedResponse = pred)	

	# # DHARMa residuals: VAST internal
	# 	vast_internal_residuals_dharma = DHARMa::createDHARMa(simulatedResponse=sim_fit + 1e-10*array(rnorm(prod(dim(sim_fit))),dim=dim(sim_fit)),
	#         observedResponse=strip_units(fit$data_list$b_i) + 1e-10*rnorm(length(fit$data_list$b_i)),
	#         fittedPredictedResponse=strip_units(fit$Report$D_i),
	#         integer=FALSE)
	# 	prop_lessthan_i = apply( sim_fit<outer(fit$data_list$b_i,rep(1,n_sim)),
	#         MARGIN=1,
	#         FUN=mean )
	#     prop_lessthanorequalto_i = apply( sim_fit<=outer(fit$data_list$b_i,rep(1,n_sim)),
	#         MARGIN=1,
	#         FUN=mean )
	#     PIT_i = runif(min=prop_lessthan_i, max=prop_lessthanorequalto_i, n=length(prop_lessthan_i) )
	#       # cbind( "Difference"=dharmaRes$scaledResiduals - PIT_i, "PIT"=PIT_i, "Original"=dharmaRes$scaledResiduals, "b_i"=x$data_list$b_i )
	#     vast_internal_residuals_dharma$scaledResiduals = PIT_i
	
	# # DHARMa residuals: VAST 	
	# 	vast_residuals_dharma = summary( fit,what = "residuals",n_samples = 500,random_seed = rs)

	# 	# plot(residuals_dharma)
	# 	# x11()
	# 	# plot(vast_internal_residuals_dharma)
	# 	# x11()
	# 	# plot(vast_residuals_dharma)

# basic QQ & residual v. predicted plot
	 	png(filename = paste0(working_dir,"dharma_qq.png"),width = 16, height = 9, units = "in",res=300)
			plot(residuals_dharma)
		dev.off()
	# residual v. predicted by model covariate
		png(filename = paste0(working_dir,"dharma_resid_year.png"),width = 9, height = 9, units = "in",res=300)
			DHARMa::plotResiduals(residuals_dharma, form = as.factor(bfish_df$year))
		dev.off()
		png(filename = paste0(working_dir,"dharma_resid_lon.png"),width = 9, height = 9, units = "in",res=300)
			DHARMa::plotResiduals(residuals_dharma, form = bfish_df$lon)
		dev.off()
		png(filename = paste0(working_dir,"dharma_resid_lat.png"),width = 9, height = 9, units = "in",res=300)
			DHARMa::plotResiduals(residuals_dharma, form = bfish_df$lat)
		dev.off()
	# test for uniformity in residuals, overdispersion, outliers
	 	png(filename = paste0(working_dir,"dharma_residual_tests.png"),width = 16, height = 9, units = "in",res=300)
			DHARMa::testResiduals(residuals_dharma)
		dev.off()
	# test for zero inflation
	 	png(filename = paste0(working_dir,"dharma_test_zi.png"),width = 9, height = 9, units = "in",res=300)
			DHARMa::testZeroInflation(residuals_dharma)
		dev.off()
	# test for over-dispersion
	 	png(filename = paste0(working_dir,"dharma_test_od.png"),width = 9, height = 9, units = "in",res=300)
			DHARMa::testDispersion(residuals_dharma,alternative="greater")
		dev.off()
	# test for spatial autocorrelation
	# DHARMa::testSpatialAutocorrelation needs unique locations
		bfish_df$spatial_group = as.numeric(as.factor(paste0(bfish_df$lon,"_",bfish_df$lat)))
		spatial_group_dt = data.table(spatial_group=bfish_df$spatial_group,x=bfish_df$lon,y=bfish_df$lat) %>%
					   .[,.(x=mean(x),y=mean(y)),by=spatial_group]
		residuals_spatial_group = DHARMa::recalculateResiduals(residuals_dharma, group = bfish_df$spatial_group)	
	 	png(filename = paste0(working_dir,"dharma_test_spcorr.png"),width = 9, height = 9, units = "in",res=300)
			DHARMa::testSpatialAutocorrelation(residuals_spatial_group, x = spatial_group_dt$x, y = spatial_group_dt$y)
		dev.off()
	# test for temporal autocorrelation
		bfish_df$temporal_group = factor(bfish_df$year,levels=as.character(sort(unique(bfish_df$year))))
		residuals_temporal_group = DHARMa::recalculateResiduals(residuals_dharma, group = bfish_df$temporal_group)	
	 	png(filename = paste0(working_dir,"dharma_test_tcorr.png"),width = 9, height = 9, units = "in",res=300)		
			DHARMa::testTemporalAutocorrelation(residuals_temporal_group, time=levels(bfish_df$temporal_group))
		dev.off()

#_____________________________________________________________________________________________________________________________
# 6) Make plots
# covariate effects
	bfish_df$pred = fit$Report$D_i
	bfish_df$pred_1 = fit$Report$R1_i
	bfish_df$pred_2 = fit$Report$R2_i
	bfish_dt = as.data.table(bfish_df)
	# year
		# p = bfish_dt %>% .[,.(weight_kg,pred,year)] %>%
		# 				 # .[,year:=as.numeric(year)] %>%
		# 				 setnames(.,c("weight_kg","pred"),c("observed","expected")) %>%
		# 				 melt(.,id.vars="year") %>%
		# 	ggplot() + 
		# 	xlab("Year") +
		# 	ylab("Response") +
		# 	geom_boxplot(aes(x = year, y = value, fill=variable))

		p = bfish_dt %>% .[,.(weight_kg,pred,year)] %>%
						 .[,year:=as.numeric(year)] %>%
						 setnames(.,c("weight_kg","pred"),c("observed","expected")) %>%
						 melt(.,id.vars="year") %>%
						 .[,.(value=mean(value)),by=.(year,variable)] %>%
						 ggplot() + 
						 xlab("Year") +
						 ylab("Response") +
						 geom_path(aes(x = year, y = value, color=variable))
			ggsave(filename=paste0("nominal_predicted_index.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

		p = bfish_dt %>% .[,.(weight_kg,pred_1,year)] %>%
						 .[,weight_kg:=ifelse(weight_kg>0,1,0)] %>%
						 .[,year:=as.numeric(year)] %>%
						 setnames(.,c("weight_kg","pred_1"),c("observed","expected")) %>%
						 melt(.,id.vars="year") %>%
						 .[,.(value=mean(value)),by=.(year,variable)] %>%
						 ggplot() + 
						 xlab("Year") +
						 ylab("Response") +
						 geom_path(aes(x = year, y = value, color=variable))
			ggsave(filename=paste0("enc.nominal_predicted_index.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
		
		p = bfish_dt %>% .[,.(weight_kg,pred_2,year)] %>%
						 .[weight_kg>0] %>%
						 .[,year:=as.numeric(year)] %>%
						 setnames(.,c("weight_kg","pred_2"),c("observed","expected")) %>%
						 melt(.,id.vars="year") %>%
						 .[,.(value=mean(value)),by=.(year,variable)] %>%
						 ggplot() + 
						 xlab("Year") +
						 ylab("Response") +
						 geom_path(aes(x = year, y = value, color=variable))
			ggsave(filename=paste0("pos.nominal_predicted_index.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
	# lon
		p = bfish_dt %>% .[,.(weight_kg,pred,lon)] %>%
						 .[,lon:=floor(lon/0.25)*0.25] %>%
 						 setnames(.,c("weight_kg","pred"),c("observed","expected")) %>%
						 .[,.(expected=mean(expected),observed=mean(observed)),by=lon] %>%
 						 melt(.,id.vars="lon") %>%
						 ggplot() + 
						 xlab("Longitude") +
						 ylab("Response") +
						 geom_line(aes(x = lon, y = value, color=variable))
					ggsave(filename=paste0("response_lon.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
		p = bfish_dt %>% .[,.(weight_kg,pred_1,lon)] %>%
						 .[,weight_kg:=ifelse(weight_kg>0,1,0)] %>%
						 .[,lon:=floor(lon/0.25)*0.25] %>%
 						 setnames(.,c("weight_kg","pred_1"),c("observed","expected")) %>%
						 .[,.(expected=mean(expected),observed=mean(observed)),by=lon] %>%
 						 melt(.,id.vars="lon") %>%
						 ggplot() + 
						 xlab("Longitude") +
						 ylab("Response") +
						 geom_line(aes(x = lon, y = value, color=variable))
					ggsave(filename=paste0("enc.response_lon.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
		p = bfish_dt %>% .[,.(weight_kg,pred_2,lon)] %>%
						 .[weight_kg>0] %>%
						 .[,lon:=floor(lon/0.25)*0.25] %>%
 						 setnames(.,c("weight_kg","pred_2"),c("observed","expected")) %>%
						 .[,.(expected=mean(expected),observed=mean(observed)),by=lon] %>%
 						 melt(.,id.vars="lon") %>%
						 ggplot() + 
						 xlab("Longitude") +
						 ylab("Response") +
						 geom_line(aes(x = lon, y = value, color=variable))
					ggsave(filename=paste0("pos.response_lon.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

		p0 = bfish_dt %>% .[,.(weight_kg,pred_1,lon)] %>%
						 .[,lon:=floor(lon/0.25)*0.25] %>%
						 .[,.(pred=mean(pred_1)),by=lon] %>%
						 .[,pred:=pred/mean(pred)]

		p = bfish_dt %>% .[,.(year,lon)] %>%
						 .[,year:=as.numeric(year)] %>%
						 .[,lon:=floor(lon/0.25)*0.25] %>%				 
						 merge(.,p0) %>%
						 .[,.(pred=mean(pred),lon=mean(lon)),by=year] %>%
						 ggplot() + 
						 xlab("Year") +
						 ylab("Influence") +
						 geom_hline(yintercept=1,linetype="dotted") +
						 geom_line(aes(x = year, y = pred))	+
						 geom_point(aes(x=year,y=pred,color=lon),size=3)					 
					ggsave(filename=paste0("enc.influ_lon.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

		p0 = bfish_dt %>% .[,.(weight_kg,pred_2,lon)] %>%
						 .[,lon:=floor(lon/0.25)*0.25] %>%
						 .[,.(pred=mean(pred_2)),by=lon] %>%
						 .[,pred:=pred/mean(pred)]

		p = bfish_dt %>% .[,.(year,lon)] %>%
						 .[,year:=as.numeric(year)] %>%
						 .[,lon:=floor(lon/0.25)*0.25] %>%				 
						 merge(.,p0) %>%
						 .[,.(pred=mean(pred),lon=mean(lon)),by=year] %>%
						 ggplot() + 
						 xlab("Year") +
						 ylab("Influence") +
						 geom_hline(yintercept=1,linetype="dotted") +
						 geom_line(aes(x = year, y = pred))	+
						 geom_point(aes(x=year,y=pred,color=lon),size=3)					 
					ggsave(filename=paste0("pos.influ_lon.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
	# lat
		p = bfish_dt %>% .[,.(weight_kg,pred,lat)] %>%
						 .[,lat:=floor(lat/0.25)*0.25] %>%
 						 setnames(.,c("weight_kg","pred"),c("observed","expected")) %>%
						 .[,.(expected=mean(expected),observed=mean(observed)),by=lat] %>%
 						 melt(.,id.vars="lat") %>%
						 ggplot() + 
						 xlab("Latitude") +
						 ylab("Response") +
						 geom_line(aes(x = lat, y = value, color=variable))
					ggsave(filename=paste0("response_lat.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
		p = bfish_dt %>% .[,.(weight_kg,pred_1,lat)] %>%
						 .[,weight_kg:=ifelse(weight_kg>0,1,0)] %>%
						 .[,lat:=floor(lat/0.25)*0.25] %>%
 						 setnames(.,c("weight_kg","pred_1"),c("observed","expected")) %>%
						 .[,.(expected=mean(expected),observed=mean(observed)),by=lat] %>%
 						 melt(.,id.vars="lat") %>%
						 ggplot() + 
						 xlab("Latitude") +
						 ylab("Response") +
						 geom_line(aes(x = lat, y = value, color=variable))
					ggsave(filename=paste0("enc.response_lat.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
		p = bfish_dt %>% .[,.(weight_kg,pred_2,lat)] %>%
						 .[weight_kg>0] %>%
						 .[,lat:=floor(lat/0.25)*0.25] %>%
 						 setnames(.,c("weight_kg","pred_2"),c("observed","expected")) %>%
						 .[,.(expected=mean(expected),observed=mean(observed)),by=lat] %>%
 						 melt(.,id.vars="lat") %>%
						 ggplot() + 
						 xlab("Latitude") +
						 ylab("Response") +
						 geom_line(aes(x = lat, y = value, color=variable))
					ggsave(filename=paste0("pos.response_lat.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

		p0 = bfish_dt %>% .[,.(weight_kg,pred_1,lat)] %>%
						 .[,lat:=floor(lat/0.25)*0.25] %>%
						 .[,.(pred=mean(pred_1)),by=lat] %>%
						 .[,pred:=pred/mean(pred)]

		p = bfish_dt %>% .[,.(year,lat)] %>%
						 .[,year:=as.numeric(year)] %>%
						 .[,lat:=floor(lat/0.25)*0.25] %>%				 
						 merge(.,p0) %>%
						 .[,.(pred=mean(pred),lat=mean(lat)),by=year] %>%
						 ggplot() + 
						 xlab("Year") +
						 ylab("Influence") +
						 geom_hline(yintercept=1,linetype="dotted") +
						 geom_line(aes(x = year, y = pred))	+
						 geom_point(aes(x=year,y=pred,color=lat),size=3)					 
					ggsave(filename=paste0("enc.influ_lat.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

		p0 = bfish_dt %>% .[,.(weight_kg,pred_2,lat)] %>%
						 .[,lat:=floor(lat/0.25)*0.25] %>%
						 .[,.(pred=mean(pred_2)),by=lat] %>%
						 .[,pred:=pred/mean(pred)]

		p = bfish_dt %>% .[,.(year,lat)] %>%
						 .[,year:=as.numeric(year)] %>%
						 .[,lat:=floor(lat/0.25)*0.25] %>%				 
						 merge(.,p0) %>%
						 .[,.(pred=mean(pred),lat=mean(lat)),by=year] %>%
						 ggplot() + 
						 xlab("Year") +
						 ylab("Influence") +
						 geom_hline(yintercept=1,linetype="dotted") +
						 geom_line(aes(x = year, y = pred))	+
						 geom_point(aes(x=year,y=pred,color=lat),size=3)					 
					ggsave(filename=paste0("pos.influ_lat.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

# predicted density
		d_dt = as.data.table(drop_units(fit$Report$D_gct)) %>%
			   .[,.(Site,Time,value)] %>%
			   setnames(.,c("Site","Time"),c("knot","year")) %>%
			   .[,knot:=as.numeric(knot)] %>%
			   .[,year:=as.numeric(year)] 

		g_dt = as.data.table(Extrapolation_List$Data_Extrap)
		g_dt$knot = spatial_list$NN_Extrap$nn.idx

		g_dt.list = lapply(unique(bfish_df$year),function(x,dt){tmp=copy(dt)%>%.[,year:=x];return(tmp)},dt=g_dt)
		g_annual_dt = rbindlist(g_dt.list) %>% .[,year:=as.numeric(year)] %>% merge(.,d_dt,by=c("year","knot"))
		
		p = g_annual_dt %>%
			ggplot() + 
			facet_wrap(~year) +
			xlab("Longitude") +
			ylab("Latitude") +
			geom_point(aes(x=Lon,y=Lat,color=value)) +	
			viridis::scale_color_viridis("Response",begin = 0.1,end = 0.8,direction = 1,option = "H",trans="log10")
		ggsave(filename=paste0("pred_bio.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

		d_dt = as.data.table(fit$Report$R1_gct) %>%
			   .[,.(Site,Time,value)] %>%
			   setnames(.,c("Site","Time"),c("knot","year")) %>%
			   .[,knot:=as.numeric(knot)] %>%
			   .[,year:=as.numeric(year)] 

		g_dt = as.data.table(Extrapolation_List$Data_Extrap)
		g_dt$knot = spatial_list$NN_Extrap$nn.idx

		g_dt.list = lapply(unique(bfish_df$year),function(x,dt){tmp=copy(dt)%>%.[,year:=x];return(tmp)},dt=g_dt)
		g_annual_dt = rbindlist(g_dt.list) %>% .[,year:=as.numeric(year)] %>% merge(.,d_dt,by=c("year","knot"))
		
		p = g_annual_dt %>%
			ggplot() + 
			facet_wrap(~year) +
			xlab("Longitude") +
			ylab("Latitude") +
			geom_point(aes(x=Lon,y=Lat,color=value)) +	
			viridis::scale_color_viridis("Response",begin = 0.1,end = 0.8,direction = 1,option = "H",trans="log10")
		ggsave(filename=paste0("pred_enc.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
		
		d_dt = as.data.table(fit$Report$R2_gct) %>%
			   .[,.(Site,Time,value)] %>%
			   setnames(.,c("Site","Time"),c("knot","year")) %>%
			   .[,knot:=as.numeric(knot)] %>%
			   .[,year:=as.numeric(year)] 

		g_dt = as.data.table(Extrapolation_List$Data_Extrap)
		g_dt$knot = spatial_list$NN_Extrap$nn.idx

		g_dt.list = lapply(unique(bfish_df$year),function(x,dt){tmp=copy(dt)%>%.[,year:=x];return(tmp)},dt=g_dt)
		g_annual_dt = rbindlist(g_dt.list) %>% .[,year:=as.numeric(year)] %>% merge(.,d_dt,by=c("year","knot"))
		
		p = g_annual_dt %>%
			ggplot() + 
			facet_wrap(~year) +
			xlab("Longitude") +
			ylab("Latitude") +
			geom_point(aes(x=Lon,y=Lat,color=value)) +	
			viridis::scale_color_viridis("Response",begin = 0.1,end = 0.8,direction = 1,option = "H",trans="log10")
		ggsave(filename=paste0("pred_pos.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)
# index
		p = as.data.table(index$log_Index_ctl) %>%
				   dcast(.,V1+V2+V3~V4) %>%
				   .[,V2:=as.numeric(V2)+2016] %>%
				   .[V3==1] %>%
				   setnames(.,c("V2","Estimate","Std. Error"),c("year","index","cv")) %>%
				   .[,.(year,index,cv)] %>%
				   .[,upper:=1.96*cv+index] %>%
				   .[,lower:=-1.96*cv+index] %>%
				   .[,index:=exp(index)] %>%
				   .[,upper:=exp(upper)] %>%
				   .[,lower:=exp(lower)] %>%
			ggplot() + 
			xlab("Year") +
			ylab("Total biomass") +
			geom_ribbon(aes(x=year,ymin=lower,ymax=upper)) +
			geom_path(aes(x=year,y=index),color="white")
		ggsave(filename=paste0("std_idx.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)


