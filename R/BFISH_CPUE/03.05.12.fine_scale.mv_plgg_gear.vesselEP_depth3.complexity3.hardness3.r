
# Nicholas Ducharme-Barth
# 12/22/2022
# Spatiotemporal model using VAST

# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(VAST)
	library(ggplot2)
	library(ggthemes)
	library(splines)
	# rgdal functions sourced directly
	# sp functions sourced directly
	# DHARMa functions sourced directly
	# sf functions sourced directly
	# ggthemes functions sourced directly
#_____________________________________________________________________________________________________________________________
# specify (coarse) model configuration
    data_flag = 2021
    link_function = "plgg" # poisson-link delta-gamma
    species = "mv"
    data_treatment = "05"
    catchability_covariates = "gear.vesselEP" # vanilla
    abundance_covariates = "depth3.complexity3.hardness3.SCALED" # vanilla
    lehi_filter = TRUE
    km_cutoff = 7.5 # make this smaller to increase the spatial resolution of the model
    fine_scale = TRUE
    bias_correct = TRUE
	residual_type = "pit" # other option is one step ahead (osa) which is sloooooow (~30 minutes)
	xval = "noxval" # "xval" # warning xval also appears to be quite slow... (~30 minutes)

    # can bring in spatial data from an existing model if spatial structure of data is identical
    load_spatial = TRUE
    
    model_name = paste(data_flag,
                link_function,
                species,data_treatment,
                catchability_covariates,
                abundance_covariates,
                lehi_filter,
                km_cutoff,
                fine_scale,
                bias_correct,
				residual_type,
                xval,
                sep="_")
    
#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"
	working_dir = paste0(proj.dir,"VAST/model_runs/",as.character(format(Sys.time(),format="%Y-%m-%d")),"/",model_name,"/")
	dir.create(working_dir,recursive=TRUE)

    load_spatial_path = paste0(proj.dir,"VAST/model_runs/2023-01-05/2021_pldg_mv_05_gearSP_v_TRUE_7.5_TRUE_TRUE_pit_noxval/")
	# xval path
    # load_xval_path = paste0(proj.dir,"VAST/xval_data/2021_single_05_TRUE_7.5_FALSE_10_123/")

#_____________________________________________________________________________________________________________________________
# 1) bring in data
	load(file=paste0(proj.dir,"Data/",data_flag,"_",data_treatment,".bfish_combined_long_dt.RData"))

	# subset to species
    if(species == "mv")
    { 
        target_species = c("prfi","etca","etco","prsi","przo","hyqu","apru")
    } else{
        target_species = species
    }
	bfish_df = bfish_combined_long_dt %>% .[species_cd %in% target_species] %>% as.data.frame(.)
	
	# remove sample with large lehi observation
    if(lehi_filter)
    {
        bfish_df =  subset(bfish_df,design_sampling_unit!="2021_Fall_32293")
    }
	
    # define catchability section
        q_data = bfish_df[,c("year","species_cd","gear_type")]
        colnames(q_data)[2] = "category"
        q_data$category = factor(q_data[,'category'],levels=c(target_species))
        q_data$gear_type = factor(q_data[,'gear_type'])

        q1_formula = ~ gear_type
        q2_formula = ~ gear_type

        continuous_q_variables = c()

	# needed to define spatial domain and for predicting on to create index
	psu_table = fread(paste0(proj.dir,"Data/BFISH PSU lookup table.csv")) %>%
				.[,.(PSU,Island,lon_deg,lat_deg,STRATA,STRATA_2020,Depth_MEDIAN_m,med_slp,med_acr,BS_pct_over_136j,pctHB,pctHS)] %>%
				.[,substrate:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][1])] %>%
				.[,slope:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][2])]
	# needed for plotting and defining the barrier feature
	hi_coast = rgdal::readOGR(dsn = paste0(proj.dir,"Data/GIS/Coastline"), layer = "Coastline")
	hi_coast_sf = sf::st_as_sf(hi_coast)

	# define abundance section
		ab_df = copy(psu_table) %>%
                        .[STRATA=="SB_H_S",STRATA:="SB_A_S"] %>%
						.[,.(lat_deg,lon_deg,STRATA,Island,Depth_MEDIAN_m,med_acr,BS_pct_over_136j)] %>%
						setnames(.,c("lat_deg","lon_deg","STRATA","Island","Depth_MEDIAN_m","med_acr","BS_pct_over_136j"),c("Lat","Lon","strata","island","depth","complexity","hardness")) %>%
						.[complexity>300,complexity:=300] %>%
                        .[is.na(complexity),complexity:=12.20] %>%
                        .[is.na(hardness),hardness:=0.094888] %>%
						.[,Year:=NA] %>%
						.[,depth_sc:=scale(depth)] %>%
						.[,complexity_sc:=scale(complexity)] %>%
						.[,hardness_sc:=scale(hardness)] %>%
						as.data.frame(.)
		
		ab1_formula = ~ bs(depth_sc,df=3) + bs(complexity_sc,df=3) + bs(hardness_sc,df=3)
        ab2_formula = ~ 0

        continuous_ab_variables = c("depth_sc","complexity_sc","hardness_sc")

	

#_____________________________________________________________________________________________________________________________
# 2) set-up spatial domain
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
    if(load_spatial)
    {
        load(file=paste0(load_spatial_path,"Extrapolation_List.RData"))
    } else {
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
    }


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
		mesh_inla = INLA::inla.mesh.2d(loc = psu_sp_eqd,boundary = mesh_boundary,max.n.strict=c(500,16),cutoff=km_cutoff)

		# define knot locations for spatial_list
			mesh_coords_eqd = sp::SpatialPoints(mesh_inla$loc[,1:2])
			sp::proj4string(mesh_coords_eqd) = crs_eqd
			mesh_coords_ll = sp::spTransform(mesh_coords_eqd, crs_ll)
			intensity_loc = mesh_coords_ll@coords 
		
        if(load_spatial)
        {
            load(file=paste0(load_spatial_path,"spatial_list.RData"))
        } else {
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
    					   map_data=hi_coast); gc()
        }


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

		save(spatial_list,file=paste0(working_dir,"spatial_list.RData"))
		save(Extrapolation_List,file=paste0(working_dir,"Extrapolation_List.RData"))
		gc()
#_____________________________________________________________________________________________________________________________
# 3) fit a basic VAST model
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
            xlim(0,NA) +
            ylim(0,NA) +
			xlab("log(Mean)") +
			ylab("log(Variance)") +
			geom_abline(slope=1,intercept=meanvar_lm$coefficients[1],linetype="dashed",color="gray70") +
			geom_abline(slope=2,intercept=meanvar_lm$coefficients[1],linetype="dashed",color="gray70") +
			geom_abline(slope=3,intercept=meanvar_lm$coefficients[1],linetype="dashed",color="gray70") +
			geom_abline(slope=meanvar_lm$coefficients[2],intercept=meanvar_lm$coefficients[1],linetype="dashed",color="black",linewidth=2) +			
			geom_point(aes(x=mean,y=var,fill=species_cd,size=N),shape=21) +
			theme_few(base_size=20) +
			geom_text(data=text_dt,aes(x=x,y=y,label=label),size=7) +
     		viridis::scale_fill_viridis("Species\ncode",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
			ggsave(filename=paste0("mean_var.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

    if(link_function == "pldg")
    {
        obs_model = c(2,4)
    } else if(link_function == "lgdg"){
        obs_model = c(2,3)
    } else if(link_function == "pldl"){
        obs_model = c(4,4)
    } else if(link_function == "plgg"){
        obs_model = c(9,4)
    }

	# make settings
		settings = make_settings(n_x=spatial_list$n_x,
								 Region="user",
								 fine_scale=fine_scale,
								 purpose="index2",
								 use_anisotropy = FALSE,
								 FieldConfig=matrix( rep("IID",6), ncol=2, nrow=3, dimnames=list(c("Omega","Epsilon","Beta"),c("Component_1","Component_2")) ),
								 Options=c("treat_nonencounter_as_zero"=TRUE ),
								 bias.correct=FALSE,
								 max_cells=Inf,
								 ObsModel=obs_model)
		settings$grid_size_km = 0.5
		settings$OverdispersionConfig = c("eta1"=1, "eta2"=1)
		settings$Options = c( settings$Options, "range_fraction"=0.01 )

	# set-up model
		dir.create(paste0(working_dir,"setup/"),recursive=TRUE)
		fit_setup = fit_model( settings=settings,
 				   				Lon_i=bfish_df$lon,
    			   				Lat_i=bfish_df$lat,
          						t_i=as.integer(bfish_df$year),
          						b_i=bfish_df$weight_kg,
          						a_i=rep(pi * (0.02760333457^2),nrow(bfish_df)), # assumed area swept from the MOUSS camera converted to km2; Ault et al 2018
	    	  					c_i = as.numeric(factor(bfish_df[,'species_cd'],levels=c(target_species)))-1,
	    	  					v_i = as.numeric(factor(bfish_df[,'platform']))-1,
								category_names=c(target_species),
	    	  					covariate_data = ab_df,
								X1_formula = ab1_formula,
								X2_formula = ab2_formula,
	    	  					catchability_data = q_data,
                                Q1_formula = q1_formula,
                                Q2_formula = q2_formula,
          						working_dir = paste0(working_dir,"setup/"),
          						newtonsteps = 0,
          						# extrapolation list args
          						projargs=slot(crs_eqd,"projargs"),input_grid=input_grid,
          						extrapolation_list = Extrapolation_List,
          						# spatial list args
    	  						Method = "Barrier",anisotropic_mesh = mesh_inla,grid_size_LL = 0.5/110,Save_Results = FALSE,LON_intensity=intensity_loc[,1],LAT_intensity=intensity_loc[,2],
    	  						spatial_list = spatial_list,build_model=TRUE,test_fit=FALSE); gc()
		
		# turn off estimation of components if poorly estimated
        # are mean estimates for L_ components going to zero
        # or are any variances blowing up
		if(length(fit_setup$parameter_estimates)==2)
		{
			fit_setup$parameter_estimates$opt
			# suggestions for problematic parameters
			which(abs(fit_setup$parameter_estimates$opt$par[which(names(fit_setup$parameter_estimates$opt$par)=="L_omega1_z")])<1e-3|abs(fit_setup$parameter_estimates$opt$par[which(names(fit_setup$parameter_estimates$opt$par)=="L_omega1_z")])>1.5e1)
			which(abs(fit_setup$parameter_estimates$opt$par[which(names(fit_setup$parameter_estimates$opt$par)=="L_epsilon1_z")])<1e-3|abs(fit_setup$parameter_estimates$opt$par[which(names(fit_setup$parameter_estimates$opt$par)=="L_epsilon1_z")])>1.5e1)
			which(abs(fit_setup$parameter_estimates$opt$par[which(names(fit_setup$parameter_estimates$opt$par)=="L_omega2_z")])<1e-3|abs(fit_setup$parameter_estimates$opt$par[which(names(fit_setup$parameter_estimates$opt$par)=="L_omega2_z")])>1.5e1)
			which(abs(fit_setup$parameter_estimates$opt$par[which(names(fit_setup$parameter_estimates$opt$par)=="L_epsilon2_z")])<1e-3|abs(fit_setup$parameter_estimates$opt$par[which(names(fit_setup$parameter_estimates$opt$par)=="L_epsilon2_z")])>1.5e1)

			modified_map = fit_setup$tmb_list$Map
			omega1_map = c(1,2,3,4,5,NA,6)
			epsilon1_map = c(1,2,NA,3,NA,4,5)
			omega2_map = c(1,2,3,NA,4,5,6)
			epsilon2_map = c(NA,1,2,3,NA,5,6)
			modified_map$L_omega1_z = factor(omega1_map,levels=1:max(omega1_map,na.rm=TRUE))
			modified_map$L_epsilon1_z = factor(epsilon1_map,levels=1:max(epsilon1_map,na.rm=TRUE))
			modified_map$L_omega2_z = factor(omega2_map,levels=1:max(omega2_map,na.rm=TRUE))
			modified_map$L_epsilon2_z = factor(epsilon2_map,levels=1:max(epsilon2_map,na.rm=TRUE))

			modified_parameters = fit_setup$tmb_list$Obj$env$parList()
			modified_parameters$L_omega1_z[which(is.na(omega1_map))] = 1e-8
			modified_parameters$L_epsilon1_z[which(is.na(epsilon1_map))] = 1e-8
			modified_parameters$L_omega2_z[which(is.na(omega2_map))] = 1e-8
			modified_parameters$L_epsilon2_z[which(is.na(epsilon2_map))] = 1e-8

		} else {
			fit_setup$parameter_estimates
			# suggestions for problematic parameters
			which(abs(fit_setup$ParHat$L_omega1_z)<1e-3|abs(fit_setup$ParHat$L_omega1_z)>1.5e1)
			which(abs(fit_setup$ParHat$L_epsilon1_z)<1e-3|abs(fit_setup$ParHat$L_epsilon1_z)>1.5e1)
			which(abs(fit_setup$ParHat$L_omega2_z)<1e-3|abs(fit_setup$ParHat$L_omega2_z)>1.5e1)
			which(abs(fit_setup$ParHat$L_epsilon2_z)<1e-3|abs(fit_setup$ParHat$L_epsilon2_z)>1.5e1)

			modified_map = fit_setup$tmb_list$Map
			omega1_map = c(1,2,3,4,5,NA,6)
			epsilon1_map = c(1,2,NA,3,NA,4,5)
			omega2_map = c(1,2,3,NA,4,5,6)
			epsilon2_map = c(NA,1,2,3,NA,5,6)
			modified_map$L_omega1_z = factor(omega1_map,levels=1:max(omega1_map,na.rm=TRUE))
			modified_map$L_epsilon1_z = factor(epsilon1_map,levels=1:max(epsilon1_map,na.rm=TRUE))
			modified_map$L_omega2_z = factor(omega2_map,levels=1:max(omega2_map,na.rm=TRUE))
			modified_map$L_epsilon2_z = factor(epsilon2_map,levels=1:max(epsilon2_map,na.rm=TRUE))

			modified_parameters = fit_setup$ParHat
			modified_parameters$L_omega1_z[which(is.na(omega1_map))] = 1e-8
			modified_parameters$L_epsilon1_z[which(is.na(epsilon1_map))] = 1e-8
			modified_parameters$L_omega2_z[which(is.na(omega2_map))] = 1e-8
			modified_parameters$L_epsilon2_z[which(is.na(epsilon2_map))] = 1e-8
		}

		fit = fit_model( settings=settings,
 				   				Lon_i=bfish_df$lon,
    			   				Lat_i=bfish_df$lat,
          						t_i=as.integer(bfish_df$year),
          						b_i=bfish_df$weight_kg,
          						a_i=rep(pi * (0.02760333457^2),nrow(bfish_df)), # assumed area swept from the MOUSS camera converted to km2; Ault et al 2018
	    	  					c_i = as.numeric(factor(bfish_df[,'species_cd'],levels=c(target_species)))-1,
	    	  					v_i = as.numeric(factor(bfish_df[,'platform']))-1,
								category_names=c(target_species),
	    	  					covariate_data = ab_df,
								X1_formula = ab1_formula,
								X2_formula = ab2_formula,
	    	  					catchability_data = q_data,
                                Q1_formula = q1_formula,
                                Q2_formula = q2_formula,
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
    	  						test_fit=TRUE); gc()
		fit$parameter_estimates

		# re-fit model using Map and Parameters from fit in order to get aggregate index
		dir.create(paste0(working_dir,"agg/"),recursive=TRUE)
		fit_agg = fit_model( settings=settings,
 				   				Lon_i=bfish_df$lon,
    			   				Lat_i=bfish_df$lat,
          						t_i=as.integer(bfish_df$year),
          						b_i=bfish_df$weight_kg,
          						a_i=rep(pi * (0.02760333457^2),nrow(bfish_df)), # assumed area swept from the MOUSS camera converted to km2; Ault et al 2018
	    	  					c_i = as.numeric(factor(bfish_df[,'species_cd'],levels=c(target_species)))-1,
	    	  					v_i = as.numeric(factor(bfish_df[,'platform']))-1,
								category_names=c(target_species),
	    	  					Expansion_cz = cbind( c(0,3,3,3,3,3,3), c(0,0,1,2,3,4,5)),
	    	  					covariate_data = ab_df,
								X1_formula = ab1_formula,
								X2_formula = ab2_formula,
	    	  					catchability_data = q_data,
                                Q1_formula = q1_formula,
                                Q2_formula = q2_formula,
          						working_dir = paste0(working_dir,"agg/"),
          						newtonsteps = 1,
          						# extrapolation list args
          						projargs=slot(crs_eqd,"projargs"),input_grid=input_grid,
          						extrapolation_list = Extrapolation_List,
          						# spatial list args
    	  						Method = "Barrier",anisotropic_mesh = mesh_inla,grid_size_LL = 0.5/110,Save_Results = FALSE,LON_intensity=intensity_loc[,1],LAT_intensity=intensity_loc[,2],
    	  						spatial_list = spatial_list,
								Map = fit$tmb_list$Map,
    	  						Parameters = fit$ParHat); gc()
			fit_agg$parameter_estimates$time_for_run
		
		# bias correction
			if(bias_correct)
			{
				A = proc.time()
				Sdreport = apply_epsilon( fit,
					ADREPORT_name = "Index_ctl",
					eps_name = "eps_Index_ctl",
					inner.control = list(sparse=TRUE, lowrank=TRUE, trace=FALSE) )
				# agg
				Sdreport_agg = apply_epsilon( fit_agg,
					ADREPORT_name = "Index_ctl",
					eps_name = "eps_Index_ctl",
					inner.control = list(sparse=TRUE, lowrank=TRUE, trace=FALSE) )
				B = proc.time()
				B-A	
			} else {
				Sdreport = fit$parameter_estimates$SD
				# agg
				Sdreport_agg = fit_agg$parameter_estimates$SD
			}


		# save outputs
			index = plot_biomass_index( fit,DirName=working_dir)
			index_agg = plot_biomass_index( fit_agg,DirName=paste0(working_dir,"agg/"))

			save(fit,file=paste0(working_dir,"fit.RData"))
			save(Sdreport,file=paste0(working_dir,"Sdreport.RData"))
			save(fit_agg,file=paste0(working_dir,"fit_agg.RData"))
			save(Sdreport_agg,file=paste0(working_dir,"Sdreport_agg.RData"))

			bfish_df$pred_weight_kg = fit$Report$D_i
			bfish_dt = as.data.table(bfish_df)		
			save(bfish_dt,file=paste0(working_dir,"bfish_dt.RData"))

		# calculate RMSE
			rmse_agg_dt = copy(bfish_dt) %>%
				  .[,.(rmse=sqrt(mean((weight_kg-pred_weight_kg)^2)),nrmse=sqrt(mean((weight_kg-pred_weight_kg)^2))/sd(weight_kg),.N)] %>%
				  .[,type:="agg"]
			rmse_spec_dt = copy(bfish_dt) %>%
					.[,.(rmse=sqrt(mean((weight_kg-pred_weight_kg)^2)),nrmse=sqrt(mean((weight_kg-pred_weight_kg)^2))/sd(weight_kg),.N),by=species_cd] %>%
					setnames(.,"species_cd","type")
			rmse_year_dt = copy(bfish_dt) %>%
					.[,.(rmse=sqrt(mean((weight_kg-pred_weight_kg)^2)),nrmse=sqrt(mean((weight_kg-pred_weight_kg)^2))/sd(weight_kg),.N),by=year] %>%
					setnames(.,"year","type")
			rmse_island_dt = copy(bfish_dt) %>%
					.[,.(rmse=sqrt(mean((weight_kg-pred_weight_kg)^2)),nrmse=sqrt(mean((weight_kg-pred_weight_kg)^2))/sd(weight_kg),.N),by=island] %>%
					setnames(.,"island","type")
			rmse_gear_dt = copy(bfish_dt) %>%
					.[,.(rmse=sqrt(mean((weight_kg-pred_weight_kg)^2)),nrmse=sqrt(mean((weight_kg-pred_weight_kg)^2))/sd(weight_kg),.N),by=gear_type] %>%
					setnames(.,"gear_type","type")
			rmse_strata_dt = copy(bfish_dt) %>%
					.[,.(rmse=sqrt(mean((weight_kg-pred_weight_kg)^2)),nrmse=sqrt(mean((weight_kg-pred_weight_kg)^2))/sd(weight_kg),.N),by=strata] %>%
					setnames(.,"strata","type")				  
			rmse_dt = rbind(rmse_agg_dt,rmse_spec_dt,rmse_year_dt,rmse_island_dt,rmse_gear_dt,rmse_strata_dt)
			fwrite(rmse_dt,file=paste0(working_dir,"rmse_dt.csv"))
			gc()


#_____________________________________________________________________________________________________________________________
# 4) Model diagnostics & plot output
	n_sim = 500
	sim_fit = matrix(NA,nrow=length(fit$data_list$b_i),ncol=n_sim)
	rs = 123
	for(i in 1:n_sim)
	{
		sim_dat = simulate_data( fit,type = 1,random_seed = list(rs+i,NULL)[[1+is.null(rs)]])
		sim_fit[,i] = sim_dat$b_i
	}

	save(sim_fit,file=paste0(working_dir,"sim_fit.RData")); gc()
	
	pred = fit$Report$D_i

	# calculate residuals
		if(residual_type == "pit")
		{
			A = proc.time()
				residuals_dharma = summary(fit,what="residuals",type=1,n_samples=500,random_seed=rs)
			B = proc.time()
			B-A
		} else if(residual_type == "osa"){
			A = proc.time()
				osa = TMBhelper::oneStepPredict_deltaModel( obj = fit$tmb_list$Obj,
				observation.name = "b_i",
				method = "cdf",
				data.term.indicator = "keep",
				deltaSupport = 0,
				trace = TRUE )

				residuals_dharma = DHARMa::createDHARMa(simulatedResponse=sim_fit,
				observedResponse=bfish_df$weight_kg,
				fittedPredictedResponse=pred,
				integer=FALSE)
				residuals_dharma$scaledResiduals = pnorm(osa$residual)
			B = proc.time()
			B-A
		}
		gc()

	# plot diagnostics

			if(residual_type == "pit")
			{
				residuals_dharma$simulatedResponse = sim_fit
				residuals_dharma$observedResponse = bfish_df$weight_kg
			}
			if(sum(is.na(pred))>0)
			{
				pred[which(is.na(pred))] = rowMeans(sim_fit[which(is.na(pred)),])
				residuals_dharma$fittedPredictedResponse = pred
			}

			# basic QQ & residual v. predicted plot
				 	png(filename = paste0(working_dir,"dharma_agg_qq.png"),width = 16, height = 9, units = "in",res=300)
						plot(residuals_dharma)
					dev.off()
				# quantile test: residual v. predicted by model covariate
					png(filename = paste0(working_dir,"dharma_agg_quant_resid.png"),width = 9, height = 9, units = "in",res=300)
						tmp_quant = DHARMa::testQuantiles(residuals_dharma)
					dev.off()
					png(filename = paste0(working_dir,"dharma_agg_quant_resid_year.png"),width = 9, height = 9, units = "in",res=300)
						tmp_quant_year = DHARMa::testQuantiles(residuals_dharma,predictor=as.factor(bfish_df$year))
						DHARMa::plotResiduals(residuals_dharma, form = as.factor(bfish_df$year))
					dev.off()
					png(filename = paste0(working_dir,"dharma_agg_quant_resid_lon.png"),width = 9, height = 9, units = "in",res=300)
						tmp_quant_lon = DHARMa::testQuantiles(residuals_dharma,predictor=bfish_df$lon)
					dev.off()
					png(filename = paste0(working_dir,"dharma_agg_quant_resid_lat.png"),width = 9, height = 9, units = "in",res=300)
						tmp_quant_lat = DHARMa::testQuantiles(residuals_dharma,predictor=bfish_df$lat)
					dev.off()
					png(filename = paste0(working_dir,"dharma_agg_quant_resid_gear.png"),width = 9, height = 9, units = "in",res=300)
						tmp_quant_gear = DHARMa::testQuantiles(residuals_dharma,predictor=as.factor(bfish_df$gear_type))
					dev.off()
				# make residual plots by factor variable
					residual_factor_vec = c("island","strata","strata_2020",
											"substrate","slope","depth_strata","depth_strata_2020",
											"complexity","hardness","season","month","platform")
					for(v in 1:length(residual_factor_vec))
					{
						png(filename = paste0(working_dir,"dharma_agg_quant_resid_",residual_factor_vec[v],".png"),width = 9, height = 9, units = "in",res=300)
							tmp_quant_factor = DHARMa::testQuantiles(residuals_dharma,predictor=as.factor(as.character(bfish_df[,residual_factor_vec[v]])))
						dev.off()
					}
					residual_continuous_vec = c("depth","time","jd","lunar_phase")
					residual_continuous_vec_pvalue = rep(NA,length(residual_continuous_vec))
					for(v in 1:length(residual_continuous_vec))
					{
						png(filename = paste0(working_dir,"dharma_agg_quant_resid_",residual_continuous_vec[v],".png"),width = 9, height = 9, units = "in",res=300)
							tmp_quant_continuous = DHARMa::testQuantiles(residuals_dharma,predictor=bfish_df[,residual_continuous_vec[v]])
						dev.off()
						residual_continuous_vec_pvalue[v] = tmp_quant_continuous$p.value
						rm(list=c("tmp_quant_continuous"))
					}
				# test for uniformity in residuals, overdispersion, outliers
				 	png(filename = paste0(working_dir,"dharma_agg_residual_tests.png"),width = 16, height = 9, units = "in",res=300)
						tmp_test = DHARMa::testResiduals(residuals_dharma)
					dev.off()
				# test for zero inflation
				 	png(filename = paste0(working_dir,"dharma_agg_test_zi.png"),width = 9, height = 9, units = "in",res=300)
						tmp_zinf = DHARMa::testZeroInflation(residuals_dharma)
					dev.off()
				# test for over-dispersion
				 	png(filename = paste0(working_dir,"dharma_agg_test_od.png"),width = 9, height = 9, units = "in",res=300)
						DHARMa::testDispersion(residuals_dharma,alternative="two.sided")
					dev.off()
				# test for spatial autocorrelation
				# DHARMa::testSpatialAutocorrelation needs unique locations
					bfish_df$spatial_group = as.numeric(as.factor(paste0(bfish_df$lon,"_",bfish_df$lat)))
					tmp_spatial_group_dt = data.table(spatial_group=bfish_df$spatial_group,x=bfish_df$lon,y=bfish_df$lat) %>%
								   .[,.(x=mean(x),y=mean(y)),by=spatial_group]
					residuals_spatial_group = DHARMa::recalculateResiduals(residuals_dharma, group = bfish_df$spatial_group)	
				 	png(filename = paste0(working_dir,"dharma_agg_test_spcorr.png"),width = 9, height = 9, units = "in",res=300)
						tmp_dharma_sp = DHARMa::testSpatialAutocorrelation(residuals_spatial_group, x = tmp_spatial_group_dt$x, y = tmp_spatial_group_dt$y)
						text(min(tmp_spatial_group_dt$x),min(tmp_spatial_group_dt$y),paste0("Moran's I p-value: ",round(tmp_dharma_sp$p.value,digits=2)),adj=c(0,0))
					dev.off()
				# test for temporal autocorrelation
					bfish_df$temporal_group = factor(bfish_df$year,levels=as.character(sort(unique(bfish_df$year))))
					residuals_temporal_group = DHARMa::recalculateResiduals(residuals_dharma, group = bfish_df$temporal_group)	
				 	png(filename = paste0(working_dir,"dharma_agg_test_tcorr.png"),width = 9, height = 9, units = "in",res=300)		
						tmp_ar = DHARMa::testTemporalAutocorrelation(residuals_temporal_group, time=levels(bfish_df$temporal_group))
					dev.off()

				residual_agg_dt = data.table(type="agg",
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
							quant_lunar_phase_pvalue=residual_continuous_vec_pvalue[4])

				# clean-up
					rm(list=c("tmp_spatial_group_dt","tmp_zinf","tmp_test","tmp_ar","tmp_dharma_sp","tmp_quant","tmp_quant_lon","tmp_quant_lat"))
		
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
				# quantile test: residual v. predicted by model covariate
					png(filename = paste0(working_dir,"dharma_",u_species[i],"_quant_resid.png"),width = 9, height = 9, units = "in",res=300)
						tmp_quant = DHARMa::testQuantiles(tmp_residuals)
					dev.off()
					png(filename = paste0(working_dir,"dharma_",u_species[i],"_quant_resid_year.png"),width = 9, height = 9, units = "in",res=300)
						tmp_quant_year = DHARMa::testQuantiles(tmp_residuals,predictor=as.factor(tmp_bfish$year))
						DHARMa::plotResiduals(tmp_residuals, form = as.factor(tmp_bfish$year))
					dev.off()
					png(filename = paste0(working_dir,"dharma_",u_species[i],"_quant_resid_lon.png"),width = 9, height = 9, units = "in",res=300)
						tmp_quant_lon = DHARMa::testQuantiles(tmp_residuals,predictor=tmp_bfish$lon)
					dev.off()
					png(filename = paste0(working_dir,"dharma_",u_species[i],"_quant_resid_lat.png"),width = 9, height = 9, units = "in",res=300)
						tmp_quant_lat = DHARMa::testQuantiles(tmp_residuals,predictor=tmp_bfish$lat)
					dev.off()
                    png(filename = paste0(working_dir,"dharma_",u_species[i],"_quant_resid_gear.png"),width = 9, height = 9, units = "in",res=300)
						tmp_quant_gear = DHARMa::testQuantiles(tmp_residuals,predictor=as.factor(tmp_bfish$gear_type))
					dev.off()
					for(v in 1:length(residual_factor_vec))
					{
                    	png(filename = paste0(working_dir,"dharma_",u_species[i],"_quant_resid_",residual_factor_vec[v],".png"),width = 9, height = 9, units = "in",res=300)
							tmp_quant_factor = DHARMa::testQuantiles(tmp_residuals,predictor=as.factor(as.character(tmp_bfish[,residual_factor_vec[v]])))
						dev.off()
					}
					tmp_continuous_vec_pvalue = rep(NA,length(residual_continuous_vec))
					for(v in 1:length(residual_continuous_vec))
					{
                    	png(filename = paste0(working_dir,"dharma_",u_species[i],"_quant_resid_",residual_continuous_vec[v],".png"),width = 9, height = 9, units = "in",res=300)
							tmp_quant_continuous = DHARMa::testQuantiles(tmp_residuals,predictor=tmp_bfish[,residual_continuous_vec[v]])
						dev.off()
						tmp_continuous_vec_pvalue[v] = tmp_quant_continuous$p.value
						rm(list=c("tmp_quant_continuous"))
					}
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
						tmp_dharma_sp = DHARMa::testSpatialAutocorrelation(tmp_residuals_spatial_group, x = tmp_spatial_group_dt$x, y = tmp_spatial_group_dt$y)
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
				residual_dt.list[[i]]$quant_pvalue=tmp_quant$p.value
				residual_dt.list[[i]]$quant_lon_pvalue=tmp_quant_lon$p.value
				residual_dt.list[[i]]$quant_lat_pvalue=tmp_quant_lat$p.value
				residual_dt.list[[i]]$quant_depth_pvalue=tmp_continuous_vec_pvalue[1]
				residual_dt.list[[i]]$quant_time_pvalue=tmp_continuous_vec_pvalue[2]
				residual_dt.list[[i]]$quant_jd_pvalue=tmp_continuous_vec_pvalue[3]
				residual_dt.list[[i]]$quant_lunar_pvalue=tmp_continuous_vec_pvalue[4]

				# clean-up
					rm(list=c("tmp_zinf","tmp_test","tmp_ar","tmp_idx","tmp_bfish","tmp_residuals","tmp_spatial_group_dt","tmp_dharma_sp","tmp_residuals_spatial_group","tmp_residuals_temporal_group","tmp_quant","tmp_quant_lon","tmp_quant_lat","tmp_continuous_vec_pvalue"))
		}
		
		residual_dt = rbind(residual_agg_dt,rbindlist(residual_dt.list))
		fwrite(residual_dt,file=paste0(working_dir,"residual_dt.csv"))
		gc()
#_____________________________________________________________________________________________________________________________
# 5) Make plots
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

				knot_coords = spatial_list$latlon_g[,2:1]
				knot_sp = sp::SpatialPoints(knot_coords)
				sp::proj4string(knot_sp) = crs_ll
				knot_sp_eqd = sp::spTransform(knot_sp, crs_eqd)
				knot_loc_dt = data.table(knot=1:nrow(knot_coords),lon=knot_sp@coords[,1],lat=knot_sp@coords[,2],lon_eqd=knot_sp_eqd@coords[,1],lat_eqd=knot_sp_eqd@coords[,2])
			predict_knot = merge(predict_knot,knot_loc_dt,by="knot")

            eta1_dt = as.data.table(fit$Report$eta1_vc) %>% 
            .[,vessel:=levels(factor(bfish_df[,'platform']))] %>%
            .[,vessel:=factor(vessel,levels=c(LETTERS[1:17],"Oscar Elton Sette","Rubber Duck","Steel Toe","Ao Shibi IV"))] %>%
            melt(.,id.vars="vessel") %>%
            setnames(.,c("variable","value"),c("species","eta1")) %>%
			.[,species:=factor(as.character(species),levels=paste0("V",1:7),labels=c("prfi","etca","etco","prsi","przo","hyqu","apru"))] %>%
            .[order(species,vessel)] 

            eta2_dt = as.data.table(fit$Report$eta2_vc) %>% 
            .[,vessel:=levels(factor(bfish_df[,'platform']))] %>%
            .[,vessel:=factor(vessel,levels=c(LETTERS[1:17],"Oscar Elton Sette","Rubber Duck","Steel Toe","Ao Shibi IV"))] %>%
            melt(.,id.vars="vessel") %>%
            setnames(.,c("variable","value"),c("species","eta2")) %>%
			.[,species:=factor(as.character(species),levels=paste0("V",1:7),labels=c("prfi","etca","etco","prsi","przo","hyqu","apru"))] %>%
            .[order(species,vessel)]
            
			eta_dt = merge(eta1_dt,eta2_dt) %>%
                      melt(.,id.vars=c("species","vessel"))

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

			# spatial_residual_dt = data.table(PSU=bfish_df$psu,year=as.numeric(as.character(bfish_df$year)),species=bfish_df$species_cd,residual=residuals_dharma$scaledResiduals)

				psu_table$knot = 1:nrow(psu_table)
			predict_psu = merge(psu_table,predict_knot[,.(knot,species,species_hw,year,density)],by="knot",allow.cartesian=TRUE) %>%
						  .[,area_km2:=0.5^2] %>%
						  merge(.,omega_dt,by=c("knot","species")) %>%
						  merge(.,epsilon_dt,by=c("year","knot","species")) %>%
						  merge(.,component_dt,by=c("year","knot","species"))# %>%
						 # merge(.,spatial_residual_dt,by=c("PSU","year","species"),all.x=TRUE)

			# filter for components that are effectively "turned-off"
			for(i in 1:length(u_species))
			{
				if(max(abs(predict_psu[species == u_species[i]]$omega1)) < 1e-8)
				{
					predict_psu = predict_psu %>%
								  .[species == u_species[i], omega1 := 0]
				}
				if(max(abs(predict_psu[species == u_species[i]]$epsilon1)) < 1e-8)
				{
					predict_psu = predict_psu %>%
								  .[species == u_species[i], epsilon1 := 0]
				}
				if(max(abs(predict_psu[species == u_species[i]]$omega2)) < 1e-8)
				{
					predict_psu = predict_psu %>%
								  .[species == u_species[i], omega2 := 0]
				}
				if(max(abs(predict_psu[species == u_species[i]]$epsilon2)) < 1e-8)
				{
					predict_psu = predict_psu %>%
								  .[species == u_species[i], epsilon2 := 0]
				}
			}

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

                p = copy(eta_dt) %>%
					ggplot() +
					ggtitle("Vessel overdispersion") +
					facet_grid(species~variable,scales="free_y") +
					xlab("Vessel") +
					ylab("Random effect") +
					geom_hline(yintercept=0) +
                    geom_segment(aes(x=vessel,xend=vessel,y=0,yend=value)) +
					geom_point(aes(x=vessel,y=value,fill=vessel),shape=21,size=3) +
					ggthemes::theme_few(base_size=20) + 
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
					viridis::scale_fill_viridis("Vessel",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
			ggsave(filename=paste0("pred_eta_agg.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)

                p = copy(eta_dt) %>%
					ggplot() +
					ggtitle("Vessel overdispersion") +
					facet_grid(~variable,scales="free_y") +
					xlab("Species") +
					ylab("Random effect") +
					geom_hline(yintercept=0) +
					geom_boxplot(aes(x=species,y=value,fill=species)) +

					ggthemes::theme_few(base_size=20) + 
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
					viridis::scale_fill_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
			ggsave(filename=paste0("pred_eta_agg_boxplot.png"), plot = p, device = "png", path = working_dir,
	  			scale = 1.25, width = 9, height = 9, units = c("in"),
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

	  			# p = copy(predict_psu) %>%
				# .[species==u_species[i]] %>%
				# .[,.(Island,PSU,year,knot,lon_eqd,lat_eqd,residual)] %>%
				# .[!is.na(residual)] %>%
				# unique(.) %>%
				# ggplot() + 
				# ggtitle(paste0(unique(predict_psu[species==u_species[i]]$species_hw)," (",toupper(u_species[i]),")")) +
				# facet_wrap(Island~year,nrow=5,scales="free") +
				# xlab("Eastings") +
				# ylab("Northings") +
				# geom_point(aes(x = lon_eqd, y = lat_eqd, fill = residual),shape=21,color="white",size=0.75) +
				# # geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5) +
				# ggthemes::theme_few() + 
				# scale_fill_gradient2(paste0(toupper(residual_type),"\nresidual"),low = "blue",mid = "gray90",high ="red",midpoint=0.5) 
				# ggsave(filename=paste0("pred_residual_psu_",u_species[i],".png"), plot = p, device = "png", path = working_dir,
	  			# scale = 1.25, width = 16, height = 9, units = c("in"),
	  			# dpi = 300, limitsize = TRUE)

	  			# p = copy(predict_psu) %>%
				# .[species==u_species[i]] %>%
				# .[,.(PSU,year,knot,lon_eqd,lat_eqd,residual)] %>%
				# .[!is.na(residual)] %>%
				# .[,.(residual=mean(residual,na.rm=TRUE)),by=.(year,knot)] %>%
				# merge(.,unique(predict_knot[,.(knot,lon_eqd,lat_eqd)]),by="knot") %>%
				# unique(.) %>%
				# ggplot() + 
				# ggtitle(paste0(unique(predict_psu[species==u_species[i]]$species_hw)," (",toupper(u_species[i]),")")) +
				# facet_wrap(~year) +
				# xlab("Eastings") +
				# ylab("Northings") +
				# geom_point(aes(x = lon_eqd, y = lat_eqd, fill = residual),color="white",size=2.25,shape=21) +
				# # geom_sf(data=hi_coast_eqd_sf,fill=NA,alpha=0.5) +
				# ggthemes::theme_few(base_size=20) + 
				# scale_fill_gradient2(paste0(toupper(residual_type),"\nresidual"),low = "blue",mid = "gray90",high ="red",midpoint=0.5) 
				# ggsave(filename=paste0("pred_residual_knot_",u_species[i],".png"), plot = p, device = "png", path = working_dir,
	  			# scale = 1.25, width = 16, height = 9, units = c("in"),
	  			# dpi = 300, limitsize = TRUE)
			}

	# plot index
		if(bias_correct)
		{
			mv_dt = as.data.table(index$Table)
			mv_dt$Estimate = Sdreport$unbiased$value[names(Sdreport$unbiased$value)=="Index_ctl"]
			mv_dt = mv_dt%>%
				.[,Category:=factor(Category,levels=paste0("Category_",1:7),labels=target_species)] %>%
				.[Stratum=="Stratum_1"] %>%
				.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
				.[,Model:="mv"] %>%
				setnames(.,"Std. Error for Estimate","SD") %>%
				.[,.(Model,Category,Time,Estimate,SD)] %>%
				.[,CV:=SD/Estimate] %>%
				.[,SD:=NULL] %>%
				.[,Estimate:=Estimate*2.20462262185] %>%
				.[,Estimate:=Estimate/1000000] %>%
				.[,l95:=exp(log(Estimate)-2*sqrt(log(CV^2+1)))] %>%
				.[,u95:=exp(log(Estimate)+2*sqrt(log(CV^2+1)))]

			mv_agg_dt = as.data.table(index_agg$Table)
			mv_agg_dt$Estimate = Sdreport_agg$unbiased$value[names(Sdreport_agg$unbiased$value)=="Index_ctl"]
			mv_agg_dt = as.data.table(index_agg$Table) %>%
				.[Category=="Category_7"] %>%
				.[,Category:="agg"] %>%
				.[Stratum=="Stratum_1"] %>%
				.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
				.[,Model:="mv"] %>%
				setnames(.,"Std. Error for Estimate","SD") %>%
				.[,.(Model,Category,Time,SD)] %>%
				merge(., mv_dt[,.(Estimate=sum(Estimate)),by=Time],by="Time") %>%
				.[,.(Model,Category,Time,Estimate,SD)] %>%
				.[,Estimate:=Estimate/2.20462262185] %>%
				.[,Estimate:=Estimate*1000000] %>%
				.[,CV:=SD/Estimate] %>%
				.[,SD:=NULL] %>%
				.[,Estimate:=Estimate*2.20462262185] %>%
				.[,Estimate:=Estimate/1000000] %>%
				.[,l95:=exp(log(Estimate)-2*sqrt(log(CV^2+1)))] %>%
				.[,u95:=exp(log(Estimate)+2*sqrt(log(CV^2+1)))]
			
		} else {
			mv_dt = as.data.table(index$Table) %>%
				.[,Category:=factor(Category,levels=paste0("Category_",1:7),labels=target_species)] %>%
				.[Stratum=="Stratum_1"] %>%
				.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
				.[,Model:="mv"] %>%
				setnames(.,"Std. Error for Estimate","SD") %>%
				.[,.(Model,Category,Time,Estimate,SD)] %>%
				.[,CV:=SD/Estimate] %>%
				.[,SD:=NULL] %>%
				.[,Estimate:=Estimate*2.20462262185] %>%
				.[,Estimate:=Estimate/1000000] %>%
				.[,l95:=exp(log(Estimate)-2*sqrt(log(CV^2+1)))] %>%
				.[,u95:=exp(log(Estimate)+2*sqrt(log(CV^2+1)))]

			mv_agg_dt = as.data.table(index_agg$Table) %>%
				.[Category=="Category_7"] %>%
				.[,Category:="agg"] %>%
				.[Stratum=="Stratum_1"] %>%
				.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
				.[,Model:="mv"] %>%
				setnames(.,"Std. Error for Estimate","SD") %>%
				.[,.(Model,Category,Time,Estimate,SD)] %>%
				.[,CV:=SD/Estimate] %>%
				.[,SD:=NULL] %>%
				.[,Estimate:=Estimate*2.20462262185] %>%
				.[,Estimate:=Estimate/1000000] %>%
				.[,l95:=exp(log(Estimate)-2*sqrt(log(CV^2+1)))] %>%
				.[,u95:=exp(log(Estimate)+2*sqrt(log(CV^2+1)))]
			}

		index_dt = rbind(mv_dt,mv_agg_dt)
		fwrite(index_dt,file=paste0(working_dir,"index_dt.csv"))

		p = copy(index_dt) %>%
		.[Category!="agg"] %>%
		ggplot() +
		ylim(0,NA) +
		ylab("Predicted biomass (millions lbs)") +
		xlab("Year") +
		geom_hline(yintercept=0) +
		geom_area(aes(x=Time,y=Estimate,group=Category,fill=Category)) +
		geom_path(data=index_dt[Category=="agg"],aes(x=Time,y=Estimate),color="black",linewidth=1.25) +
		geom_point(data=index_dt[Category=="agg"],aes(x=Time,y=Estimate),color="black",size=4) +
		geom_segment(data=index_dt[Category=="agg"],aes(x=Time,xend=Time,y=l95,yend=u95),color="black",linewidth=1.25) +
		# geom_ribbon(data=index_dt[Category=="agg"],aes(x=Time,ymin=l95,ymax=u95),color="black",linetype="dashed",fill=NA,size=1.25) +
		viridis::scale_color_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
		viridis::scale_fill_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
		theme_few(base_size=20)
		ggsave(filename=paste0("index_by_species_stack.png"), plot = p, device = "png", path = working_dir,
					scale = 1, width = 9, height = 9, units = c("in"),
					dpi = 300, limitsize = TRUE)


		p = copy(index_dt) %>%
		.[Category!="agg"] %>%
			ggplot() +
			ylim(0,NA) +
			ylab("Predicted biomass (millions lbs)") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Category),alpha=0.25) +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Category),linewidth=1.5) +
			viridis::scale_color_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
			ggsave(filename=paste0("index_by_species.png"), plot = p, device = "png", path = working_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)

		nominal_dt = data.table(Category=factor(bfish_df[,'species_cd'],levels=c(target_species)),
								Time=fit$data_list$t_i+2016,
								b_i=bfish_df$weight_kg) %>%
					 .[,.(Estimate=mean(b_i)),by=.(Category,Time)] %>%
					 .[,Estimate:=Estimate/mean(Estimate),by=Category]

		p = copy(index_dt) %>%
		.[Category!="agg"] %>%
			.[,Estimate:=Estimate/mean(Estimate),by=Category] %>%
					.[,l95:=exp(log(Estimate)-2*sqrt(log(CV^2+1)))] %>%
					.[,u95:=exp(log(Estimate)+2*sqrt(log(CV^2+1)))] %>%	  
			ggplot() +
			ylim(0,NA) +
			ylab("Relative abundance") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_hline(yintercept=1,linetype="dashed") +
			geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Category),alpha=0.25) +
			geom_path(data=nominal_dt,aes(x=Time,y=Estimate),linewidth=1.5,color="black",linetype="dotted") +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Category),linewidth=1.5) +
			viridis::scale_color_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
			ggsave(filename=paste0("index_by_species_relative.png"), plot = p, device = "png", path = working_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)
	
		# influence plots
	# spatial effects are not standardized out, they are assumed to impact abundance
        qeffect_1 = fit$data_list$Q1_ik 
		for(i in 1:nrow(qeffect_1))
		{
			qeffect_1[i,] = qeffect_1[i,] * fit$ParHat$lambda1_k
		}
        dimnames(qeffect_1) = dimnames(fit$data_list$Q1_ik)
        qeffectnames_1 = names(fit$effects$catchability_data_full)[-c(1:2,ncol(fit$effects$catchability_data_full))]

        format_qeffect_1 = matrix(NA,nrow=nrow(qeffect_1),ncol=length(qeffectnames_1))
        colnames(format_qeffect_1) = qeffectnames_1
        for(i in 1:length(qeffectnames_1))
        {
            tmp_qeffect = qeffect_1[,grep(qeffectnames_1[i],colnames(qeffect_1),fixed=TRUE)]
            if(is.null(dim(tmp_qeffect)))
            {
                tmp_qeffect = matrix(tmp_qeffect,nrow = length(tmp_qeffect),ncol=1)
            }
            format_qeffect_1[,i] = rowSums(tmp_qeffect)
            rm(list="tmp_qeffect")
        }

        qeffect_1_dt = as.data.table(cbind(bfish_df[,c("year","species_cd")],format_qeffect_1)) %>%
                        .[,component:="1st"]

        qeffect_2 = fit$data_list$Q2_ik 
       for(i in 1:nrow(qeffect_2))
		{
			qeffect_2[i,] = qeffect_2[i,] * fit$ParHat$lambda2_k
		}
	    dimnames(qeffect_2) = dimnames(fit$data_list$Q2_ik)
        qeffectnames_2 = names(fit$effects$catchability_data_full)[-c(1:2,ncol(fit$effects$catchability_data_full))]

        format_qeffect_2 = matrix(NA,nrow=nrow(qeffect_2),ncol=length(qeffectnames_2))
        colnames(format_qeffect_2) = qeffectnames_2
        for(i in 1:length(qeffectnames_2))
        {
            tmp_qeffect = qeffect_2[,grep(qeffectnames_2[i],colnames(qeffect_2),fixed=TRUE)]
            if(is.null(dim(tmp_qeffect)))
            {
                tmp_qeffect = matrix(tmp_qeffect,nrow = length(tmp_qeffect),ncol=1)
            }
            format_qeffect_2[,i] = rowSums(tmp_qeffect)
            rm(list="tmp_qeffect")
        }

        qeffect_2_dt = as.data.table(cbind(bfish_df[,c("year","species_cd")],format_qeffect_2)) %>%
                        .[,component:="2nd"]

        qeffect_dt = rbind(qeffect_1_dt,qeffect_2_dt) %>%
                     melt(.,id.vars=c("year","species_cd","component")) %>%
                     setnames(.,"species_cd","species")
        
        qinflu_dt = copy(qeffect_dt) %>%
                    .[,rho:=mean(value),by=.(component,species,variable)] %>%
				   .[,.(delta_sum=sum(value-rho),.N),by=.(component,species,variable,year)] %>%
				   .[,delta:=delta_sum/N] %>%
				   .[,influ:=exp(delta)] %>%
				   .[,avg_influ:=exp(mean(delta))-1,by=.(component,species,variable)]

        eta1_influ_dt = as.data.table(bfish_df) %>%
        .[,.(year,platform,species_cd)] %>%
        setnames(.,c("platform","species_cd"),c("vessel","species")) %>%
        merge(.,eta1_dt[,.(species,vessel,eta1)],by=c("vessel","species"),all.x=TRUE) %>%
        .[,component:="1st"] %>%
        .[,vessel:=NULL] %>%
        setnames(.,"eta1","eta") %>%
     	  melt(.,id.vars=c("component","species","year")) %>%
          .[,rho:=mean(value),by=.(component,species,variable)] %>%
				   .[,.(delta_sum=sum(value-rho),.N),by=.(component,species,variable,year)] %>%
				   .[,delta:=delta_sum/N] %>%
				   .[,influ:=exp(delta)] %>%
				   .[,avg_influ:=exp(mean(delta))-1,by=.(component,species,variable)]

        eta2_influ_dt = as.data.table(bfish_df) %>%
        .[,.(year,platform,species_cd)] %>%
        setnames(.,c("platform","species_cd"),c("vessel","species")) %>%
        merge(.,eta2_dt[,.(species,vessel,eta2)],by=c("vessel","species"),all.x=TRUE) %>%
        .[,component:="2nd"] %>%
        .[,vessel:=NULL] %>%
        setnames(.,"eta2","eta") %>%
     	  melt(.,id.vars=c("component","species","year")) %>%
          .[,rho:=mean(value),by=.(component,species,variable)] %>%
				   .[,.(delta_sum=sum(value-rho),.N),by=.(component,species,variable,year)] %>%
				   .[,delta:=delta_sum/N] %>%
				   .[,influ:=exp(delta)] %>%
				   .[,avg_influ:=exp(mean(delta))-1,by=.(component,species,variable)]


		influ1_dt = data.table(component="1st",
								  species=factor(bfish_df[,'species_cd'],levels=c(target_species)),
								  year=fit$data_list$t_i+2016,
								  knot=spatial_list$knot_i) %>%
								  merge(.,omega1_dt,by=c("species","knot")) %>%
								  merge(.,epsilon1_dt,by=c("species","knot","year")) %>%
								  setnames(.,c("omega1","epsilon1"),c("omega","epsilon")) %>%
								  melt(.,id.vars=c("component","species","knot","year"))

		influ2_dt = data.table(component="2nd",
								  species=factor(bfish_df[,'species_cd'],levels=c(target_species)),
								  year=fit$data_list$t_i+2016,
								  knot=spatial_list$knot_i) %>%
								  merge(.,omega2_dt,by=c("species","knot")) %>%
								  merge(.,epsilon2_dt,by=c("species","knot","year")) %>%
								  setnames(.,c("omega2","epsilon2"),c("omega","epsilon")) %>%
								  melt(.,id.vars=c("component","species","knot","year"))
		influ_dt = rbind(influ1_dt,influ2_dt) %>%
				   .[,rho:=mean(value),by=.(component,species,variable)] %>%
				   .[,.(delta_sum=sum(value-rho),.N),by=.(component,species,variable,year)] %>%
				   .[,delta:=delta_sum/N] %>%
				   .[,influ:=exp(delta)] %>%
				   .[,avg_influ:=exp(mean(delta))-1,by=.(component,species,variable)] %>%
                    rbind(.,qinflu_dt) %>%
                	rbind(.,eta1_influ_dt) %>%
					rbind(.,eta2_influ_dt)
		fwrite(influ_dt,file=paste0(working_dir,"influ_dt.csv"))

		p = copy(influ_dt) %>% 
            .[!(variable %in% c("omega","epsilon"))] %>%
			ggplot() +
			ylab("Influence") +
			xlab("Year") +
			facet_grid(component~species) +
			geom_hline(yintercept=1,linetype="dashed") +
			geom_line(aes(x=year,y=influ,group=variable,color=variable),linewidth=1.5) +
			viridis::scale_color_viridis("Variable",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
			ggsave(filename=paste0("influ_by_species.png"), plot = p, device = "png", path = working_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)
    
    # effects plots
        factor_effects_1 = qeffectnames_1[which(!(qeffectnames_1 %in% continuous_q_variables))]
        plot_effect1_dt.list = as.list(rep(NA,length(factor_effects_1)))
        for(i in 1:length(factor_effects_1))
        {
            tmp_cols = factor_effects_1[i]
            tmp_pe_dt = qeffect_1_dt[,..tmp_cols]
            colnames(tmp_pe_dt) = "value"
            plot_effect1_dt.list[[i]] = tmp_pe_dt %>%
                        cbind(bfish_df[,c("year","species_cd",factor_effects_1[i])],.) %>%
                        as.data.table(.) %>%
                        setnames(.,c("species_cd",factor_effects_1[i]),c("species","variable")) %>%
                        .[,component:="1st"] %>%
                        .[,covariate:=factor_effects_1[i]] %>%
                        .[,.(component,species,covariate,variable,value)] %>%
                        unique(.) %>%
						.[,species:=factor(species,levels=c(target_species))] %>%
                        .[,value_link:=exp(value)]

        }    

        plot_effect1_dt = rbindlist(plot_effect1_dt.list)
        
        factor_effects_2 = qeffectnames_2[which(!(qeffectnames_2 %in% continuous_q_variables))]
        plot_effect2_dt.list = as.list(rep(NA,length(factor_effects_2)))
        for(i in 1:length(factor_effects_2))
        {
            tmp_cols = factor_effects_2[i]
            tmp_pe_dt = qeffect_2_dt[,..tmp_cols]
            colnames(tmp_pe_dt) = "value"
            plot_effect2_dt.list[[i]] = tmp_pe_dt %>%
                        cbind(bfish_df[,c("year","species_cd",factor_effects_2[i])],.) %>%
                        as.data.table(.) %>%
                        setnames(.,c("species_cd",factor_effects_2[i]),c("species","variable")) %>%
                        .[,component:="2nd"] %>%
                        .[,covariate:=factor_effects_2[i]] %>%
                        .[,.(component,species,covariate,variable,value)] %>%
                        unique(.) %>%
						.[,species:=factor(species,levels=c(target_species))] %>%
                        .[,value_link:=exp(value)]
        }    

        plot_effect2_dt = rbindlist(plot_effect2_dt.list)
        plot_effect_dt = rbind(plot_effect1_dt,plot_effect2_dt)

	    p = copy(plot_effect_dt) %>%
            .[component=="1st"] %>%
			ggplot() +
			ylab("Effect") +
			xlab("Variable") +
			facet_wrap(covariate~species,scales="free_x",nrow=uniqueN(plot_effect_dt$covariate)) +
			geom_hline(yintercept=1,linetype="dashed") +
            geom_segment(aes(x=variable,xend=variable,y=1,yend=value_link)) +
			geom_point(aes(x=variable,y=value_link,fill=value_link),shape=21,size=5) +
            theme_few(base_size=20) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            viridis::scale_fill_viridis("Effect",begin = 0.1,end = 0.8,direction = 1,option = "H")
			ggsave(filename=paste0("effect_by_species_1st.png"), plot = p, device = "png", path = working_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)
	    p = copy(plot_effect_dt) %>%
            .[component=="2nd"] %>%
			ggplot() +
			ylab("Effect") +
			xlab("Variable") +
			facet_wrap(covariate~species,scales="free_x",nrow=uniqueN(plot_effect_dt$covariate)) +
			geom_hline(yintercept=1,linetype="dashed") +
            geom_segment(aes(x=variable,xend=variable,y=1,yend=value_link)) +
			geom_point(aes(x=variable,y=value_link,fill=value_link),shape=21,size=5) +
            theme_few(base_size=20) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            viridis::scale_fill_viridis("Effect",begin = 0.1,end = 0.8,direction = 1,option = "H")
			ggsave(filename=paste0("effect_by_species_2nd.png"), plot = p, device = "png", path = working_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)


        if(length(continuous_q_variables)>0)
        {
            factor_effects_1 = qeffectnames_1[which((qeffectnames_1 %in% continuous_q_variables))]
            plot_effect1_dt.list = as.list(rep(NA,length(factor_effects_1)))
            for(i in 1:length(factor_effects_1))
            {
                tmp_cols = factor_effects_1[i]
                tmp_pe_dt = qeffect_1_dt[,..tmp_cols]
                colnames(tmp_pe_dt) = "value"
                plot_effect1_dt.list[[i]] = tmp_pe_dt %>%
                            cbind(bfish_df[,c("year","species_cd",factor_effects_1[i])],.) %>%
                            as.data.table(.) %>%
                            setnames(.,c("species_cd",factor_effects_1[i]),c("species","variable")) %>%
                            .[,component:="1st"] %>%
                            .[,covariate:=factor_effects_1[i]] %>%
                            .[,.(component,species,covariate,variable,value)] %>%
                            unique(.) %>%
                            .[,species:=factor(species,levels=c(target_species))] %>%
                            .[,value_link:=exp(value)]

            }    

            plot_effect1_dt = rbindlist(plot_effect1_dt.list)
            
            factor_effects_2 = qeffectnames_2[which((qeffectnames_2 %in% continuous_q_variables))]
            plot_effect2_dt.list = as.list(rep(NA,length(factor_effects_2)))
            for(i in 1:length(factor_effects_2))
            {
                tmp_cols = factor_effects_2[i]
                tmp_pe_dt = qeffect_2_dt[,..tmp_cols]
                colnames(tmp_pe_dt) = "value"
                plot_effect2_dt.list[[i]] = tmp_pe_dt %>%
                            cbind(bfish_df[,c("year","species_cd",factor_effects_2[i])],.) %>%
                            as.data.table(.) %>%
                            setnames(.,c("species_cd",factor_effects_2[i]),c("species","variable")) %>%
                            .[,component:="2nd"] %>%
                            .[,covariate:=factor_effects_2[i]] %>%
                            .[,.(component,species,covariate,variable,value)] %>%
                            unique(.) %>%
                            .[,species:=factor(species,levels=c(target_species))] %>%
                            .[,value_link:=exp(value)]
            }    

            plot_effect2_dt = rbindlist(plot_effect2_dt.list)
            plot_effect_dt = rbind(plot_effect1_dt,plot_effect2_dt) %>%
                             .[order(component,species,covariate,variable)]

            p = copy(plot_effect_dt) %>%
            .[component=="1st"] %>%
			ggplot() +
			ylab("Effect") +
			xlab("Variable") +
			facet_grid(covariate~species,scales="free_x") +
			geom_hline(yintercept=1,linetype="dashed") +
			geom_path(aes(x=variable,y=value_link,color=species)) +
            theme_few(base_size=20) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            viridis::scale_fill_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H")
			ggsave(filename=paste0("continuous_effect_by_species_1st.png"), plot = p, device = "png", path = working_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)
	    p = copy(plot_effect_dt) %>%
            .[component=="2nd"] %>%
			ggplot() +
			ylab("Effect") +
			xlab("Variable") +
			facet_grid(covariate~species,scales="free_x") +
			geom_hline(yintercept=1,linetype="dashed") +
			geom_path(aes(x=variable,y=value_link,color=species)) +
            theme_few(base_size=20) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            viridis::scale_fill_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H")
			ggsave(filename=paste0("continuous_effect_by_species_2nd.png"), plot = p, device = "png", path = working_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)

        }
    
    # plot abundance covariates
			ab_df_plot = ab_df
			ab_df_plot = ab_df_plot[,-which(colnames(ab_df_plot)%in%c("depth_sc","complexity_sc","hardness_sc"))]
			colnames(ab_df_plot)[which(colnames(ab_df_plot)%in%c("depth","complexity","hardness"))] = paste0(colnames(ab_df_plot)[which(colnames(ab_df_plot)%in%c("depth","complexity","hardness"))],"_sc")
          
		    abundance_effect_dt.list = as.list(rep(NA,length(continuous_ab_variables)))
            for(i in 1:length(continuous_ab_variables))
            {
                abundance_effect_dt.list[[i]] = as.list(rep(NA,length(target_species)))
                tmp_cols = continuous_ab_variables[i]
                for(s in 1:length(target_species))
                {
                    abundance_effect_dt.list[[i]][[s]] = as.data.table(ab_df_plot) %>%
					.[,..tmp_cols] %>%
                    .[,value:=fit$data_list$X1_gctp[,s,1,grep(continuous_ab_variables[i],dimnames(fit$data_list$X1_gctp)[[4]])] %*% t(t(fit$ParHat$gamma1_cp[s,grep(continuous_ab_variables[i],dimnames(fit$data_list$X1_gctp)[[4]])]))] %>%
                    .[,species_cd:=target_species[s]] %>%
                    setnames(.,tmp_cols,"variable") %>%
                    .[,value_link:=exp(value)] %>%
                    .[,category:=tmp_cols] %>%
                    .[,component:="1st"]

                }
                abundance_effect_dt.list[[i]] = rbindlist(abundance_effect_dt.list[[i]])
                rm(list="tmp_cols")
            }    

            plot_ae1_dt = rbindlist(abundance_effect_dt.list)
             p = copy(plot_ae1_dt) %>%
            .[component=="1st"] %>%
			ggplot() +
			ylab("Effect") +
			xlab("Variable") +
			facet_wrap(category~species_cd,scales="free",ncol=7) +
			geom_hline(yintercept=1,linetype="dashed") +
			geom_line(aes(x=variable,y=value_link,color=species_cd)) +
            theme_few(base_size=20) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            viridis::scale_color_viridis("Species",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
			ggsave(filename=paste0("continuous_abundance_effect_by_species_1st.png"), plot = p, device = "png", path = working_dir,
						scale = 1, width = 16, height = 9, units = c("in"),
						dpi = 300, limitsize = TRUE)
#_____________________________________________________________________________________________________________________________
# 6) cross-validation

	A = proc.time()
	if(xval == "xval")
	{
		k_folds = length(list.files(load_xval_path))/4
		valid_vec = rep(NA,k_folds)
		xval_dt.list = as.list(rep(NA,k_folds))

		for(k in 1:k_folds)
		{
			load(file=paste0(load_xval_path,k,"_xval_extrapolation_list.RData"))
        	load(file=paste0(load_xval_path,k,"_xval_spatial_list.RData"))
        	load(file=paste0(load_xval_path,k,"_xval_train_df.RData"))
        	load(file=paste0(load_xval_path,k,"_xval_test_df.RData"))

			if(target_species != "mv")
			{
				tmp_weight_kg = xval_train_df[,target_species]
				xval_train_df$weight_kg = tmp_weight_kg
				xval_train_df$species_cd = target_species
				rm(list="tmp_weight_kg")
			}
			
			xval_working_dir = paste0(working_dir,"xval/",k,"/")
			dir.create(xval_working_dir,recursive = TRUE)

			xval_fit = fit_model( settings=settings,
 				   				Lon_i=xval_train_df$lon,
    			   				Lat_i=xval_train_df$lat,
          						t_i=as.integer(xval_train_df$year),
          						b_i=xval_train_df$weight_kg,
          						a_i=rep(pi * (0.02760333457^2),nrow(xval_train_df)), # assumed area swept from the MOUSS camera converted to km2; Ault et al 2018
	    	  					c_i = as.numeric(factor(xval_train_df[,'species_cd'],levels=c(target_species)))-1,
	    	  					category_names=c(target_species),
	    	  					covariate_data = NULL,
	    	  					catchability_data = NULL,
          						working_dir = xval_working_dir,
          						newtonsteps = 1,
          						# extrapolation list args
          						extrapolation_list = xval_extrapolation_list,
          						# spatial list args
    	  						spatial_list = xval_spatial_list,
    	  						test_fit=TRUE)
			
			if("xval_fit" %in% ls())
			{
				valid_vec[k] = k

				if(target_species != "mv")
				{
					tmp_weight_kg = xval_test_df[,target_species]
					xval_test_df$weight_kg = tmp_weight_kg
					xval_test_df$species_cd = target_species
					rm(list="tmp_weight_kg")
				}

				xval_predict_working_dir = paste0(working_dir,"xval/",k,"/predict/")
				dir.create(xval_predict_working_dir,recursive = TRUE)

				xval_fit$input_args$spatial_args_input = list(
														n_x=xval_spatial_list$n_x,
														Lon_i=xval_test_df$lon,
														Lat_i=xval_test_df$lat,
														Extrapolation_List=xval_extrapolation_list,
														Method = "Barrier",
														anisotropic_mesh = xval_spatial_list$MeshList$anisotropic_mesh,
														grid_size_km = 0.5,
														grid_size_LL = 0.5/110,
														fine_scale = xval_spatial_list$fine_scale,
														Save_Results = FALSE,
														LON_intensity=xval_spatial_list$latlon_x[,"Lon"],
														LAT_intensity=xval_spatial_list$latlon_x[,"Lat"],
														map_data=hi_coast
    	  												)

				xval_predict = predict(xval_fit,
                  what="D_i",
                  Lat_i=xval_test_df$lat,
                  Lon_i=xval_test_df$lon,
                  t_i=as.integer(xval_test_df$year),
                  a_i=rep(pi * (0.02760333457^2),nrow(xval_test_df)),
                  c_iz = as.numeric(factor(xval_test_df[,'species_cd'],levels=c(target_species)))-1,
                  do_checks = TRUE,
                  working_dir = xval_predict_working_dir )

				  xval_test_df$pred_weight_kg = xval_predict
				  xval_test_df$partition = k

				  xval_dt.list[[k]] = as.data.table(xval_test_df) %>%
				  					  .[,.(partition,model_sampling_unit,year,lon,lat,species_cd,weight_kg,pred_weight_kg)]


				  rm(list=c("xval_predict_working_dir","xval_predict"))
				  gc();

			}

			# clean-up
        	rm(list=c("xval_fit","xval_working_dir","xval_extrapolation_list","xval_spatial_list","xval_train_df","xval_test_df"))
			gc();
		}

		xval_dt = rbindlist(xval_dt.list[na.omit(valid_vec)])
		xval_rmse_dt = copy(xval_dt) %>%
					   .[,.(rmse=sqrt(mean((weight_kg-pred_weight_kg)^2)),nrmse=sqrt(mean((weight_kg-pred_weight_kg)^2))/sd(weight_kg),.N),by=.(partition)]
		fwrite(xval_dt,file=paste0(working_dir,"xval_dt.csv"))
		fwrite(xval_rmse_dt,file=paste0(working_dir,"xval_rmse_dt.csv"))

	}
	B = proc.time()
	B-A
