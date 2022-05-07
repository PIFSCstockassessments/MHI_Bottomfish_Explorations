

# Nicholas Ducharme-Barth
# 04/14/2022
# Set up simple spatiotemporal model example using sdmTMB
# 1) Bring data (including spatial data)
# 2) Set-up barrier mesh (including conversion to equal distant projection)
# 3) Fit single species (e.g. prfi - Opakapaka) model using Tweedie distribution
# 4) Examine diagnostics
# 5) Make predictions on PSU grid & examine outputs
# 6) Calculate index

# This combines topics covered in sdmTMB vignettes (https://pbs-assess.github.io/sdmTMB/index.html)

# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(sdmTMB)
	library(ggplot2)
	# rgdal functions sourced directly
	# sp functions sourced directly
	# DHARMa functions sourced directly
	# sf functions sourced directly
	# ggthemes functions sourced directly
#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"
#_____________________________________________________________________________________________________________________________
# define data_flag
	# data_flag = "" # only loads data up through 2020
	data_flag = "2021_" # includes data through 2021
#_____________________________________________________________________________________________________________________________
# 1) bring in data
	load(file=paste0(proj.dir,"Data/",data_flag,"bfish_combined_long_dt.RData"))
	# subset to single species
	bfish_df = bfish_combined_long_dt %>% .[species_cd == "prfi"] %>% as.data.frame(.)
	# needed to define spatial domain and for predicting on to create index
	psu_table = fread(paste0(proj.dir,"Data/BFISH PSU lookup table.csv")) %>%
				.[,.(PSU,Island,lon_deg,lat_deg,STRATA,STRATA_2020,Depth_MEDIAN_m,med_slp,med_acr,BS_pct_over_136j,pctHB,pctHS)] %>%
				.[,substrate:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][1])] %>%
				.[,slope:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][2])]
	# needed for plotting and defining the barrier feature
	hi_coast = rgdal::readOGR(dsn = paste0(proj.dir,"Data/GIS/Coastline"), layer = "Coastline")
	hi_coast_sf = sf::st_as_sf(hi_coast)
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

	# define mesh
		# setting the boundary helps restrict the placement of knots to where the data is
		# the smaller convex is, the more intricate the boundary can become and the larger resolution is needed
		# also a more intricate boundary will result in longer mesh generation time and a greater number of knots
		# this takes a few seconds
		A = proc.time()
		mesh_boundary = INLA::inla.nonconvex.hull(cbind(psu_table$lon_eqd, psu_table$lat_eqd), convex = -0.005,resolution=205)
		B = proc.time()
		print(paste0("mesh boundary took ",round((B-A)[3],digits=2)," seconds to generate"))
		# in inla.mesh.2d the max.n.strict argument restricts the number of knots more than setting max.n. However, the cutoff distance between knots is what really controls the number of knots.
		# for this example leave cutoff == 7.5 (in km) to keep the number of knots manageable, but for full analysis can decrease to << 4.
		mesh_inla = INLA::inla.mesh.2d(loc = psu_sp_eqd,boundary = mesh_boundary,max.n.strict=c(500,16),cutoff=7.5)
		print(paste0(mesh_inla$n," <- number of knots in spatial mesh"))  
		mesh_sdmTMB = make_mesh(bfish_df, c("lon_eqd", "lat_eqd"), mesh = mesh_inla)
		# range fraction controls the degree to which the correlation is reduced across barriers
		mesh_barrier = add_barrier_mesh(mesh_sdmTMB,hi_coast_eqd_sf,range_fraction = 0.01,plot=FALSE)

		# add plot to show the mesh structure
		# pull-out needed quantities for plotting
			t.sub=1:nrow(mesh_barrier$mesh$graph$tv)
		    idx = cbind(mesh_barrier$mesh$graph$tv[t.sub,c(1:3,1), drop=FALSE], NA)
		    x = mesh_barrier$mesh$loc[t(idx), 1]
	        y = mesh_barrier$mesh$loc[t(idx), 2]
		    idx_normal = idx[-mesh_barrier$barrier_triangles,]
	        x_normal = mesh_barrier$mesh$loc[t(idx_normal), 1]
	        y_normal = mesh_barrier$mesh$loc[t(idx_normal), 2]
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
	 		
	 		par(mar=c(5,5,1,1))
			sp::plot(hi_coast_eqd,axes=TRUE,col=land_col,ylim=c(-300,300),xlab="Eastings (km)",ylab="Northings (km)",cex=1.5,cex.axis=1.5,cex.lab=1.5,las=1)
			lines(x,y,col=barrier_col)
			points(x,y,pch=16,cex=0.5,col=barrier_col)
			lines(x_normal,y_normal,col=normal_col)
			lines(x_boundary,y_boundary,col=boundary_col,lwd=3)
			points(x_normal,y_normal,pch=16,cex=0.5,col=normal_col)
			legend("bottomleft",legend=c("normal knot","barrier knot","normal edge","barrier edge","knot boundary","land"),lwd=c(NA,NA,1,1,3,1),pch=c(16,16,NA,NA,NA,NA),col=c(normal_col,barrier_col,normal_col,barrier_col,boundary_col,NA),fill=c(NA,NA,NA,NA,NA,land_col),border=c(NA,NA,NA,NA,NA,"black"),bty="n",cex=1.5)

#_____________________________________________________________________________________________________________________________
# 3) fit a basic sdmTMB model
# Error structure: Tweedie
# Fixed effects: Year
# Random effects: spatial + spatiotemporal
	# this takes a few seconds
	A = proc.time()
	fit = sdmTMB(
	  data = bfish_df, 
	  formula = weight_kg ~ 0 + as.factor(year),
	  time = "year",
	  mesh = mesh_barrier, 
	  family = tweedie(link = "log"))
	B = proc.time()
	B-A
	print(paste0(round((B-A)[3],digits=2)," seconds to fit spatiotemporal model"))

#_____________________________________________________________________________________________________________________________
# 4) Model diagnostics & plot output
	# need to simulate responses to calculate DHARMa style residuals
	# can use built in functions from sdmTMB package but DHARMa package has some built in tests/plots
	set.seed(123)
	sim_fit = simulate(fit, nsim = 500)
	pred = fit$family$linkinv(predict(fit)$est)
	residuals_dharma =  DHARMa::createDHARMa(simulatedResponse = sim_fit,observedResponse = bfish_df$weight_kg,fittedPredictedResponse = pred)	
	# basic QQ & residual v. predicted plot
		plot(residuals_dharma)
	# residual v. predicted by model covariate
		DHARMa::plotResiduals(residuals_dharma, form = as.factor(bfish_df$year))
		DHARMa::plotResiduals(residuals_dharma, form = bfish_df$lon)
		DHARMa::plotResiduals(residuals_dharma, form = bfish_df$lat)
	# test for uniformity in residuals, overdispersion, outliers
		DHARMa::testResiduals(residuals_dharma)
	# test for zero inflation
		DHARMa::testZeroInflation(residuals_dharma)
	# test for spatial autocorrelation
	# DHARMa::testSpatialAutocorrelation needs unique locations
		bfish_df$spatial_group = as.numeric(as.factor(paste0(bfish_df$lon,"_",bfish_df$lat)))
		spatial_group_dt = data.table(spatial_group=bfish_df$spatial_group,x=bfish_df$lon,y=bfish_df$lat) %>%
					   .[,.(x=mean(x),y=mean(y)),by=spatial_group]
		residuals_spatial_group = DHARMa::recalculateResiduals(residuals_dharma, group = bfish_df$spatial_group)	
		DHARMa::testSpatialAutocorrelation(residuals_spatial_group, x = spatial_group_dt$x, y = spatial_group_dt$y)
	# test for temporal autocorrelation
		bfish_df$temporal_group = factor(bfish_df$year,levels=as.character(sort(unique(bfish_df$year))))
		residuals_temporal_group = DHARMa::recalculateResiduals(residuals_dharma, group = bfish_df$temporal_group)	
		DHARMa::testTemporalAutocorrelation(residuals_temporal_group, time=levels(bfish_df$temporal_group))
		
		bfish_df$temporal_group =factor(bfish_df$year_continuous,levels=as.character(sort(unique(bfish_df$year_continuous))))
		residuals_temporal_group = DHARMa::recalculateResiduals(residuals_dharma, group = bfish_df$temporal_group)	
		DHARMa::testTemporalAutocorrelation(residuals_temporal_group, time=levels(bfish_df$temporal_group))

#_____________________________________________________________________________________________________________________________
# 5) Make predictions across model domain (all PSUs across years)
		# expand psu_table by year
			psu_table.list = lapply(unique(bfish_df$year),function(x,dt){tmp=copy(dt)%>%.[,year:=x];return(tmp)},dt=psu_table)
			psu_table_annual = rbindlist(psu_table.list)

		# predict on annual psu table
			predict_psu = predict(fit, newdata = psu_table_annual)
			# conduct boot-strapping from joint precision matrix to get estimates of uncertainty
			predict_sim = predict(fit, newdata = psu_table_annual, nsim=500)
			predict_psu_sim = cbind(psu_table_annual,predict_sim) %>%
							  melt(.,id.vars=colnames(psu_table_annual))

		# plot estimated response
			predict_psu %>%
			ggplot() + 
			facet_wrap(~year) +
			xlab("Longitude") +
			ylab("Latitude") +
			geom_point(aes(x = lon_eqd, y = lat_eqd, color = exp(est)),size=0.05) +
			geom_sf(data=hi_coast_eqd_sf,fill=land_col,alpha=0.5) +
			ggthemes::theme_few(base_size=20) + 
			viridis::scale_color_viridis("Est. response",begin = 0.1,end = 0.8,direction = 1,option = "H",trans="log10")

		# plot estimated response uncertainty
			predict_psu_sim %>% .[,.(mean=mean(exp(value)),sd=sd(exp(value))),by=.(year,lon_eqd,lat_eqd)] %>% .[,cv:=sd/mean] %>%
			ggplot() + 
			facet_wrap(~year) +
			xlab("Longitude") +
			ylab("Latitude") +
			geom_point(aes(x = lon_eqd, y = lat_eqd, color = cv),size=0.05) +
			geom_sf(data=hi_coast_eqd_sf,fill=land_col,alpha=0.5) +
			ggthemes::theme_few(base_size=20) + 
			viridis::scale_color_viridis("CV (response)",begin = 0.1,end = 0.8,direction = 1,option = "H",trans="sqrt")

		# Estimated fixed effects
			predict_psu %>%
			ggplot() + 
			facet_wrap(~year) +
			xlab("Longitude") +
			ylab("Latitude") +
			geom_point(aes(x = lon_eqd, y = lat_eqd, color = exp(est_non_rf)),size=0.05) +
			geom_sf(data=hi_coast_eqd_sf,fill=land_col,alpha=0.5) +
			ggthemes::theme_few(base_size=20) + 
			viridis::scale_color_viridis("Est. fixed effects",begin = 0.1,end = 0.8,direction = 1,option = "H",trans="log10")

		# Estimated spatial effects
			predict_psu %>% .[,.(lon_eqd,lat_eqd,omega_s)] %>% unique(.) %>%
			ggplot() + 
			xlab("Longitude") +
			ylab("Latitude") +
			geom_point(aes(x = lon_eqd, y = lat_eqd, color = omega_s),size=0.05) +
			geom_sf(data=hi_coast_eqd_sf,fill=land_col,alpha=0.5) +
			ggthemes::theme_few(base_size=20) + 
			scale_color_gradient2("Est. spatial random effects",low="blue",high="red")

		# Estimated spatiotemporal effects
			predict_psu %>%
			ggplot() + 
			facet_wrap(~year) +
			xlab("Longitude") +
			ylab("Latitude") +
			geom_point(aes(x = lon_eqd, y = lat_eqd, color = epsilon_st),size=0.05) +
			geom_sf(data=hi_coast_eqd_sf,fill=land_col,alpha=0.5) +
			ggthemes::theme_few(base_size=20) + 
			scale_color_gradient2("Est. spatiotemporal random effects",low="blue",high="red")

#_____________________________________________________________________________________________________________________________
# 6) Calculate index
		index_dt = predict_psu_sim %>% .[,.(value=mean(exp(value))),by=.(year,variable)] %>%
							.[,.(est=quantile(value,probs=0.5),lwr=quantile(value,probs=0.025),upr=quantile(value,probs=0.975),log_est=mean(log(value)),se=sd(log(value))),by=year]
		index_regional_dt = predict_psu_sim %>% .[,.(value=mean(exp(value))),by=.(year,Island,variable)] %>%
							.[,.(est=quantile(value,probs=0.5),lwr=quantile(value,probs=0.025),upr=quantile(value,probs=0.975),log_est=mean(log(value)),se=sd(log(value))),by=.(year,Island)] %>%
							.[,Island:=factor(Island,levels=c("Niihau","Kauai","Oahu","Maui Nui","Big Island"),labels=c("Ni'ihau","Kaua'i","O'ahu","Maui Nui","Hawai'i"))]


		# plot index with log-normal error
		index_dt %>%
		ggplot() +
		ggtitle("BFISH relative abundance: Opakapaka") +
		ylim(0,NA) +
		geom_hline(yintercept=0) +
		xlab("Year") +
		ylab("Avg. biomass (kg)") +
		geom_ribbon(aes(x=as.numeric(year),ymin=lwr,ymax=upr),alpha=0.5) +
		geom_path(aes(x=as.numeric(year),y=est),size=1.25) +
		ggthemes::theme_few(base_size=20)

		# plot regional indices with log-normal error
		index_regional_dt %>%
		ggplot() +
		facet_wrap(~Island,scales="free_y") +
		ggtitle("BFISH relative abundance: Opakapaka") +
		ylim(0,NA) +
		geom_hline(yintercept=0) +
		xlab("Year") +
		ylab("Avg. biomass (kg)") +
		geom_ribbon(aes(x=as.numeric(year),ymin=lwr,ymax=upr,fill=Island),alpha=0.25) +
		geom_path(aes(x=as.numeric(year),y=est,color=Island),size=1.25) +
		viridis::scale_color_viridis("Island group",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete="TRUE") +
		viridis::scale_fill_viridis("Island group",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete="TRUE") +
		ggthemes::theme_few(base_size=20)

