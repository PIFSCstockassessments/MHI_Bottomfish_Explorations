

# Nicholas Ducharme-Barth
# 07/20/2022
# Define VAST spatial_mesh

# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(VAST)
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
# 2) set-up extrapolation_list
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
		extrapolation_list = make_extrapolation_info(Region="user",projargs=slot(crs_eqd,"projargs"),input_grid=input_grid)
	# add strata definitions for each Island group, STRATA and STRATA_2020
		a_el_orig = extrapolation_list$a_el[,1]
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

		extrapolation_list$a_el = a_el_tmp
		units(extrapolation_list$a_el) = "km^2"
		save(extrapolation_list,file=paste0(proj.dir,"VAST/extrapolation_list/extrapolation_list.bfish_default.RData"))

#_____________________________________________________________________________________________________________________________
# 3) set-up spatial_list: define mesh boundary
# iterate over options: fine_scale and knot density
	# define input mesh
	# setting the boundary helps restrict the placement of knots to where the data is
	# the closer to zero convex is, the more intricate the boundary can become and the larger resolution is needed
	# also a more intricate boundary will result in longer mesh generation time and a greater number of knots
	# this takes a few seconds
	A = proc.time()
		mesh_boundary.a = INLA::inla.nonconvex.hull(cbind(psu_table$lon_eqd, psu_table$lat_eqd), convex = -0.005,resolution=205) # original boundary definition
	B = proc.time()
	print(paste0("mesh boundary took ",round((B-A)[3],digits=2)," seconds to generate"))
	A = proc.time()
		mesh_boundary.b = INLA::inla.nonconvex.hull(cbind(psu_table$lon_eqd, psu_table$lat_eqd), convex = -0.00375,resolution=272)
	B = proc.time()
	print(paste0("mesh boundary took ",round((B-A)[3],digits=2)," seconds to generate"))
	A = proc.time()
		mesh_boundary.c = INLA::inla.nonconvex.hull(cbind(psu_table$lon_eqd, psu_table$lat_eqd), convex = -0.0025,resolution=405)
	B = proc.time()
	print(paste0("mesh boundary took ",round((B-A)[3],digits=2)," seconds to generate"))
	A = proc.time()
		mesh_boundary.d = INLA::inla.nonconvex.hull(cbind(psu_table$lon_eqd, psu_table$lat_eqd), convex = -0.0015,resolution=672)
	B = proc.time()
	print(paste0("mesh boundary took ",round((B-A)[3],digits=2)," seconds to generate"))

		save(mesh_boundary.a,file=paste0(proj.dir,"VAST/spatial_list/mesh_boundary.a.bfish_default.RData"))
		save(mesh_boundary.b,file=paste0(proj.dir,"VAST/spatial_list/mesh_boundary.b.bfish_default.RData"))
		save(mesh_boundary.c,file=paste0(proj.dir,"VAST/spatial_list/mesh_boundary.c.bfish_default.RData"))
		save(mesh_boundary.d,file=paste0(proj.dir,"VAST/spatial_list/mesh_boundary.d.bfish_default.RData"))


	# make plot of boundary
		boundary_coords = function(boundary_obj)
		{
			idx_boundary = boundary_obj$idx
			# add NA at each non-consecutive index
			non_consec_idx = c(1,which(diff(idx_boundary[,2])<0)+1)
			new_idx_boundary = matrix(c(NA,1),nrow=1,ncol=2)
			for(j in 1:(length(non_consec_idx)-1))
			{
				new_idx_boundary = rbind(new_idx_boundary,idx_boundary[non_consec_idx[j]:non_consec_idx[j+1],],c(idx_boundary[non_consec_idx[j+1],1],NA))
			}

			# remove duplicate indices for plotting
			    new_idx_boundary = unique(new_idx_boundary)
			    x_boundary = boundary_obj$loc[t(new_idx_boundary), 1]
			    y_boundary = boundary_obj$loc[t(new_idx_boundary), 2]

			return(data.frame(x=x_boundary,y=y_boundary)) 
		} 
		
		png(filename = paste0(proj.dir,"/Plot/BFISH_CPUE/VAST/spatial/mesh_boundary_multipanel.png"),width = 16, height = 9, units = "in",res=300)
		 	par(mar=c(5,5,1,1),mfrow=c(2,2))
			sp::plot(hi_coast_eqd,axes=TRUE,col="gray90",ylim=c(-300,300),xlab="Eastings (km)",ylab="Northings (km)",cex=1.5,cex.axis=1.5,cex.lab=1.5,las=1)
			points(psu_table$lon_eqd, psu_table$lat_eqd,pch=3,cex=0.1,col="gray50")
			lines(boundary_coords(mesh_boundary.a),col="blue",lwd=2)

			sp::plot(hi_coast_eqd,axes=TRUE,col="gray90",ylim=c(-300,300),xlab="Eastings (km)",ylab="Northings (km)",cex=1.5,cex.axis=1.5,cex.lab=1.5,las=1)
			points(psu_table$lon_eqd, psu_table$lat_eqd,pch=3,cex=0.1,col="gray50")
			lines(boundary_coords(mesh_boundary.b),col="orange",lwd=2)

			sp::plot(hi_coast_eqd,axes=TRUE,col="gray90",ylim=c(-300,300),xlab="Eastings (km)",ylab="Northings (km)",cex=1.5,cex.axis=1.5,cex.lab=1.5,las=1)
			points(psu_table$lon_eqd, psu_table$lat_eqd,pch=3,cex=0.1,col="gray50")
			lines(boundary_coords(mesh_boundary.c),col="green",lwd=2)

			sp::plot(hi_coast_eqd,axes=TRUE,col="gray90",ylim=c(-300,300),xlab="Eastings (km)",ylab="Northings (km)",cex=1.5,cex.axis=1.5,cex.lab=1.5,las=1)
			points(psu_table$lon_eqd, psu_table$lat_eqd,pch=3,cex=0.1,col="gray50")
			lines(boundary_coords(mesh_boundary.d),col="hotpink",lwd=2)
			legend("bottomleft",legend=c("PSU domain","mesh boundary - a","mesh boundary - b","mesh boundary - c","mesh boundary - d","land"),lwd=c(NA,3,3,3,3,1),pch=c(16,NA,NA,NA,NA,NA),col=c("gray50","blue","orange","green","hotpink",NA),fill=c(NA,NA,NA,NA,NA,"gray90"),border=c(NA,NA,NA,NA,NA,"black"),bty="n",cex=1.5)
		dev.off()

		png(filename = paste0(proj.dir,"/Plot/BFISH_CPUE/VAST/spatial/mesh_boundary.png"),width = 16, height = 9, units = "in",res=300)
		 	par(mar=c(5,5,1,1))
			sp::plot(hi_coast_eqd,axes=TRUE,col="gray90",ylim=c(-300,300),xlab="Eastings (km)",ylab="Northings (km)",cex=1.5,cex.axis=1.5,cex.lab=1.5,las=1)
			points(psu_table$lon_eqd, psu_table$lat_eqd,pch=3,cex=0.1,col="gray50")
			lines(boundary_coords(mesh_boundary.a),col="blue",lwd=2)
			lines(boundary_coords(mesh_boundary.b),col="orange",lwd=2)
			lines(boundary_coords(mesh_boundary.c),col="green",lwd=2)
			lines(boundary_coords(mesh_boundary.d),col="hotpink",lwd=2)
			legend("bottomleft",legend=c("PSU domain","mesh boundary - a","mesh boundary - b","mesh boundary - c","mesh boundary - d","land"),lwd=c(NA,3,3,3,3,1),pch=c(16,NA,NA,NA,NA,NA),col=c("gray50","blue","orange","green","hotpink",NA),fill=c(NA,NA,NA,NA,NA,"gray90"),border=c(NA,NA,NA,NA,NA,"black"),bty="n",cex=1.5)
		dev.off()

#_____________________________________________________________________________________________________________________________
# 4) set-up spatial_list
# iterate over options: fine_scale and knot density
	spatial_options_df = expand.grid(mesh_boundary=c("a","b","c"),km_cutoff=c(6,4,2),fine_scale=c(TRUE,FALSE))

	load(file=paste0(proj.dir,"VAST/extrapolation_list/extrapolation_list.bfish_default.RData"))
	load(file=paste0(proj.dir,"VAST/spatial_list/mesh_boundary.a.bfish_default.RData"))
	load(file=paste0(proj.dir,"VAST/spatial_list/mesh_boundary.b.bfish_default.RData"))
	load(file=paste0(proj.dir,"VAST/spatial_list/mesh_boundary.c.bfish_default.RData"))

	for(i in 1:nrow(spatial_options_df))
	{

		# in inla.mesh.2d the max.n.strict argument restricts the number of knots more than setting max.n. However, the cutoff distance between knots is what really controls the number of knots.
		mesh_inla = INLA::inla.mesh.2d(loc = psu_sp_eqd,boundary = get(paste0("mesh_boundary.",spatial_options_df$mesh_boundary[i])),max.n.strict=c(500,16),cutoff=spatial_options_df$km_cutoff[i])
		save(mesh_inla,file=paste0(proj.dir,"VAST/spatial_list/mesh_inla.",spatial_options_df$mesh_boundary[i],".",spatial_options_df$km_cutoff[i],".bfish_default.RData"))

		# define knot locations for spatial_list
			mesh_coords_eqd = sp::SpatialPoints(mesh_inla$loc[,1:2])
			sp::proj4string(mesh_coords_eqd) = crs_eqd
			mesh_coords_ll = sp::spTransform(mesh_coords_eqd, crs_ll)
			intensity_loc = mesh_coords_ll@coords 
		
		spatial_list = make_spatial_info(n_x=nrow(mesh_coords_ll@coords),
    					   Lon_i=bfish_df$lon,
    					   Lat_i=bfish_df$lat,
    					   Extrapolation_List=extrapolation_list,
    					   Method = "Barrier",
    					   anisotropic_mesh = mesh_inla,
    					   grid_size_km = 0.5,
    					   grid_size_LL = 0.5/110,
    					   fine_scale = spatial_options_df$fine_scale[i],
    					   Save_Results = FALSE,
    					   LON_intensity=intensity_loc[,1],
    					   LAT_intensity=intensity_loc[,2],
    					   map_data=hi_coast)
		save(spatial_list,file=paste0(proj.dir,"VAST/spatial_list/spatial_list.",spatial_options_df$mesh_boundary[i],".",spatial_options_df$km_cutoff[i],".",spatial_options_df$fine_scale[i],".bfish_default.RData"))
		
		# add plot to show the mesh structure
		# pull-out needed quantities for plotting
			t.sub=1:nrow(spatial_list$MeshList$anisotropic_mesh$graph$tv)
		    idx = cbind(spatial_list$MeshList$anisotropic_mesh$graph$tv[t.sub,c(1:3,1), drop=FALSE], NA)
		    x = spatial_list$MeshList$anisotropic_mesh$loc[t(idx), 1]
	        y = spatial_list$MeshList$anisotropic_mesh$loc[t(idx), 2]
		    idx_normal = idx[-spatial_list$MeshList$anisotropic_mesh_triangles_over_land,]
	        x_normal = spatial_list$MeshList$anisotropic_mesh$loc[t(idx_normal), 1]
	        y_normal =spatial_list$MeshList$anisotropic_mesh$loc[t(idx_normal), 2]

 			normal_col = "black"
 			barrier_col = "gray70"
 			boundary_col = "hotpink"
 			land_col = "gray90"
	 		
			png(filename = paste0(proj.dir,"/Plot/BFISH_CPUE/VAST/spatial/spatial_mesh.",spatial_options_df$mesh_boundary[i],".",spatial_options_df$km_cutoff[i],".",spatial_options_df$fine_scale[i],".png"),width = 16, height = 9, units = "in",res=300)
		 		par(mar=c(5,5,1,1))
				sp::plot(hi_coast_eqd,axes=TRUE,col=land_col,ylim=c(-300,300),xlab="Eastings (km)",ylab="Northings (km)",cex=1.5,cex.axis=1.5,cex.lab=1.5,las=1)
				lines(x,y,col=barrier_col)
				points(x,y,pch=16,cex=0.5,col=barrier_col)
				lines(x_normal,y_normal,col=normal_col)
				lines(boundary_coords(get(paste0("mesh_boundary.",spatial_options_df$mesh_boundary[i]))),col=boundary_col,lwd=2)
				points(x_normal,y_normal,pch=16,cex=0.5,col=normal_col)
				legend("bottomleft",legend=c("normal knot","barrier knot","normal edge","barrier edge","knot boundary","land"),lwd=c(NA,NA,1,1,3,1),pch=c(16,16,NA,NA,NA,NA),col=c(normal_col,barrier_col,normal_col,barrier_col,boundary_col,NA),fill=c(NA,NA,NA,NA,NA,land_col),border=c(NA,NA,NA,NA,NA,"black"),bty="n",cex=1.5)
			dev.off()
	}
