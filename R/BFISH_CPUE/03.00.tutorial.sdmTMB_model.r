

# Nicholas Ducharme-Barth
# 04/14/2022
# Set up simple spatiotemporal model example using sdmTMB
# 1) Bring data (including spatial data)
# 2) Set-up barrier mesh (including conversion to equal distant projection)
# 3) Fit single species (e.g. prfi - Opakapaka) model using Tweedie distribution
# 4) Examine diagnostics
# 5) Make predictions
# 6)

# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(sdmTMB)
	# rgdal functions sourced directly
	# sp functions sourced directly

#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"

#_____________________________________________________________________________________________________________________________
# 1) bring in data
	load(file=paste0(proj.dir,"Data/bfish_combined_long_dt.RData"))
	# subset to single species
	bfish_df = bfish_combined_long_dt %>% .[species_cd == "prfi"] %>% as.data.frame(.)
	# needed to define spatial domain and for predicting on to create index
	psu_table = fread(paste0(proj.dir,"Data/BFISH PSU lookup table.csv")) %>%
				.[,.(PSU,Island,lon_deg,lat_deg,STRATA,STRATA_2020,Depth_MEDIAN_m,med_slp,med_acr,BS_pct_over_136j,pctHB,pctHS)] %>%
				.[,substrate:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][1])] %>%
				.[,slope:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][2])]
	# needed for plotting and defining the barrier feature
	hi_coast = rgdal::readOGR(dsn = paste0(proj.dir,"Data/GIS/Coastline"), layer = "Coastline")

#_____________________________________________________________________________________________________________________________
# 2) set-up spatial domain
	# convert to equal distant projection
		# get original lat-lon crs
		crs_ll = sp::CRS(sp::proj4string(hi_coast))
		# use two-point equi-distant projection
		crs_eqd = sp::CRS("+proj=tpeqd +lat_1=20.45 +lon_1=-158.9 +lat_2=20.45 +lon_2=-156.2 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
	  	hi_coast_eqd = sp::spTransform(hi_coast, crs_eqd)
	  	hi_coast_eqd_sf = sf::st_as_sf(hi_coast_eqd) 

	# get equi-distant coordinates for samples & psu
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
		print(paste0("mesh boundary took ",(B-A)[3]," seconds to generate"))
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
	        # remove duplicate indices
	        new_idx_boundary = unique(new_idx_boundary)
	        x_boundary = mesh_boundary$loc[t(new_idx_boundary), 1]
	        y_boundary = mesh_boundary$loc[t(new_idx_boundary), 2] 
 		
	 		par(mar=c(5,5,1,1))
			sp::plot(hi_coast_eqd,axes=TRUE,col="gray90",ylim=c(-300,300),xlab="Eastings (km)",ylab="Northings (km)",cex=1.5,cex.axis=1.5,cex.lab=1.5,las=1)
			lines(x,y,col="red")
			points(x,y,pch=16,cex=0.5,col="red")
			lines(x_normal,y_normal,col="black")
			lines(x_boundary,y_boundary,col="blue",lwd=3)
			points(x_normal,y_normal,pch=16,cex=0.5,col="black")
			legend("bottomleft",legend=c("normal knot","barrier knot","normal edge","barrier edge","knot boundary","land"),lwd=c(NA,NA,1,1,3,1),pch=c(16,16,NA,NA,NA,NA),col=c("black","red","black","red","blue",NA),fill=c(NA,NA,NA,NA,NA,"gray90"),border=c(NA,NA,NA,NA,NA,"black"),bty="n",cex=1.5)

#_____________________________________________________________________________________________________________________________
# 3) fit a basic model
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
	print(paste0((B-A)[3]," seconds to fit spatiotemporal model"))
