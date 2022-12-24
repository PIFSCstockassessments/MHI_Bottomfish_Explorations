

# Nicholas Ducharme-Barth
# 12/22/2022
# Set up single species cross validated data sets
# data set, Extrapolation_List, and spatial_list

# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(ggplot2)
	library(ggthemes)
    library(VAST)
	# rgdal functions sourced directly
	# sp functions sourced directly
	# DHARMa functions sourced directly
	# sf functions sourced directly
	# ggthemes functions sourced directly
#_____________________________________________________________________________________________________________________________
# specify (coarse) model configuration
    data_flag = 2021
    species = "mv" # or "mv"
    data_treatment = "05"
    lehi_filter = TRUE
    km_cutoff = 7.5 # make this smaller to increase the spatial resolution of the model
    fine_scale = FALSE
    k_folds = 10
    seed = 123

    model_name = paste(data_flag,
                species,
                data_treatment,
                lehi_filter,
                km_cutoff,
                fine_scale,
                k_folds,
                seed,
                sep="_")

#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"
	working_dir = paste0(proj.dir,"VAST/xval_data/",model_name,"/")
	dir.create(working_dir,recursive=TRUE)

#_____________________________________________________________________________________________________________________________
# 1) bring in data

    if(species == "single")
    {
        load(file=paste0(proj.dir,"Data/",data_flag,"_",data_treatment,".bfish_combined_wide_dt.RData"))
        main_bfish_df = bfish_combined_wide_dt %>% as.data.frame(.)
    } else if(species == "mv"){
        load(file=paste0(proj.dir,"Data/",data_flag,"_",data_treatment,".bfish_combined_long_dt.RData"))
        main_bfish_df = bfish_combined_long_dt %>% as.data.frame(.)   
    }

	# remove sample with large lehi observation
    if(lehi_filter)
    {
        main_bfish_df =  subset(main_bfish_df,design_sampling_unit!="2021_Fall_32293")
    }

#_____________________________________________________________________________________________________________________________
# 2) partition
    set.seed(seed)
    u_sample_dt = data.table(u_sample = unique(main_bfish_df$model_sampling_unit))
    u_sample_dt$partition = sample(1:k_folds,size=nrow(u_sample_dt),replace=TRUE)

#_____________________________________________________________________________________________________________________________
# 3) iterate over partitions
# save data, Extrapolation_List, and spatial_list

    # processing that only needs to be done once
    	# needed to define spatial domain and for predicting on to create index
        psu_table = fread(paste0(proj.dir,"Data/BFISH PSU lookup table.csv")) %>%
                    .[,.(PSU,Island,lon_deg,lat_deg,STRATA,STRATA_2020,Depth_MEDIAN_m,med_slp,med_acr,BS_pct_over_136j,pctHB,pctHS)] %>%
                    .[,substrate:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][1])] %>%
                    .[,slope:=sapply(STRATA,function(x)strsplit(x,"_")[[1]][2])]
        # needed for plotting and defining the barrier feature
        hi_coast = rgdal::readOGR(dsn = paste0(proj.dir,"Data/GIS/Coastline"), layer = "Coastline")
        hi_coast_sf = sf::st_as_sf(hi_coast)
        # convert to equal distant projection
            # get original lat-lon crs
            crs_ll = sp::CRS(sp::proj4string(hi_coast))
            # use two-point equi-distant projection
            crs_eqd = sp::CRS("+proj=tpeqd +lat_1=20.45 +lon_1=-158.9 +lat_2=20.45 +lon_2=-156.2 +datum=WGS84 +ellps=WGS84 +units=km +no_defs")
            hi_coast_eqd = sp::spTransform(hi_coast, crs_eqd)
            hi_coast_eqd_sf = sf::st_as_sf(hi_coast_eqd) 
        	psu_coords = cbind(psu_table$lon_deg, psu_table$lat_deg)
            psu_sp = sp::SpatialPoints(psu_coords)
            sp::proj4string(psu_sp) = crs_ll
            psu_sp_eqd = sp::spTransform(psu_sp, crs_eqd)
            psu_table$lon_eqd = psu_sp_eqd@coords[,1]
            psu_table$lat_eqd = psu_sp_eqd@coords[,2]
        # define extrapolation grid based on PSUs
		    input_grid = cbind(psu_coords[,2],psu_coords[,1],0.5^2)
		    colnames(input_grid) = c("Lat","Lon","Area_km2")
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

for(k in 1:k_folds)
{
    train_samples = u_sample_dt[partition!=k]$u_sample
    test_samples = u_sample_dt[partition==k]$u_sample

    bfish_df = subset(main_bfish_df,model_sampling_unit %in% train_samples)
    test_df = subset(main_bfish_df,model_sampling_unit %in% test_samples)

    # convert to equi-distant coordinates for samples & psu
		bfish_coords = cbind(bfish_df$lon, bfish_df$lat)
		bfish_sp = sp::SpatialPoints(bfish_coords)
		sp::proj4string(bfish_sp) = crs_ll
		bfish_sp_eqd = sp::spTransform(bfish_sp, crs_eqd)
		bfish_df$lon_eqd = bfish_sp_eqd@coords[,1]
		bfish_df$lat_eqd = bfish_sp_eqd@coords[,2]
    # define Extrapolation_List
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
    # define spatial_list
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
    
    # save
        xval_extrapolation_list = Extrapolation_List
        xval_spatial_list = spatial_list
        xval_train_df = bfish_df
        xval_test_df = test_df

        save(xval_extrapolation_list,file=paste0(working_dir,k,"_xval_extrapolation_list.RData"))
        save(xval_spatial_list,file=paste0(working_dir,k,"_xval_spatial_list.RData"))
        save(xval_train_df,file=paste0(working_dir,k,"_xval_train_df.RData"))
        save(xval_test_df,file=paste0(working_dir,k,"_xval_test_df.RData"))

    # clean-up
        rm(list=c("train_samples","test_samples","bfish_df","test_df",
            "bfish_coords","bfish_sp","bfish_sp_eqd",
            "Extrapolation_List","a_el_orig","a_el_tmp",
            "spatial_list",
            "xval_extrapolation_list","xval_spatial_list","xval_train_df","xval_test_df"))
}

