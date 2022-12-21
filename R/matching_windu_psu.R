library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(tidyverse) 
library(sf)
library(data.table)
install.packages("RANN")
library(RANN)

wind_data <- nc_open(here("NetCDF files", "wrf_hi_WRF_Hawaii_Regional_Atmospheric_Model_best_2017_F.nc"))
print(wind_data) #look at variables and dimensions available in the dataset

lat.dat <- ncvar_get(wind_data, "lat")
long.dat <- ncvar_get(wind_data, "lon")
head(lat.dat) #units = "degrees_north"
head(long.dat) #units = "degrees_east"

time.dat <- ncvar_get(wind_data, "time")

wind_u <- ncvar_get(wind_data, "Uwind")
dim(wind_u) #should match long, lat, and time dims (213, 138, 2099)
wind_v <- ncvar_get(wind_data, "Vwind")

nc_close(wind_data) 

## time is hours since 2010-05-14 00:00:00.000 UTC
#time of prediction
new_time <- as.POSIXct(time.dat*60*60, origin = "2010-05-14", tz = "GMT") #convert hours to seconds with time.dat*60*60 and origin is date that the time is a measure from 
class(new_time)
head(new_time)


##psu shape file
psu <- st_read("./Data/GIS/Shapefiles/BFISH_PSU.shp")

camera_psu <- read.csv("Camera_PSU.csv")
rf_psu <- read.csv("Res_Fish_PSU.csv")
uni_psu <- unique(c(camera_psu$PSU, rf_psu$PSU))

psu <- psu %>% 
  filter(PSU %in% uni_psu) %>% 
  distinct()

#convert psu to dataframe then back to sf specfiying to use lat_deg and lon_deg as coordinates so the match up with wind crs
psu_sf <- psu %>% 
  as.data.frame() %>% 
  dplyr::select(OBJECTID, Island, PSU, lat_deg, lon_deg) %>% 
  st_as_sf(., coords = c("lon_deg", "lat_deg"), crs = "EPSG:4326")

psu.coords <- st_coordinates(psu_sf)




## find the closest lon and lat points from the wind data and psu
## nn.idx gives the index of which long or lat value is the nearest neighbor to the psu.coords row (1,2,3...2414) 
lon.nn <- nn2(as.numeric(long.dat), psu.coords[,1], k = 1, searchtype = "priority")$nn.idx
lat.nn <- nn2(as.numeric(lat.dat), psu.coords[,2], k = 1, searchtype = "priority")$nn.idx
#number of matrices in the array (number of sampling occaisons)
n.mat <- dim(wind_u)[3]
wind.dat <- c()
windu.fall <- list()
for(t in 1:length(n.mat)){
  
  ## get one time sampling period from wind data set
  windu.slice <- wind_u[, , t] 
  
  for(i in 1:length(lon.nn)){
    wind.dat[i] = windu.slice[lon.nn[i], lat.nn[i]]
    
  }
  
  windu.fall[[t]] <- data.frame(time = new_time[t], 
                           lon = long.dat[lon.nn],
                           lat = lat.dat[lat.nn],
                           windu = wind.dat,
                           psu.lon = psu.coords[,1],
                           psu.lat = psu.coords[,2])
}

## output is a list of dataframes of UWind data. Each dataframe in list is from 1 sampling time and within each dataframe is the lon and lat coords where the wind value was sampled from, the wind value, and the lon and lat coords of the psu closest to that wind sampling location.
## Next step is to match up psu sampling times with wind sampling times

### Create time column in rfc_psu df to match with CDF data
psu_time <- bind_rows(rf_psu, camera_psu) %>% 
  mutate(
    #add a zero in front of times with only 3 numbers (AM)
    start_time = ifelse(nchar(SAMPLE_START_TIME) == 3, 
                        str_pad(SAMPLE_START_TIME, 4, side = "left", pad = "0"), 
                        SAMPLE_START_TIME),
    #paste date and time together and convert to type POSIXl* format in HST time zone
    time = as.POSIXlt(paste(SAMPLE_DATE, start_time), 
                      format = "%Y-%m-%d %H%M", 
                      tz="US/Hawaii"),
    #subtracting 10 hours (10*60*60 = 36000) from time to convert to UTC to match wind time data
    time = time - 36000
  ) %>% 
  separate(time, into = c("date", "time.gmt"), sep = " ", remove = FALSE)




setkey(website, name, join_time)

psu_2017 <- psu_time %>% 
  filter(str_detect(BFISH, "2017_F")) %>% 
  dplyr::select(time, start_time, time.gmt, PSU, lat_deg, lon_deg) %>% 
  filter(!is.na(time)) %>% ## removed 152 rows that had no sample date or time
  mutate(join_time = time) %>% 
  as.data.table()

windu.fall.dt <- as.data.table(windu.fall)
windu.fall.dt$join_time <- windu.fall.dt$time

time.dt <- as.data.table(data.frame(time = new_time, join_time = new_time))

setkey(psu_2017, join_time)  
setkey(time.dt, join_time)  

matched.time.df <- time.dt[psu_2017, roll = T]

## Now can match a row of matched.time.df to a value in new_time and get the index. That would give you the index for wind.fall list to match with the psu sampling time and location








