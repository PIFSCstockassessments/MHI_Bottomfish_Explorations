
dat<-read.csv("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\CPUE standardization\\Wind\\wrf_hi_6258_f657_ef41.csv") #more than 30 min to load
                                                                                                                    #should look into geting these as binary or something....didn't marc say previously?
dim(dat)

names(dat)


#################################
#as Netcdf binary###############
################################

rm(list = ls())

library(chron)
library(RNetCDF)


wind_data=open.nc("C:\\Users\\John.Syslo\\Documents\\2023 Deep 7\\CPUE standardization\\Wind\\wrf_hi_WRF_Hawaii_Regional_Atmospheric_Model_best_2020.nc")
print.nc(wind_data)
var.inq.nc(wind_data,time)

time_data=read.nc(wind_data)$time
time_run_data=read.nc(wind_data)$time_run

#time units are ="hours since 2010-05-14 00:00:00.000 UTC" ;
#explore time dimension
dim(time_data)
#[1] 8785
dim(time_data)/365 #about every 24 hrs, as described
#[1] 24.06849
length(unique(time_data))
#[1] 8785
length(unique(time_run_data)) #run about once per day, safe to assume evrything is from a 24-hour forecast (they go up to 7-days)?
#[1] 359
length((time_run_data)) 
#[1] 8785


#use POSIX date values, convert 2010-05-14 00:00:00.000 UTC to seconds since January 1, 1970
start<-1273795200
#class(start)=c('POSIXt','POSIXct')


#time of prediction
new_time<-start+time_data*60*60 #add thenumber of seconds since start time
class(new_time)=c('POSIXt','POSIXct')
wind_data$time_2<-new_time
  
#head(wind_data$time_2)

#run time
new_time_run<-start+time_run_data*60*60
class(new_time_run)=c('POSIXt','POSIXct')
wind_data$time_run_2<-new_time_run

#head(wind_data$time_run_2)


#need to re-evaluate the following.  how to convert between u-vector and v-vector to to wind speed and direction...
# ws = sqrt(u2+v2 )
#direction = atan2(y,x) or atan2(v,u)

#wind_data$u_data=read.nc(wind_data)$Uwind
#wind_data$v_data=read.nc(wind_data)$Vwind

u_data=read.nc(wind_data)$Uwind
v_data=read.nc(wind_data)$Vwind

#creat new variables for speed and direction
#speed<-sqrt(u_data^2+v_data^2 )
#wind_data$speed<-speed

wind_data$speed<-sqrt(wind_data$Uwind^2+wind_data$Vwind^2 )

#THIS ONE IS CORRECT
wind_data$direction<-(270-atan2(v_data,u_data)*180/pi) %% 360      #https://www.eol.ucar.edu/content/wind-direction-quick-reference #"feels" more accurate, but not sure about basis for 270 -
#https://stackoverflow.com/questions/8673137/calculating-wind-direction-from-u-and-v-components-of-the-wind-using-lapply-or-i:

#Check out:
  
#  https://www.eol.ucar.edu/content/wind-direction-quick-reference

#In short, you want to use atan2 to take care of the different quadrants, but don't forget that the angles are defined differently in meteorology! We measure the wind direction in degrees clockwise from North, whereas atan2-type functions generally work in radians from the
#X-direction (i.e. East). So you want to use something like:

#WDIR= 270-atan2(V,U)*180/pi

#other references exist about 0 being along the x-axis...hence 270


head(wind_data$direction)
summary(wind_data$direction)
hist(wind_data$direction) #slightly more north than I would expect (without 270), but this is 1 year of data....



#code from previous wind data###############


#xdir_data=read.nc(wind_data)$u
#ydir_data=read.nc(wind_data)$v
#speed_data=read.nc(wind_data)$w



#this will need to be overhauled
#time.detail=month.day.year(time_data)
#vals_y=which(time.detail$year==month.day.year(time_data/24, c(1,1,1978))$year)
#vals_m=which(time.detail$month==month.day.year(time_data/24, c(1,1,1978))$month)
#vals_d=which(time.detail$day==month.day.year(time_data/24, c(1,1,1978))$day)
#date_loc=intersect(intersect(vals_y,vals_m),vals_d)

#try to figure out the dimnames
#character(file.inq.nc(wind_data)$ndims) #just blank
