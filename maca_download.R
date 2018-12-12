# check out ento: ssh ento -l tyson -p 732
# shell call to copy macaV2 files (each 1.8G)
# scp -P 732 tyson@ento:/home/macav2metdata/IPSL_rcp85/macav2metdata_tas*_2046_2050_CONUS_daily.nc /home/tyson/REPO/photovoltinism/data/maca

library(raster)
library(ncdf4)
source('CDL_funcs.R')
source('species_params.R')

# download files from ENTO
system("scp -P 732 tyson@ento:/home/macav2metdata/IPSL_rcp85/macav2metdata_tas*_2016_2020_CONUS_daily.nc /home/tyson/REPO/photovoltinism/data/maca")

tmin <- brick("data/maca/macav2metdata_tasmin_IPSL-CM5A-MR_r1i1p1_rcp85_2046_2050_CONUS_daily.nc", 
              varname = "air_temperature", lvar = 3, level = 4)
tmax <- brick("data/maca/macav2metdata_tasmax_IPSL-CM5A-MR_r1i1p1_rcp85_2046_2050_CONUS_daily.nc", 
              varname = "air_temperature", lvar = 3, level = 4)

nc <- nc_open("data/maca/macav2metdata_tasmin_IPSL-CM5A-MR_r1i1p1_rcp85_2046_2050_CONUS_daily.nc")
v3 <- nc$var[[1]]
lonsize <- v3$varsize[1]
latsize <- v3$varsize[2] 
endcount <- v3$varsize[3] 

lon <- ncvar_get(nc,"air_temperature")

tmin <- tmin[[1]]
tmax <- tmax[[1]]
# longitude in degrees EAST
# time is days since 1900-01-01 00:00:00
# temperature in K? -273.15
# easier to adjust insect temperature threshold to kelvin than do extra addition on each raster

dd <- TriDD(tmax, tmin, LDT = 10 + 273.15, UDT = 37.8 + 273.15)





