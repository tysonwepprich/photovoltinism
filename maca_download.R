# check out ento: ssh ento -l tyson -p 732
# shell call to copy macaV2 files (each 1.8G)
# scp -P 732 tyson@ento:/home/macav2metdata/IPSL_rcp85/macav2metdata_tas*_2046_2050_CONUS_daily.nc /home/tyson/REPO/photovoltinism/data/maca

library(raster)
library(ncdf4)
library(RCurl)
source('CDL_funcs.R')
source('species_params.R')

# download files from ENTO
system("scp -P 732 tyson@ento:/home/macav2metdata/IPSL_rcp85/macav2metdata_tas*_2016_2020_CONUS_daily.nc /home/tyson/REPO/photovoltinism/data/maca")

# failures #####
d = pipe( 'ssh tyson@ento:732 "cat data.txt"' )
d = pipe( 'ssh -p 732 -l tyson -t -t ento:/home/macav2metdata/IPSL_rcp85/macav2metdata_tasmin_IPSL-CM5A-MR_r1i1p1_rcp85_2046_2050_CONUS_daily.nc')

tmin <- RCurl::scp("tyson@ento:732", 
                         "/home/macav2metdata/IPSL_rcp85/macav2metdata_tasmin_IPSL-CM5A-MR_r1i1p1_rcp85_2046_2050_CONUS_daily.nc", 
                         user = "tyson")
### 

tmin <- brick("data/maca/macav2metdata_tasmin_IPSL-CM5A-MR_r1i1p1_rcp85_2046_2050_CONUS_daily.nc", 
              varname = "air_temperature", lvar = 3, level = 4)
tmax <- brick("data/maca/2046/macav2metdata_tasmax_IPSL-CM5A-MR_r1i1p1_rcp85_2046_2050_CONUS_daily.nc", 
              varname = "air_temperature", lvar = 3, level = 4)

nc <- nc_open("tyson@ento: /home/macav2metdata/IPSL_rcp85/macav2metdata_tasmin_IPSL-CM5A-MR_r1i1p1_rcp85_2046_2050_CONUS_daily.nc")
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

macafile <- "data/maca/macav2metdata_tasmin_IPSL-CM5A-MR_r1i1p1_rcp85_2046_2050_CONUS_daily.nc"


# Function to split 5-year data into yearly tmin/tmax files for use in lifecycle model
SplitMACA <- function(macafile){
  ras <- brick(macafile, 
               varname = "air_temperature", lvar = 3, level = 4)
  
  yrs <- gregexpr(pattern = "[0-9]{4}", text = macafile)
  yr1 <- as.numeric(substr(macafile, start = yrs[[1]][1], stop = yrs[[1]][1] + 3))
  yr2 <- as.numeric(substr(macafile, start = yrs[[1]][2], stop = yrs[[1]][2] + 3))
  var <- sub(pattern = "as", replacement = "", x = stringr::str_split_fixed(macafile, pattern = "_", n = 3)[, 2], fixed = TRUE)
  
  # # account for leap years
  # days <- rep(floor(dim(ras)[3] / 5), length(c(yr1:yr2)))
  # days[which(c(yr1:yr2) %% 4 == 0)] <- 366
  
  for (y in yr1:yr2){
    newname <- paste0("MACAV2_", var, "_", y, ".grd")
    newras <- ras[[which(lubridate::year(as.Date(getZ(ras))) == y)]]
    dir.create(paste0("data/maca/", y))
    writeRaster(x = newras, filename = paste("data", "maca", y, newname, sep = "/"), format = "raster", overwrite = TRUE)
  }
}

SplitMACA(macafile)


