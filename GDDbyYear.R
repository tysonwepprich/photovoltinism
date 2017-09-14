# this outputs average GDD per day over 2007-2013, 
# to match models in Grevstad & Coop 2015



library(sp)
library(rgdal)
library(raster)
source('CDL_funcs.R')

rasterOptions(overwrite = FALSE, 
              chunksize = 1e+08,
              maxmemory = 1e+09)

# for parallel simulations with control over seed for reproducibility
library(doRNG)
library(foreach) # for parallelized loops
library(doMC) 

region_param <- "CONUS"

LDT <- 10
UDT <- 37.8

base_path <- "/data/PRISM/"
years <- 2007 #c(2007:2013)
ncores <- length(years)
registerDoMC(cores = ncores)

foreach(yr = years, .packages= "raster")   %dopar% {
  prism_path <- paste(base_path, yr, sep = "/")
  
  #Search pattern for PRISM daily temperature grids. Load them for processing.
  pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
  tminfiles <- list.files(path = prism_path, pattern=pattern, all.files=FALSE, full.names=TRUE)
  r <- raster(tminfiles[1])
  tminstack <- stack(r)
  tminstack@layers <- sapply(tminfiles, function(x) { r@file@name=x; r } ) 
  
  pattern = paste("(PRISM_tmax_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
  tmaxfiles <- list.files(path = prism_path, pattern=pattern, all.files=FALSE, full.names=TRUE)
  r <- raster(tmaxfiles[1])
  tmaxstack <- stack(r)
  tmaxstack@layers <- sapply(tmaxfiles, function(x) { r@file@name=x; r } ) 
  

  REGION <- switch(region_param,
                   "CONUS"        = extent(-125.0,-66.5,24.0,50.0),
                   "NORTHWEST"    = extent(-125.1,-103.8,40.6,49.2),
                   "OR"           = extent(-124.7294, -116.2949, 41.7150, 46.4612),
                   "TEST"         = extent(-124, -122.5, 44, 45))
  
  tminstack <- crop(tminstack, REGION)
  tmaxstack <- crop(tmaxstack, REGION)

  # system.time({
  GDD <- overlay(x = tmaxstack, y = tminstack, 
                 fun = function(x, y){
                   Tmp1=6*((x-LDT)*(x-LDT))/(x-y)
                   Tmp2=6*((x-UDT)*(x-UDT))/(x-y)
                   Cond(x < LDT,0,
                        Cond(y >= UDT,UDT-LDT,
                             Cond((x < UDT) & (y <= LDT), Tmp1/12,
                                  Cond((y <= LDT) & (x >= UDT), (Tmp1-Tmp2)/12,
                                       Cond((y > LDT) & (x >= UDT), 6*(x+y-2*LDT)/12 - (Tmp2/12),
                                            Cond((y > LDT) & (x < UDT), 6*(x+y-2*LDT)/12,0))))))
                 },
                 filename = paste("dailygdd", yr, sep = "_"),
                 recycle = FALSE,
                 overwrite = TRUE)
  # })

}


rasfiles <- list.files(pattern = "dailygdd")
rasfiles <- rasfiles[grep(pattern = ".grd", x = rasfiles, fixed = TRUE)]

# average GDD per day/raster
# account for leap years

# overlay? calc? mean?

# maybe extract values as matrix, make array, then apply mean over dimension
test1 <- stack(rasfiles[7])
test2 <- brick(rasfiles[2])[[1:365]] # remove leap year

test3 <- mean(test1, test2, na.rm = TRUE)

test <- as.array(brick(rasfiles[1]))

theATs <- lapply(rasfiles, FUN = function(x) brick(x)[[1:365]])
res <- Reduce("+", theATs, accumulate = TRUE)
# last rasterbrick in res in final accumulated sum
meanGDD <- res[[length(res)]] / length(res)
writeRaster(x = meanGDD, filename = "meanGDD_07_13")