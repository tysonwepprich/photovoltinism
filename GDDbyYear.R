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
# library(doMC) 
library(doSNOW) # for WINDOWS

region_param <- "EAST"

LDT <- 10
UDT <- 37.8

# base_path <- "/data/PRISM/"
base_path <- "prismDL"

years <- c(2014:2017)
# LINUX
# ncores <- length(years)
# registerDoMC(cores = ncores)

# WINDOWS
ncores <- length(years)
cl <- makeCluster(ncores)
registerDoSNOW(cl)

foreach(yr = years, .packages= "raster")   %dopar% {
  prism_path <- paste(base_path, yr, sep = "/")
  
  #Search pattern for PRISM daily temperature grids. Load them for processing.
  pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
  tminfiles <- list.files(path = prism_path, pattern=pattern, 
                          all.files=FALSE, full.names=TRUE, recursive = TRUE)
  dates <- regexpr(pattern = "[0-9]{8}", text = tminfiles)
  fileorder <- order(regmatches(tminfiles, dates))
  tminfiles <- tminfiles[fileorder]
  r <- raster(tminfiles[1])
  tminstack <- stack(r)
  tminstack@layers <- sapply(tminfiles, function(x) { r@file@name=x; r } ) 
  
  pattern = paste("(PRISM_tmax_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
  tmaxfiles <- list.files(path = prism_path, pattern=pattern, 
                          all.files=FALSE, full.names=TRUE, recursive = TRUE)
  dates <- regexpr(pattern = "[0-9]{8}", text = tmaxfiles)
  fileorder <- order(regmatches(tmaxfiles, dates))
  tmaxfiles <- tmaxfiles[fileorder]
  r <- raster(tmaxfiles[1])
  tmaxstack <- stack(r)
  tmaxstack@layers <- sapply(tmaxfiles, function(x) { r@file@name=x; r } ) 
  

  REGION <- assign_extent(region_param)
  
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
                 filename = paste("dailygdd", yr, region_param, sep = "_"),
                 recycle = FALSE,
                 overwrite = TRUE)
  # })

} # end foreach loop
stopCluster(cl) #WINDOWS

# 
# # plot prism rasterbrick
# # testing
# res <- tmaxstack
# names(res) <- paste("d", formatC(1:nlayers(res), width = 3, format = "d", flag = "0"), sep = "")
# 
# gdd <- raster::extract(res, y = sites[, 2:3])
# gdd <- cbind(sites, gdd)
# gdd <- gdd %>% 
#   tidyr::gather(key = "DOY", value = "GDD", d001:d294) %>% 
#   dplyr::mutate(DOY = as.numeric(gsub(pattern = "d", replacement = "", x = .$DOY))) %>% 
#   group_by(ID) %>% 
#   arrange(DOY) %>% 
#   dplyr::mutate(Accum_GDD = cumsum(GDD))
# 
# plt <- ggplot(gdd, aes(x = DOY, y = GDD, group = ID)) +
#   geom_line(size = 2) +
#   facet_wrap(~ID, scales = "free_y")
# plt




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
