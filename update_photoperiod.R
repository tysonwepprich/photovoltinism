# 3 minute runtime with 4 cores

# TODO:
# - automate new directory creation and naming, run this script multiple times
# - switch in diapause prediction (1 - ...) based on Fritzi or Dan regression coefs
#####
# packages, options, functions loaded
#####
library(sp)
library(rgdal)
library(raster)
library(lubridate)
library(mixsmsn)
library(dplyr)

# for parallel simulations with control over seed for reproducibility
# need different packages for windows computers
library(doRNG)
library(foreach) # for parallelized loops
if(.Platform$OS.type == "unix"){
  library(doMC)    # parallel backend for foreach, only for linux/mac
}else if(.Platform$OS.type == "windows"){
  library(doSNOW)
}
# Run model, but only with photoperiod update for speed
# The way I built the diapause decision, model results
# already track lifecycle through the end of GDD accumulation.
# To update photoperiod response, just need to take model
# results and rerun sens stage x photoperiod cue = diapause.

# DOESN'T WORK WITH WEIGHTED RESULTS
# NEED SEPARATE SUBSTAGE SIMS I THINK

source('CDL_funcs.R')
region_param <- "NW_SMALL"
REGION <- assign_extent(region_param = region_param)
# Substages with skewed t
arg1 = list(mu = 97.94, sigma2 = 2241.7, shape = 3.92, nu = 9.57)
x = seq(50, 350, length.out = 1000)
y = mixsmsn:::dt.ls(x, loc = arg1$mu, sigma2 = arg1$sigma2, shape = arg1$shape, nu = arg1$nu)
inputdist <- data.frame(x = x, y = y) %>% 
  arrange(x) %>% 
  mutate(CDF = cumsum(y/sum(y)))
substages <- SubstageDistrib(dist = inputdist, numstage = 7, perc = .99)
substages$means <- substages$means + 15

newdirs <- c("GCA_NWSMALL_2014", "GCA_NWSMALL_2015")
sitedirs <- "South_sol" # c("North_sol", "South", "South_sol")
params <- expand.grid(newdirs, sitedirs)
names(params) <- c("newname", "sitedir")

for (n in 1:nrow(params)){
  # newname <- "GCA_NWSMALL_NCDL_2016"
  newname <- as.character(params$newname[n])
  sitedir <- as.character(params$sitedir[n])
  coefs <- switch(substring(sitedir, 1, 1),
                  "N" = c(-56.9745, 3.5101),
                  "S" = c(-60.3523, 3.8888))
  if(length(grep(pattern = "sol", x = sitedir, fixed = TRUE)) == 1){
    solstice <- "after"
  }else{
    solstice <- "ignore"
  }
  
  gdd_file <- paste(newname, "dailygdd.grd", sep = "/")
  GDD <- brick(gdd_file)
  
  template <- GDD[[1]]
  template[!is.na(template)] <- 0
  template <- crop(template, REGION)
  
  # input sensitive stage
  # solstice <- "ignore" # "after" #
  sens <- "LS3"
  #GCA
  # coefs <- c(-56.9745, 3.5101) # Northern
  # coefs <- c(-60.3523, 3.8888) # Southern
  #DCA
  # coefs <- c(53.82741, -3.939534) # Big Bend State Park
  # coefs <- c(102.64569, -7.170386) # Delta
  # coefs <- c(83.85416, -5.967620) # Gold Butte
  # coefs <- c(98.44406, -6.910930) # Lovelock
  # coefs <- c(31.22615, -2.448147) # Topock Marsh
  # sitedir <- "South"
  
  returnwd <- getwd()
  setwd(newname)
  f <-list.files()
  rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]
  rasfiles <- rasfiles[-grep(pattern = "gdd", x = rasfiles, fixed = TRUE)]
  # for each sim with unique information to save
  ls <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,1])
  maps <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,2])
  sims <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,3])
  sims <- gsub(pattern = ".grd", replacement = "", x = sims)
  sims <- sims[which(sims != "weighted")]
  cdl_b0 <- coefs[1]
  cdl_b1 <- coefs[2]
  mapcode <- maps # for this case, no split map
  
  
  # parallel backend for foreach loop
  if(.Platform$OS.type == "unix"){
    ncores <- length(ls) * length(maps) / 2
    registerDoMC(cores = ncores)
  }else if(.Platform$OS.type == "windows"){
    ncores <- parallel::detectCores() - 1
    cl <- makeCluster(ncores)
    registerDoSNOW(cl)
  }
  
  system.time({
    foreach(sim = sims,
            .packages = "raster",
            .export = c("rasfiles", "template", "sitedir", 
                        "cdl_b0", "cdl_b1", "solstice")) %dopar%{
                          
                          # for (sim in sims){
                          oldfile <- rasfiles[grep(pattern = sim, x = rasfiles, fixed = TRUE)]
                          oldfile <- oldfile[grep(pattern = sens, x = oldfile, fixed = TRUE)]
                          
                          sens_stage <- brick(oldfile)
                          sens_stage <- sens_stage + template
                          LS4 <- sens_stage[[1]]
                          LS4[!is.na(LS4)] <- 0
                          if (exists("LS4stack")){
                            rm(LS4stack)
                          } 
                          
                          for (d in 1:nlayers(sens_stage)){
                            # doy <- lubridate::yday(lubridate::ymd(d))
                            doy <- as.numeric(gsub("layer.", replacement = "", x = names(sens_stage)[d]))
                            # can choose whether diapause induced before summer solstice
                            if (solstice == "after"){
                              if (doy > 170){
                                photo <- RasterPhoto(template, doy, perc_twilight = 25)
                                prop_diap <- 1 - exp(cdl_b0 + cdl_b1 * photo) /
                                  (1 + exp(cdl_b0 + cdl_b1 * photo))
                                tmpLS4 <- Cond(sens_stage[[d]] == 1, prop_diap, LS4)
                                LS4 <- Cond(prop_diap < LS4, LS4, tmpLS4)
                              }
                            }else{
                              photo <- RasterPhoto(template, doy, perc_twilight = 25)
                              prop_diap <- 1 - exp(cdl_b0 + cdl_b1 * photo) /
                                (1 + exp(cdl_b0 + cdl_b1 * photo))
                              tmpLS4 <- Cond(sens_stage[[d]] == 1, prop_diap, LS4)
                              LS4 <- Cond(prop_diap < LS4, LS4, tmpLS4)
                            }
                            
                            if (!exists("LS4stack")){
                              LS4stack <- stack(LS4)
                            } else {
                              LS4stack <- addLayer(LS4stack, LS4)
                            }
                          }
                          
                          LS4File <- writeRaster(LS4stack, filename = paste(sitedir,"/", "LS4_", "001_", sim, sep = ""),
                                                 overwrite = TRUE, datatype = "FLT8S")
                        }
  })
  
  
  stopCluster(cl) #WINDOWS
  setwd(returnwd)
}
