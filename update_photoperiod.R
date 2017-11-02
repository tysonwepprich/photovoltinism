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
region_param <- "SOUTHWEST"
gdd_file <- "dailygdd.grd"
REGION <- assign_extent(region_param = region_param)
newname <- "DCA_SW_2016"

# Substages with skewed t
arg1 = list(mu = 97.94, sigma2 = 2241.7, shape = 3.92, nu = 9.57)
x = seq(50, 350, length.out = 1000)
y = mixsmsn:::dt.ls(x, loc = arg1$mu, sigma2 = arg1$sigma2, shape = arg1$shape, nu = arg1$nu)
inputdist <- data.frame(x = x, y = y) %>% 
  arrange(x) %>% 
  mutate(CDF = cumsum(y/sum(y)))
substages <- SubstageDistrib(dist = inputdist, numstage = 7, perc = .99)
substages$means <- substages$means + 157

GDD <- brick(gdd_file)

template <- GDD[[1]]
template[!is.na(template)] <- 0
template <- crop(template, REGION)

# input sensitive stage
sens <- "LS3"
# coefs <- c(53.82741, -3.939534) # Big Bend State Park
coefs <- c(102.64569, -7.170386) # Delta
# coefs <- c(83.85416, -5.967620) # Gold Butte
# coefs <- c(98.44406, -6.910930) # Lovelock
# coefs <- c(31.22615, -2.448147) # Topock Marsh
Lifestage <- Lifestage + template

LS4 <- Lifestage[[1]]
LS4[!is.na(LS4)] <- 0

for (d in 1:nlayers(Lifestage)){
  # doy <- lubridate::yday(lubridate::ymd(d))
  doy <- as.numeric(gsub("layer.", replacement = "", x = names(Lifestage)[d]))
  photo <- RasterPhoto(template, doy, perc_twilight = 25)
  
  # sens_mask <- Cond(Lifestage %in% photo_sens, 1, 0)
  prop_diap <- 1 - exp(cdl_b0 + cdl_b1 * photo) /
    (1 + exp(cdl_b0 + cdl_b1 * photo))
  # tmpLS4 <- Cond(sens_mask == 1, prop_diap, LS4)
  
  tmpLS4 <- Lifestage[[d]] * prop_diap
  # LS4 <- tmpLS4 + LS4
  # need to account for prop_diap declining up until solstice
  # only let proportion diapausing increase over season
  LS4 <- Cond(tmpLS4 < LS4, LS4, LS4 + tmpLS4)
  if (!exists("LS4stack")){
    LS4stack <- stack(LS4)
  } else {
    LS4stack <- addLayer(LS4stack, LS4)
  }
}

plot(LS4stack[[seq(100, 290, 40)]])



returnwd <- getwd()
setwd(newname)

f <-list.files()
rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]

# for each sim with unique information to save
ls <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,1])
maps <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,2])
sims <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,3])
sims <- gsub(pattern = ".grd", replacement = "", x = sims)
sims <- sims[which(sims != "weighted")]

# parallel backend for foreach loop
if(.Platform$OS.type == "unix"){
  ncores <- length(ls) * length(maps) / 2
  registerDoMC(cores = ncores)
}else if(.Platform$OS.type == "windows"){
  ncores <- parallel::detectCores() / 2
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
}

# run inside foreach now
# tmppath <- paste0("~/REPO/photovoltinism/rastertmp/", "run", newname)
# dir.create(path = tmppath, showWarnings = FALSE)
# #sets temp directory
# rasterOptions(tmpdir=file.path(tmppath)) 

# this loop takes a lot of time

system.time({
  foreach(i = ls,
          .packages= "raster",
          .export = c("rasfiles", "substages")) %:% 
    foreach(m = maps,
            .packages = "raster",
            .export = c("rasfiles", "substages")) %dopar%{
              # .inorder = FALSE) %do%{
              # foreach(i = ls, .packages = "raster") %dopar% {
              tmppath <- paste0("~/REPO/photovoltinism/rastertmp/", i, "run", m)
              dir.create(path = tmppath, showWarnings = FALSE)
              #sets temp directory
              rasterOptions(tmpdir=file.path(tmppath)) 
              
              fs <- sort(rasfiles[grep(pattern = i, x = rasfiles, fixed = TRUE)])
              fs <- sort(fs[grep(pattern = m, x = fs, fixed = TRUE)])
              template <- brick(fs[1])[[1]]
              template[!is.na(template)] <- 0
              
              ll <- replicate(nlayers(brick(fs[1])), template)
              blank <- brick(ll)
              for (j in 1:length(fs)){
                ras_weighted <- brick(fs[j]) * substages[j, 2] # weights for each substage size
                blank <- overlay(blank, ras_weighted, fun=function(x,y) x + y)
              }
              outras <- writeRaster(blank, filename = paste(i, m, "weighted", sep = "_"),
                                    overwrite = TRUE)
              removeTmpFiles(h = 0)
              unlink(tmppath, recursive = TRUE)
              
            }
})

stopCluster(cl) #WINDOWS
# setwd(returnwd)