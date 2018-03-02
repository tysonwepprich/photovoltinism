# Condensed DDRP model
# Universal code for 3 biocontrol species
# Any number of lifestages allowed
library(profvis)

profvis({
  testtime <- system.time({
#####
# packages, options, functions loaded
#####
library(sp)
library(rgdal)
library(raster)
library(lubridate)
library(mixsmsn)
library(dplyr)
library(stringr)
library(purrr)

# for parallel simulations 
library(foreach) # for parallelized loops
library(doParallel) # this package can be used across operating systems

# adjust raster options for your computer
rasterOptions(overwrite = FALSE, 
              chunksize = 1e+07,
              maxmemory = 1e+08)

# load collection of functions for this model
source('CDL_funcs.R') 

# directory with a year of PRISM tmax and tmin files
prism_path <- "prismDL/2014"
# prism_path <- "/data/PRISM/2014"
# prism_path <- "/data/PRISM/"



#####
#input parameters
#####
# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model extext
yr <- 2014
start_doy  <- 1
end_doy    <- 365
region_param <- "NORTHWEST"

# life cycle parameters
# stgorder   <- c("OA","E","L","P","A","F")
stgorder   <- c("OA","E","L","P", "A") # reassigned as 1, 2, 3, 4, 5 in model
# Degree day thresholds
# need to match length and order of stgorder
stage_ldt <- rep(10, 5)
stage_udt <- rep(37.8, 5)
stage_dd <- c(167.6, 93.3, 136.4, 137.7, 125.9)

# diapause response to photoperiod
photo_sens <- 5 #c(-1, 3) # integer life stages for now
CDL_mu        <- 16.23 # 15.52
model_CDL  <- 1 # if 1, model photoperiod decision
CDL_log    <- 1 # if 1, model CDL from logistic regression results

# Logistic regression photoperiod
# Careful to check if regression models diapause percent or reproduce percent as output
# from Fritzi's lab data, two different populations
coefs <- c(-56.9745, 3.5101) # Northern
# coefs <- c(-60.3523, 3.8888) # Southern

# GDD data and calculation
gdd_data <- "calculate" # c("load", "calculate")
# gdd_file <- "dailygdd_2017_WEST.grd"
calctype   <-"triangle"

# introducing individual variation, tracked with simulation for each substage
vary_indiv <- 1 # turn on indiv. variation
if (vary_indiv == 1){
  nsim <- 4 # number of substages
}else{
  nsim <- 1
}

#####
# derived parameters
#####
REGION <- assign_extent(region_param = region_param)

nday <- length(start_doy:end_doy)

# Take empirical distribution and calculate substages
# OW oviposition distribution
# # TRY 1 with beta
# eggdist <- dbeta(x = seq(0, 1, length.out = 1000), 
#                  shape1 = 3.888677, shape2 = 2.174208)
# inputdist <- data.frame(x = seq(59.6, 223.3677, length.out = 1000),
#                         y = eggdist)
# inputdist$CDF <- cumsum(inputdist$y) / sum(inputdist$y, na.rm = TRUE)
# 
# substages <- SubstageDistrib(dist = inputdist, numstage = nsim, perc = .99)

# Try 2 with skewed t
arg1 = list(mu = 97.94, sigma2 = 2241.7, shape = 3.92, nu = 9.57)
x = seq(50, 350, length.out = 1000)
y = mixsmsn:::dt.ls(x, loc = arg1$mu, sigma2 = arg1$sigma2, shape = arg1$shape, nu = arg1$nu)
inputdist <- data.frame(x = x, y = y) %>% 
  arrange(x) %>% 
  mutate(CDF = cumsum(y/sum(y)))
substages <- SubstageDistrib(dist = inputdist, numstage = 4, perc = .99)

# To get observations to fit for overwinter adults and F1 eggs, 
# substages$means <- substages$means + 15


# try to run sims in parallel and also split map if large REGION
if (region_param == "CONUS"){
  ncores <- nsim * 2 # trying out half the cores needed, see if writeRaster better
}else{
  ncores <- nsim
}

# avoid memory issues on laptop
if(ncores > (parallel::detectCores() / 2)){
  ncores <- parallel::detectCores() / 2
}
cl <- makePSOCKcluster(ncores)
registerDoParallel(cl)

# still here for day naming in loop, but not necessary
# figure out how to cut
pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
files <- list.files(path = prism_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = TRUE)
#Check that there are enough files for the year
length(files)
numlist <- vector()
for (file in files) {
  num = strsplit(file,split="_")[[1]][5]
  numlist <-c(numlist,num)
}

sortedlist <- sort(numlist)
template <- raster(files[1])
template <- crop(raster(files[1]), REGION)
# template <- crop(aggregate(raster(files[1]), fact = 2), REGION)

template[!is.na(template)] <- 0



# If GDD not pre-calculated
pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
tminfiles <- list.files(path = prism_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = TRUE)
tminfiles <- ExtractBestPRISM(tminfiles, yr)
r <- crop(raster(tminfiles[1]), REGION)
tminstack <- stack(r)
for (i in tminfiles[-1]){
  tminstack <- addLayer(tminstack, crop(raster(i), REGION))
}
# Robert Hijmans hack, not useful with crop
# tminstack@layers <- sapply(tminfiles, function(x) { r@file@name=x; r } ) 

pattern = paste("(PRISM_tmax_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
tmaxfiles <- list.files(path = prism_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = TRUE)
tmaxfiles <- ExtractBestPRISM(tmaxfiles, yr)
r <- crop(raster(tmaxfiles[1]), REGION)
tmaxstack <- stack(r)
for (i in tmaxfiles[-1]){
  tmaxstack <- addLayer(tmaxstack, crop(raster(i), REGION))
}

# tmaxstack <- stack(r)
# tmaxstack@layers <- sapply(tmaxfiles, function(x) { r@file@name=x; r } )

# really slow on full stack
# tminstack <- crop(tminstack, REGION)
# tmaxstack <- crop(tmaxstack, REGION)


# if (region_param == "CONUS"){
#   # splits template into list of smaller map rasters to run in parallel
#   # only doing this for CONUS, because benefit lost for smaller map sizes
#   SplitMap <- SplitRas(template, ppside = 2, TRUE, FALSE)
#   # runs_par <- expand.grid(1:nsim, 1:length(SplitMap))
#   # names(runs_par) <- c("sim", "map")
# }else{
SplitMap <- list(template)
runs_par <- expand.grid(1:nsim, 1:length(SplitMap))
names(runs_par) <- c("sim", "map")
# }

# parameters of required degree-days and two coefficients of photoperiod response model
params <- cbind(substages$weights, substages$means, matrix(rep(c(stage_dd[-1], coefs), nrow(substages)), nrow = nrow(substages), byrow = TRUE))


# create directory for model output
newname <- paste("run", gsub("-", "", 
                             gsub(" ", "_", 
                                  gsub(":", "", Sys.time(), fixed = TRUE),
                                  fixed = TRUE),
                             fixed = TRUE), sep = "")
dir.create(newname)

# if errors with some sims/maps, rerun foreach loop below again for them
# use folders left in raster tmp directory to know which didn't work

#Run model
############################

# system.time({
  outfiles <- foreach(sim = 1:nsim,
                      # outfiles <- foreach(sim = 5, # if some runs don't work, rerun individually
                      .packages= "raster",
                      .export = c("SplitMap", "tmaxstack", "tminstack",
                                  "newname", "params"),
                      .inorder = FALSE) %:% 
    foreach(map = 1:length(SplitMap),
            # foreach(map = 3, # if some runs don't work, rerun individually
            .packages = "raster",
            .export = c("SplitMap", "tmaxstack", "tminstack",
                        "newname", "params"),
            .inorder = FALSE) %dopar%{
              # .inorder = FALSE) %do%{
              
              # profvis({
              # creates unique filepath for temp directory to store rasters
              tmppath <- paste0("~/REPO/photovoltinism/rastertmp/", "run", map, sim)
              dir.create(path = tmppath, showWarnings = FALSE)
              #sets temp directory
              rasterOptions(tmpdir=file.path(tmppath)) 
              
              template <- SplitMap[[map]]
              # tmpGDD <- crop(GDD, template)
              # tmpGDD <- tmpGDD[[start_doy:end_doy]]
              
              #### All Stages at once
              # varying ldt and udt determined above in parameter file
              stage_dd <- params[sim, 1:length(stgorder)]
              cdl_b0 <- params[sim, length(stgorder)+1]
              cdl_b1 <- params[sim, length(stgorder)+2]
              
              dd_accum <- template
              lifestage <- template + 1 # Need for OW designated as stage 1
              numgen <- template
              diap <- template
              # }) # profvis
              # loop over days
              sublist <- sortedlist[start_doy:end_doy]
              for (d in sublist) {
                # profvis({
                print(d)
                index <- which(sublist == d)
                # if (.Platform$OS.type == "windows") flush.console()
                # Sys.sleep(1)
                
                # profvis({
                tmin <- tminstack[[index]]
                tmax <- tmaxstack[[index]]
                
                if(region_param == "CONUS"){
                  tmin <- crop(tmin, template)
                  tmax <- crop(tmax, template)
                }
                
                
                # photoperiod for this day across raster
                # could be done outside for loop for speed up
                if (model_CDL == 1){
                  doy <- lubridate::yday(lubridate::ymd(d))
                  photo <- RasterPhoto(template, doy, perc_twilight = 25)
                }
                
                ls_ldt <- setValues(lifestage, stage_ldt[getValues(lifestage)])
                ls_udt <- setValues(lifestage, stage_udt[getValues(lifestage)])
                ls_dd  <- setValues(lifestage, stage_dd[getValues(lifestage)])
                
                dd_tmp <- overlay(tmax, tmin, ls_ldt, ls_udt, fun = TriDD)
                
                #Accumulate degree days
                # dd.stat <- cellStats(dd_tmp,stat='max',na.rm=TRUE)
                # if (dd.stat > 0) {
                  dd_accum <- dd_accum + dd_tmp
                  #Calculate lifestage progression: for dd accum in correct lifestage, is it >= egg DD threshold?
                  progress <- dd_accum >= ls_dd
                  lifestage <- lifestage + progress
                  #Reset the DDaccum cells to zero for cells that progressed to next lifestage
                  dd_accum <- dd_accum - (progress * ls_dd)
                  ovip <- Cond(lifestage == (which(stgorder == "A") + 1), 1, 0)
                  numgen <- numgen + ovip
                  lifestage <- Cond(ovip == 1, which(stgorder == "E"), lifestage)
                # }
                
                if (model_CDL == 1){
                  # add logistic regression variation in CDL response
                  # keep running lifecycle as if all directly developing
                  # but add new raster to track percent of population actually
                  # in diapause following logistic regression
                  sens_mask <- Cond(lifestage %in% photo_sens, 1, 0)
                  prop_diap <- 1 - exp(cdl_b0 + cdl_b1 * photo) /
                    (1 + exp(cdl_b0 + cdl_b1 * photo))
                  tmpdiap <- Cond(sens_mask == 1, prop_diap, diap)
                  # need to account for prop_diap declining up until solstice
                  # only let proportion diapausing increase over season
                  diap <- Cond(prop_diap < diap, diap, tmpdiap)
                  
                  # memory mgmt, but slow step in computation
                  # rm(sens_mask, prop_diap, tmpdiap)
                  # gc()
                }
                
                
                # Add section to combine simulations into weighted values before stacking
                # Get one raster per day per 3 outputs first, fewer files to deal with later
                # But how to parallelize?
                
                # stack each days raster
                if (!exists("LS_stack")){
                  # hope that this will keep this large files on disk, out of memory load
                  LS_stack <- stack(lifestage)
                } else {
                  LS_stack <- addLayer(LS_stack, lifestage)
                  
                  # #alternatives
                  # # too slow
                  #  LS0stack <- stack(LS0stack, LS0)
                  # # faster but dangerous, doesn't check extent in each file
                  #  LS0stack@layers[[index]] <- LS0
                }
                
                if (!exists("diap_stack")){
                  diap_stack <- stack(diap)
                } else {
                  diap_stack <- addLayer(diap_stack, diap)
                }    
                
                
                if (!exists("numgen_stack")){
                  numgen_stack <- stack(numgen)
                } else {
                  numgen_stack <- addLayer(numgen_stack, numgen)
                }
                
                # }, interval = 0.005) #profvis
              }  # daily loop
              
              
              # if parallel by map chunk
              # each band is a day, first index is map chunk index
              mapcode <- formatC(map, width = 3, format = "d", flag = "0")
              
              NumGenFile <- writeRaster(numgen_stack,
                                        filename = paste(newname, "/NumGen_", mapcode, "_sim", sim, sep = ""),
                                        overwrite = TRUE, datatype = "INT1U")
              LSFile <- writeRaster(LS_stack, filename = paste(newname, "/LS_", mapcode, "_sim", sim, sep = ""),
                                    overwrite = TRUE, datatype = "INT1U")
              LS4File <- writeRaster(diap_stack, filename = paste(newname, "/diap_", mapcode, "_sim", sim, sep = ""),
                                     overwrite = TRUE, datatype = "FLT4S")
              
              # memory management
              rm(numgen_stack, LS_stack, diap_stack)
              #removes entire temp directory without affecting other running processes
              unlink(tmppath, recursive = TRUE)
              gc()
                # }, interval = 0.005) #profvis
                
              # keep track of written output files
              out <- lapply(X = list(NumGenFile,
                                     LSFile,
                                     LS4File), filename)
            } # close foreach loop
# }) #system.time

# if(.Platform$OS.type == "windows"){
  stopCluster(cl)
# }
}) # system.time
}) #profvis