# Take Gerricke's FCM code and convert to Galerucella parameters
# Add photoperiod variation and substages

# TODO
# -ARRANGE PRISM FILES IN ORDER AS CHECK (2017 off with naming differences)
# -Output has strange shoreline edges with no results. Not NA, but
# also does not advance lifestages. Issue with PRISM or raster?
# -Pest Event Maps
# -Diapause stage leading to overwintering stage (using same raster?)
# -Parallel amended to use doParallel package wrapper for cross platform
# -Frost events at start/end of season. Track annual variation in season length.
# -Twilight assumption plays big role in potential voltinism


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

rasterOptions(overwrite = FALSE, 
              chunksize = 1e+07,
              maxmemory = 1e+08)

source('CDL_funcs.R') # load collection of functions for this model

prism_path <- "prismDL/2017"
# prism_path <- "/data/PRISM/2014"
# prism_path <- "/data/PRISM/"



#####
#input parameters
#####
# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model extext
start_doy  <- 1
end_doy    <- 294
region_param <- "WEST"
# life cycle parameters
stgorder   <- c("OA","E","L","P","A","F")
photo_sens <- 3 #c(-1, 3) # integer life stages for now
CDL_mu        <- 16.23 # 15.52
model_CDL  <- 1 # if 1, model photoperiod decision
CDL_log    <- 1 # if 1, model CDL from logistic regression results
owstage    <- "OA"
# Logistic regression photoperiod
# from Fritzi's lab data, two different populations
coefs <- c(-56.9745, 3.5101) # Northern
# coefs <- c(-60.3523, 3.8888) # Southern
# Degree day thresholds
# LDT = lower development threshold, temp at which growth = 0 (using PRISM tmean)
eggLDT     <- 10
eggUDT     <- 37.8  
larvaeLDT  <- 10
larvaeUDT  <- 37.8  
pupaeLDT   <- 10
pupaeUDT   <- 37.8
adultLDT   <- 10 
adultUDT   <- 37.8
# Degree day requirements for life stages
# DD = degree days, number of cumulative heat units to complete that lifestage
OWadultDD_mu  <- 167.6 #based on Oregon counts
  #or 108 based on McAvoy 1997 phenology of 3 years
eggDD_mu = 93.3
larvaeDD_mu = 136.4 # 46.9 + 45.8 + 43.7 instars
pupDD_mu = 137.7 
adultDD_mu = 125.9 #time to oviposition
# GDD data and calculation
gdd_data <- "load" # "calculate"
gdd_file <- "dailygdd_2017_WEST.grd"
calctype   <-"triangle"
# introducing individual variation, tracked with simulation for each substage
vary_indiv <- 1 # turn on indiv. variation
if (vary_indiv == 1){
  nsim <- 7 # number of substages
}else{
  nsim <- 1
}

#####
# derived parameters
#####
REGION <- switch(region_param,
                 "CONUS"        = extent(-125.0,-66.5,24.0,50.0),
                 "NORTHWEST"    = extent(-125.1,-103.8,40.6,49.2),
                 "OR"           = extent(-124.7294, -116.2949, 41.7150, 46.4612),
                 "TEST"         = extent(-124, -122.5, 44, 45),
                 "WEST"         = extent(-125.14, -109, 37, 49.1))

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
substages <- SubstageDistrib(dist = inputdist, numstage = 7, perc = .99)
# To get observations to fit for overwinter adults and F1 eggs, 
# overwinter pre-oviposition period is only 50 deg days
substages$means <- substages$means + 15


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

if(.Platform$OS.type == "unix"){
  registerDoMC(cores = ncores)
}else if(.Platform$OS.type == "windows"){
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
}

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
template[!is.na(template)] <- 0


# preload PRISM data into RasterBrick
if (gdd_data == "load"){
  GDD <- brick(gdd_file)
  DT_same <- TRUE
}

if (gdd_data == "calculate"){
  
  #Order by date starting with first downloaded date. Up to user to download date range
  #in line with species' biofix. Assumes Jan 01  biofix
  #Sorting is necessary in most recent year with provisional data and different filename structure
  #Search pattern for PRISM daily temperature grids. Load them for processing.
  pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
  files <- list.files(path = prism_path, pattern=pattern, all.files=FALSE, full.names=TRUE)
  #Check that there are enough files for the year
  length(files)
  numlist <- vector()
  for (file in files) {
    num = strsplit(file,split="_")[[1]][5]
    numlist <-c(numlist,num)
  }
  
  sortedlist <- sort(numlist)
  #Initialize Tracking Rasters
  #Read in first raster as a template and intialize tracking rasters
  template <- raster(files[1])
  template <- crop(raster(files[1]), REGION)
  template[!is.na(template)] <- 0
  # dataType(template) <- "INT1U"
  
  
  # if all thresholds same across lifestages, save time by calculated DD now
  if (isTRUE(all.equal(eggLDT, larvaeLDT, pupaeLDT, adultLDT)) &
      isTRUE(all.equal(eggUDT, larvaeUDT, pupaeUDT, adultUDT))){
    
    # load tmin/tmax, crop, calculate gdd for all days
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
    
    tminstack <- crop(tminstack, REGION)
    tmaxstack <- crop(tmaxstack, REGION)
    
    #adapted triangle DD function for RasterBrick
    #TODO: include other DDcalc functions    
    
    LDT <- eggLDT
    UDT <- eggUDT
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
                   filename = "dailygdd",
                   recycle = FALSE,
                   overwrite = TRUE)
    # })
    names(GDD) <- sortedlist
    DT_same <- TRUE
  }
  
} # close gdd_data == calculate


if (region_param == "CONUS"){
  # splits template into list of smaller map rasters to run in parallel
  # only doing this for CONUS, because benefit lost for smaller map sizes
  SplitMap <- SplitRas(template, ppside = 2, TRUE, FALSE)
  # runs_par <- expand.grid(1:nsim, 1:length(SplitMap))
  # names(runs_par) <- c("sim", "map")
}else{
  SplitMap <- list(template)
  # runs_par <- expand.grid(1:nsim, 1:length(SplitMap))
  # names(runs_par) <- c("sim", "map")
}

if (vary_indiv == 1){
  eggDD = eggDD_mu
  larvaeDD = larvaeDD_mu
  pupDD = pupDD_mu
  adultDD = adultDD_mu
  OWadultDD = substages[, 1] # if using SubstageDistrib function
  cdl_b0 <- coefs[1]
  cdl_b1 <- coefs[2]
  # CDL <- CDL_mu
  params <- data.frame(eggDD, larvaeDD, pupDD, adultDD, OWadultDD, cdl_b0, cdl_b1)
}else{ # doesn't make sense to do this with foreach if no variation, though!
  eggDD = eggDD_mu
  larvaeDD = larvaeDD_mu
  pupDD = pupDD_mu
  adultDD = adultDD_mu
  OWadultDD = OWadultDD_mu
  CDL <- CDL_mu
  params <- data.frame(eggDD, larvaeDD, pupDD, adultDD, OWadultDD, CDL)
}

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

system.time({
  outfiles <- foreach(sim = 1:nsim,
                      # outfiles <- foreach(sim = 5, # if some runs don't work, rerun individually
                      .packages= "raster",
                      .export = c("SplitMap", "GDD", 
                                  "newname", "params"),
                      .inorder = FALSE) %:% 
    foreach(map = 1:length(SplitMap),
            # foreach(map = 3, # if some runs don't work, rerun individually
            .packages = "raster",
            .export = c("SplitMap", "GDD", 
                        "newname", "params"),
            .inorder = FALSE) %dopar%{
              # .inorder = FALSE) %do%{

              # creates unique filepath for temp directory to store rasters
              tmppath <- paste0("~/REPO/photovoltinism/rastertmp/", "run", map, sim)
              dir.create(path = tmppath, showWarnings = FALSE)
              #sets temp directory
              rasterOptions(tmpdir=file.path(tmppath)) 
              
              template <- SplitMap[[map]]
              tmpGDD <- crop(GDD, template)
              tmpGDD <- tmpGDD[[start_doy:end_doy]]
              
              # old nested foreach loop
              # foreach(sim = 1:nsim, .packages= "raster") %dopar% {
              
              # Initialize all tracking rasters as zero with the template
              eggDD <- params$eggDD[sim]
              larvaeDD <- params$larvaeDD[sim]
              pupDD <- params$pupDD[sim]
              adultDD <- params$adultDD[sim]
              OWadultDD <- params$OWadultDD[sim]
              # CDL <- params$CDL[sim]
              cdl_b0 <- params$cdl_b0
              cdl_b1 <- params$cdl_b1
              
              DDaccum <- template
              # ddtotal <- template
              Lifestage <- template
              Lifestage <- Lifestage - 1  # Need for OW designated as stage -1
              NumGen <- template
              #Lifestage: [0] = egg, [1] = larvae, [2] = pupae, [3] = adult
              
              #### Init Lifestage tracking
              if (owstage == "OE") {
                LSOW0 <- Lifestage == -1 
              } else if (owstage == "OL") {
                LSOW1 <- Lifestage == -1 
              } else if (owstage == "OP") {
                LSOW2 <- Lifestage == -1 
              } else if (owstage == "OA") {
                LSOW3 <- Lifestage == -1 
              }
              LS0 <- Lifestage == 0
              LS1 <- Lifestage == 1
              LS2 <- Lifestage == 2
              LS3 <- Lifestage == 3
              LS4 <- Lifestage == 4
              
              # loop over days
              sublist <- sortedlist[start_doy:end_doy]
              for (d in sublist) {
                print(d)
                index <- which(sublist == d)
                if (.Platform$OS.type == "windows") flush.console()
                Sys.sleep(1)
                
                # If GDD not pre-calculated
                # Read in day's PRISM raster files
                if (isFALSE(DT_same)){
                  pattern = paste("(PRISM_tmin_)(.*)(",d,")(_bil.bil)$", sep="")
                  temp <- list.files(path = prism_path, pattern=pattern,all.files=FALSE, full.names=TRUE)
                  tmin <- crop(raster(temp), template)
                  pattern = paste("(PRISM_tmax_)(.*)(",d,")(_bil.bil)$", sep="")
                  temp <- list.files(path = prism_path, pattern=pattern,all.files=FALSE, full.names=TRUE)
                  tmax <- crop(raster(temp), template)
                  rm(pattern,temp)
                  
                  # tmean not in GRUB PRISM data, use approximation
                  tmean <- (tmax + tmin) / 2
                }
                
                # photoperiod for this day across raster
                if (model_CDL == 1){
                  doy <- lubridate::yday(lubridate::ymd(d))
                  photo <- RasterPhoto(template, doy, perc_twilight = 25)
                }
                
                #### Loop through Stages; order of stages now read from SPP param file ####
                
                for (i in stgorder) {  # Handle stages in the model
                  ####  MAIN STEPS FOR EGG STAGE ####
                  if (i == "E" | i == "OE") {   # Egg Stage
                    if(isFALSE(DT_same)){
                      if(calctype=="average") { #devel DDs (zero values for temps below LDT)
                        dd0tmp <- AvgDD(tmax,tmin,eggLDT,eggUDT)
                      } else if(calctype=="triangle") {
                        dd0tmp <- TriDD(tmax,tmin,eggLDT,eggUDT)
                      } else { # assume (calctype=="simple") 
                        dd0tmp <- SimpDD(tmean,eggLDT)
                      }
                    }else if(isTRUE(DT_same)){
                      # extract layer from RasterBrick
                      dd0tmp <- tmpGDD[[names(tmpGDD)[index]]]
                    }
                    
                    if (i == "OE") { 
                      dd0 <- dd0tmp * LSOW0 
                    } else if (i == "E") { 
                      dd0 <- dd0tmp[[1]] * LS0 
                    }
                    
                    #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
                    #count cells with dd>0
                    dd.stat <- cellStats(dd0,stat='max',na.rm=TRUE)
                    if (dd.stat > 0) {
                      DDaccum <- DDaccum + dd0
                      #Calculate lifestage progression: for dd accum in correct lifestage, is it >= egg DD threshold?
                      if (i == "OE") {
                        progressOW0 <- (DDaccum * LSOW0) >= OWeggDD
                        Lifestage <- Cond(LSOW0 == 1 & progressOW0 == 1, 1, Lifestage)
                        #Reset the DDaccum cells to zero for cells that progressed to next lifestage
                        DDaccum <- Cond(progressOW0 == 1,(DDaccum - OWeggDD) * LSOW0, DDaccum)
                      } else if (i == "E") {
                        progress0 <- (DDaccum * LS0) >= eggDD
                        Lifestage <- Cond(LS0 == 1 & progress0 == 1, 1, Lifestage)
                        #Reset the DDaccum cells to zero for cells that progressed to next lifestage
                        DDaccum <- Cond(progress0 == 1,(DDaccum - eggDD) * LS0, DDaccum)
                      }
                    }
                    
                    # memory management
                    rm(dd0, dd0tmp)
                    gc()
                    
                    
                    ####  MAIN STEPS FOR LARVAL STAGE ####
                  } else if (i == "L" | i == "OL") {  # Larval Stage
                    #developmental degree days
                    if(isFALSE(DT_same)){
                      if(calctype=="average") {
                        dd1tmp <- AvgDD(tmax,tmin,larvaeLDT,larvaeUDT)
                      } else if(calctype=="triangle") {
                        dd1tmp <- TriDD(tmax,tmin,larvaeLDT,larvaeUDT)
                      } else { # assume (calctype=="simple") 
                        dd1tmp <- SimpDD(tmean,larvaeLDT)
                      }
                      # ddtotal <- ddtotal + dd1tmp #LEN : Accumulate total degree days for the year for larvae
                    }else if(isTRUE(DT_same)){
                      # extract layer from RasterBrick
                      dd1tmp <- tmpGDD[[names(tmpGDD)[index]]]
                    }
                    
                    if (i == "OL") { 
                      dd1 <- dd1tmp * LSOW1
                    }else if (i == "L") { 
                      dd1 <- dd1tmp * LS1
                    }
                    
                    #Accumulate degree days, if dd1 > 0 otherwise exclusion masks get applied to this generation.
                    dd.stat <- cellStats(dd1,stat='max',na.rm=TRUE)
                    if (dd.stat > 0) {
                      DDaccum <- DDaccum + dd1
                    }
                    #Calculate lifestage progression: for dd accum in correct lifestage, is it >= egg DD threshold?
                    if (i == "OL") {
                      progressOW1 <- (DDaccum * LSOW1) >= OWlarvaeDD
                      Lifestage <- Cond(LSOW1 == 1 & progressOW1 == 1, 2, Lifestage)
                      #Reset the DDaccum cells to zero for cells that progressed to next lifestage
                      DDaccum <- Cond(progressOW1 == 1,(DDaccum - OWlarvaeDD) * LSOW1, DDaccum)
                    } else if (i == "L") {
                      progress1 <- (DDaccum * LS1) >= larvaeDD
                      Lifestage <- Cond(LS1 == 1 & progress1 == 1, 2, Lifestage)
                      DDaccum <- Cond(progress1 == 1,(DDaccum - larvaeDD) * LS1, DDaccum)
                    }
                    
                    # memory management
                    rm(dd1, dd1tmp)
                    gc()
                    
                    ####  MAIN STEPS FOR PUPAL STAGE ####
                  } else if (i == "P" | i == "OP") {   # Pupal Stage
                    if(isFALSE(DT_same)){
                      #developmental degree days
                      if(calctype=="average") {
                        dd2tmp <- AvgDD(tmax,tmin,pupaeLDT,pupaeUDT)
                      } else if(calctype=="triangle") {
                        dd2tmp <- TriDD(tmax,tmin,pupaeLDT,pupaeUDT)
                      } else { # assume (calctype=="simple") 
                        dd2tmp <- SimpDD(tmean,pupaeLDT)
                      }
                    }else if(isTRUE(DT_same)){
                      # extract layer from RasterBrick
                      dd2tmp <- tmpGDD[[names(tmpGDD)[index]]]
                    }
                    
                    if (i == "OP") { 
                      dd2 <- dd2tmp * LSOW2
                    } else if (i == "P") { 
                      dd2 <- dd2tmp * LS2
                    }
                    
                    #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
                    dd.stat <- cellStats(dd2,stat='max',na.rm=TRUE)
                    if (dd.stat > 0) {
                      DDaccum <- DDaccum + dd2
                    }
                    #Calculate lifestage progression: for dd accum in correct lifestage, is it >= pupae DD threshold?
                    if (i == "OP") {
                      progressOW2 <- (DDaccum * LSOW2) >= OWpupDD
                      Lifestage <- Cond(LSOW2 == 1 & progressOW2 == 1, 3, Lifestage)
                      #Reset the DDaccum cells to zero for cells that progressed to next lifestage
                      DDaccum <- Cond(progressOW2 == 1,(DDaccum - OWpupDD) * LSOW2, DDaccum)
                    } else if (i == "P") {
                      progress2 <- (DDaccum * LS2) >= pupDD
                      Lifestage <- Cond(LS2 == 1 & progress2 == 1, 3, Lifestage)
                      DDaccum <- Cond(progress2 == 1,(DDaccum - pupDD) * LS2, DDaccum)
                    }
                    
                    # memory management
                    rm(dd2, dd2tmp)
                    gc()
                    
                    ####  MAIN STEPS FOR ADULT STAGE ####
                  } else if (i == "A" | i == "OA") {  # Adult stage, or time to 50% oviposition
                    if(isFALSE(DT_same)){
                      #developmental degree days
                      if(calctype=="average") {
                        dd3tmp <- AvgDD(tmax,tmin,adultLDT,adultUDT)
                      } else if(calctype=="triangle") {
                        dd3tmp <- TriDD(tmax,tmin,adultLDT,adultUDT)
                      } else { # assume (calctype==simple) 
                        dd3tmp <- SimpDD(tmean,adultLDT)
                      }
                    }else if(isTRUE(DT_same)){
                      # extract layer from RasterBrick
                      dd3tmp <- tmpGDD[[names(tmpGDD)[index]]]
                    }
                    
                    if (i == "OA") { 
                      dd3 <- dd3tmp * LSOW3 
                    } else if (i == "A") { 
                      dd3 <- dd3tmp * LS3 
                    }
                    
                    #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
                    dd.stat <- cellStats(dd3,stat='max',na.rm=TRUE)
                    if (dd.stat > 0) {
                      DDaccum <- DDaccum + dd3
                    }
                    
                    #Calculate lifestage progression: for dd accum in correct lifestage, is it >= adult DD threshold?
                    if (i == "OA") {
                      progressOW3 <- (DDaccum * LSOW3) >= OWadultDD
                      DDaccum <- Cond(progressOW3 == 1,(DDaccum - OWadultDD) * LSOW3, DDaccum)
                      progressOW3[is.na(progressOW3)] <- template[is.na(progressOW3)]
                      Lifestage <- Cond(LSOW3 == 1 & progressOW3 == 1, 0, Lifestage)
                      Lifestage[is.na(Lifestage)] <- template[is.na(Lifestage)]
                      #Increment NumGen + 1
                      NumGen <- NumGen + progressOW3
                      #writeRaster(NumGen,paste("NumGen",d,sep=""), format="GTiff",overwrite=TRUE)
                    } else if (i == "A") {
                      progress3 <- (DDaccum * LS3) >= adultDD
                      #Reset the DDaccum cells to zero for cells that progressed to next lifestage
                      DDaccum <- Cond(progress3 == 1,(DDaccum - adultDD) * LS3, DDaccum)
                      #Remove masking effect from progress counter so it doesn't perpetuate through NumGen raster
                      progress3[is.na(progress3)] <- template[is.na(progress3)]
                      Lifestage <- Cond(LS3 == 1 & progress3 == 1,0, Lifestage)
                      Lifestage[is.na(Lifestage)] <- template[is.na(Lifestage)]
                      #Increment NumGen + 1
                      NumGen <- NumGen + progress3
                    }
                    
                    # memory management
                    if (i == "OA"){
                      rm(progressOW3)
                    } else if (i == "A"){
                      rm(progress3)
                    }
                    rm(dd3, dd3tmp)
                    gc()
                    
                    ####  MAIN STEPS FOR END OF EACH DAY ####
                  } else if (i == "F") { # end of the day placeholder 
                    
                    if (model_CDL == 1){
                      # # First try
                      # # diapause decision with no variation
                      # # add one more integer to Lifestage, 4 == diapause
                      # sens_mask <- Cond(Lifestage %in% photo_sens, 1, 0)
                      # diap_mask <- Cond(photo < CDL, sens_mask, 0)
                      # rm(sens_mask)
                      # Lifestage <- Cond(diap_mask == 1, 4, Lifestage)
                      # rm(diap_mask)
                      # LS4 <- Lifestage == 4
                      # gc()
                      
                      # Second Try
                      # add logistic regression variation in CDL response
                      # keep running lifecycle as if all directly developing
                      # but add new raster to track percent of population actually
                      # in diapause following logistic regression
                      sens_mask <- Cond(Lifestage %in% photo_sens, 1, 0)
                      prop_diap <- 1 - exp(cdl_b0 + cdl_b1 * photo) /
                        (1 + exp(cdl_b0 + cdl_b1 * photo))
                      tmpLS4 <- Cond(sens_mask == 1, prop_diap, LS4)
                      # need to account for prop_diap declining up until solstice
                      # only let proportion diapausing increase over season
                      LS4 <- Cond(prop_diap < LS4, LS4, tmpLS4)
                      rm(sens_mask, prop_diap, tmpLS4)
                      gc()
                    }
                    
                    # update lifestage masks
                    if (owstage == "OE") {
                      LSOW <- LSOW0 <- Lifestage == -1
                    } else if (owstage == "OL") {
                      LSOW <- LSOW1 <- Lifestage == -1
                    } else if (owstage == "OP") {
                      LSOW <- LSOW2 <- Lifestage == -1
                    } else if (owstage == "OA") {
                      LSOW <- LSOW3 <- Lifestage == -1
                    }
                    LS0 <- Lifestage == 0
                    LS1 <- Lifestage == 1
                    LS2 <- Lifestage == 2
                    LS3 <- Lifestage == 3
                    
                    # stack each days raster
                    if (!exists("LS0stack")){
                      # hope that this will keep this large files on disk, out of memory load
                      LSOWstack <- stack(LSOW)
                      LS0stack <- stack(LS0)
                      LS1stack <- stack(LS1)
                      LS2stack <- stack(LS2)
                      LS3stack <- stack(LS3)
                      LS4stack <- stack(LS4)
                    } else {
                      LSOWstack <- addLayer(LSOWstack, LSOW)
                      LS0stack <- addLayer(LS0stack, LS0)
                      LS1stack <- addLayer(LS1stack, LS1)
                      LS2stack <- addLayer(LS2stack, LS2)
                      LS3stack <- addLayer(LS3stack, LS3)
                      LS4stack <- addLayer(LS4stack, LS4)
                      
                      # #alternatives
                      # # too slow
                      #  LS0stack <- stack(LS0stack, LS0)
                      # # faster but dangerous, doesn't check extent in each file
                      #  LS0stack@layers[[index]] <- LS0
                    }
                    
                    if (!exists("NumGenstack")){
                      NumGenstack <- stack(NumGen)
                    } else {
                      NumGenstack <- addLayer(NumGenstack, NumGen)
                    }
                    
                  }  # lifestage F (end of day calcs)
                }  # lifestage for loop
              }  # daily loop
              
              
              # if parallel by map chunk
              # each band is a day, first index is map chunk index
              mapcode <- formatC(map, width = 3, format = "d", flag = "0")
              
              NumGenFile <- writeRaster(NumGenstack,
                                        filename = paste(newname, "/NumGen_", mapcode, "_sim", sim, sep = ""),
                                        overwrite = TRUE, datatype = "INT1U")
              LSOWFile <- writeRaster(LSOWstack, filename = paste(newname, "/LSOW_", mapcode, "_sim", sim, sep = ""),
                                      overwrite = TRUE, datatype = "INT1U")
              LS0File <- writeRaster(LS0stack, filename = paste(newname, "/LS0_", mapcode, "_sim", sim, sep = ""),
                                     overwrite = TRUE, datatype = "INT1U")
              LS1File <- writeRaster(LS1stack, filename = paste(newname, "/LS1_", mapcode, "_sim", sim, sep = ""),
                                     overwrite = TRUE, datatype = "INT1U")
              LS2File <- writeRaster(LS2stack, filename = paste(newname, "/LS2_", mapcode, "_sim", sim, sep = ""),
                                     overwrite = TRUE, datatype = "INT1U")
              LS3File <- writeRaster(LS3stack, filename = paste(newname, "/LS3_", mapcode, "_sim", sim, sep = ""),
                                     overwrite = TRUE, datatype = "INT1U")
              LS4File <- writeRaster(LS4stack, filename = paste(newname, "/LS4_", mapcode, "_sim", sim, sep = ""),
                                     overwrite = TRUE, datatype = "FLT8S")
              
              # memory management
              rm(DDaccum, LSOWstack, LS0stack, LS1stack, LS2stack, LS3stack, LS4stack, NumGenstack)
              #removes entire temp directory without affecting other running processes
              unlink(tmppath, recursive = TRUE)
              gc()
              
              # keep track of written output files
              out <- lapply(X = list(NumGenFile,
                                     LSOWFile,
                                     LS0File, LS1File, LS2File,
                                     LS3File, LS4File), filename)
            } # close foreach loop
}) #system.time

if(.Platform$OS.type == "windows"){
  stopCluster(cl)
}


# Run model, but only with photoperiod update for speed
# The way I built the diapause decision, model results
# already track lifecycle through the end of GDD accumulation.
# To update photoperiod response, just need to take model
# results and rerun sens stage x photoperiod cue = diapause.

# DOESN'T WORK WITH WEIGHTED RESULTS
# NEED SEPARATE SUBSTAGE SIMS I THINK

source('CDL_funcs.R')
region_param <- "WEST"
gdd_file <- "dailygdd_2017_WEST.grd"

REGION <- switch(region_param,
                 "CONUS"        = extent(-125.0,-66.5,24.0,50.0),
                 "NORTHWEST"    = extent(-125.1,-103.8,40.6,49.2),
                 "OR"           = extent(-124.7294, -116.2949, 41.7150, 46.4612),
                 "TEST"         = extent(-124, -122.5, 44, 45),
                 "WEST"         = extent(-125.14, -109, 37, 49.1))


GDD <- brick(gdd_file)

template <- GDD[[1]]
template[!is.na(template)] <- 0
template <- crop(template, REGION)

# input sensitive stage
# values are 0-1 proportion of population in stage
Lifestage <- brick("GCA_WEST_SCDL_2017/LS3_001_weighted.grd")
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

################################  
