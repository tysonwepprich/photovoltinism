# Take Gerricke's FCM code and convert to Galerucella parameters
# Test out photoperiod and substages with constrained extexts

# SLOW, 
# ways to speed up, reduce memory load
# 1. datatype for raster to integer instead of floating point, 
# would need to multiply PRISM data to remove decimals
# 2. multithread by year

# NW, 10 days, 7 substages
# 36: 332
# 16: 272
# 1: 289
# 1 (DT_same): 169
# 1 (DT_same, sims parallel): 34 ---- with 365 days, only 710!

# TODO
# 1. Write file names with 3 digits so sorted the right way
# formatC(SiteID, width = 3, format = "d", flag = "0")
# 2. Better strategy for parallel, now with map chunks there's no benefit 
# of running 36 cores vs 16 in some tests. 
# 3. Output has strange shoreline edges with no results. Not NA, but
# also does not advance lifestages. Issue with PRISM or raster?
# 4. Multiple years for GDD to match with 2015 paper, possibly store
# all Tmax/Tmin PRISM data in one RasterBrick to speed cropping
# 5. Photoperiod threshold flexibility/logistic regression probability
# 6. Pest Event Maps
# 7. Tracking partial generations due to diapause for each substage, 
# could be done with one raster per substage with % population removed?
# 8. Diapause stage leading to overwintering stage (using same raster?)

library(sp)
library(rgdal)
library(raster)

rasterOptions(overwrite = FALSE, 
              chunksize = 1e+08,
              maxmemory = 1e+09)

# for parallel simulations with control over seed for reproducibility
library(doRNG)
library(foreach) # for parallelized loops
library(doMC)    # parallel backend for foreach


# setwd("F:/PRISM/2014") # remove this
prism_path <- "/data/PRISM/2014"
# prism_path <- "/data/PRISM/"

# # TODO loop with multiple years
# # OR take years and extract "average" GDD accumulation
# prism_path <- "/data/PRISM/"
# years <- c(2007:2013)
# for (yr in years){
#   path <- paste(prism_path, yr, sep = "/")
#   
# }


source('CDL_funcs.R')
#Daily PRISM Statistics (degree days and thresholds) for pest development models


####Pest Specific, Multiple Life Stage Phenology Model Parameters:
#LDT = lower development threshold, temp at which growth = 0 (using PRISM tmean)
#DD = degree days, number of cumulative heat units to complete that lifestage
start_doy  <- 1
end_doy    <- 365
stgorder   <- c("OA","E","L","P","A","F")
photo_sens <- 3 #c(-1, 3) # integer life stages for now
CDL_mu        <- 14.25
model_CDL  <- 1 # if 1, model photoperiod decision
owstage    <- "OA"
OWadultDD_mu  <- 100 #108 # text OW stage dev 39 DD "post diapause"
calctype   <-"triangle"
eggLDT     <- 10
eggUDT     <- 37.8  #
larvaeLDT  <- 10
larvaeUDT  <- 37.8  #upper dev. threshold-need to verify
pupaeLDT   <- 10
pupaeUDT   <- 37.8
adultLDT   <- 10 #for oviposition
adultUDT   <- 37.8
eggDD_mu = 100 #93.3
larvaeDD_mu = 140 #136.4 # 46.9 + 45.8 + 43.7 instars
pupDD_mu = 140 #137.7 
adultDD_mu = 130 #125.9 #time to oviposition
region_param <- "NORTHWEST"

# introducing individual variation, tracked with simulations
nday <- length(start_doy:end_doy)
nsim <- 7 # for now, use square number
DD_sd <- 10
CDL_sd <- .1
vary_indiv <- 1 # turn on indiv. variation
# Gaussian quadrature to assign number/weights for substages
substages <- gauss.hermite(11, iterlim = 50)[3:9, ] # remove ends with near zero weights

# try to run sims in parallel rather than splitting map
ncores <- nsim
registerDoMC(cores = ncores)



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

# preload PRISM data into RasterBrick

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


#Order by date starting with first downloaded date. Up to user to download date range
#in line with species' biofix. Assumes Jan 01  biofix
#Sorting is necessary in most recent year with provisional data and different filename structure
sortedlist <- sort(numlist)
filelist <- vector()
#Initialize Tracking Rasters
#Read in first raster as a template and intialize tracking rasters
template <- raster(files[1])
REGION <- switch(region_param,
                 "CONUS"        = extent(-125.0,-66.5,24.0,50.0),
                 "NORTHWEST"    = extent(-125.1,-103.8,40.6,49.2),
                 "OR"           = extent(-124.7294, -116.2949, 41.7150, 46.4612),
                 "TEST"         = extent(-124, -122.5, 44, 45))
template <- crop(raster(files[1]), REGION)
template[!is.na(template)] <- 0

# crop tmin/tmax stacks
# slower to brick, then crop surprisingly
# system.time({
#   tminbrick <- brick(tminstack)
#   tminbrick <- crop(tminbrick, REGION)
# })


tminstack <- crop(tminstack, REGION)
tmaxstack <- crop(tmaxstack, REGION)

# if all thresholds same across lifestages, save time by calculated DD now
if (isTRUE(all.equal(eggLDT, larvaeLDT, pupaeLDT, adultLDT)) &
    isTRUE(all.equal(eggUDT, larvaeUDT, pupaeUDT, adultUDT))){
      
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


# new try: split raster into smaller chunks to run in parallel
# each simulation just accumulates, so only one raster/lifestage saved at once per chunk

# # splits template into list of smaller map rasters
# SplitMap <- SplitRas(template, ppside = floor(sqrt(ncores)), TRUE, FALSE)

# draw parameters outside of parallel loops so all maps use same for the sims

set.seed(1234)
# Draw new params for each simulation run
if (vary_indiv == 1){
  # eggDD = rnorm(nsim, mean = eggDD_mu, sd = DDsd)
  # larvaeDD = rnorm(nsim, mean = larvaeDD_mu, sd = DDsd)
  # pupDD = rnorm(nsim, mean = pupDD_mu, sd = DDsd)
  # adultDD = rnorm(nsim, mean = adultDD_mu, sd = DDsd)
  # OWadultDD = rnorm(nsim, mean = OWadultDD_mu, sd = DDsd)
  # CDL <- rnorm(nsim, mean = CDL_mu, sd = CDLsd)
  eggDD = eggDD_mu
  larvaeDD = larvaeDD_mu
  pupDD = pupDD_mu
  adultDD = adultDD_mu
  OWadultDD = substages[, 1] * DD_sd + OWadultDD_mu
  # CDL <- substages[, 1] * CDL_sd + CDL_mu
  CDL <- CDL_mu
}else{ # doesn't make sense to do this with foreach if no variation, though!
  eggDD = eggDD_mu
  larvaeDD = larvaeDD_mu
  pupDD = pupDD_mu
  adultDD = adultDD_mu
  OWadultDD = OWadultDD_mu
  CDL <- CDL_mu
}
params <- data.frame(eggDD, larvaeDD, pupDD, adultDD, OWadultDD, CDL)
# params <- data.frame(eggDD, larvaeDD, pupDD, adultDD, OWadultDD)

# idea was to cross CDL variation with GDD variation, but still 
# runs into issue of correlations between life stage speed.
# also, will still be too many simulations if we want speed
# params <- do.call("rbind", replicate(n = length(CDL), params, simplify = FALSE))
# params$CDL <- rep(CDL, each = length(eggDD))

newname <- paste("run", gsub("-", "", 
                             gsub(" ", "_", 
                                  gsub(":", "", Sys.time(), fixed = TRUE),
                                  fixed = TRUE),
                             fixed = TRUE), sep = "")
dir.create(newname)

system.time({
# foreach(map = 1:length(SplitMap), .packages= "raster") %dopar% {
#   print(map)
#   template <- SplitMap[[map]]
  foreach(sim = 1:nsim, .packages= "raster") %dopar% {
  
  #Initialize all tracking rasters as zero with the template
  
  # loop over simulations
  
  # for (sim in 1:nsim){
    print(sim)
    eggDD <- params$eggDD[sim]
    larvaeDD <- params$larvaeDD[sim]
    pupDD <- params$pupDD[sim]
    adultDD <- params$adultDD[sim]
    OWadultDD <- params$OWadultDD[sim]
    CDL <- params$CDL[sim]
    
    DDaccum <- template
    ddtotal <- template
    Lifestage <- template
    Lifestage <- Lifestage - 1  # Need for OW designated as stage -1
    NumGen <- template
    #rm(template)
    #Lifestage: [0] = egg, [1] = larvae, [2] = pupae, [3] = adult
    
    #### Init OWSTAGE maps #### 
    if (owstage == "OE") {
      OWeggP <- template # proportion in OW egg stage
    } else if (owstage == "OL") {
      OWlarvP <- template # proportion in OW larval stage
    } else if (owstage == "OP") {
      OWpupP <- template # proportion in OW pup stage
    } else if (owstage == "OA") {
      OWadultP <- template # proportion in OW adult stage
    } 
    
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
    
    
    
    # loop over days
    sublist <- sortedlist[start_doy:end_doy]
    for (d in sublist) {
      print(d)
      index <- which(sublist == d)
      if (.Platform$OS.type == "windows") flush.console()
      Sys.sleep(1)
      
      if (isFALSE(DT_same)){
      #Read in that day's PRISM raster files
      # pattern = paste("(PRISM_tmean_)(.*)(",d,")(_bil.bil)$", sep="")
      # temp <- list.files(pattern=pattern,all.files=FALSE, full.names=TRUE)
      # tmean <- raster(temp)
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
        photo <- RasterPhoto(template, doy)
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
            dd0tmp <- GDD[[paste("X", d, sep = "")]]
          }
          
          if (i == "OE") { dd0 <- dd0tmp * LSOW0 
          } else if (i == "E") { 
            dd0 <- dd0tmp * LS0 
          }
          
          
          #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
          #count cells with dd>0
          dd.stat <- cellStats(dd0,stat='max',na.rm=TRUE)
          if (dd.stat > 0) {
            DDaccum <- DDaccum + dd0
            # if (pems) {
            #   # set DOY in PEM if curr gen, not set yet, DDs > eventDDs, and Stage is eggs, else keep same value
            #   if (PEMnumgens > 0) {
            #     PEMe1 <- Cond(NumGen==1 & PEMe1 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe1) 
            #   }
            #   if (PEMnumgens > 1) {
            #     PEMe2 <- Cond(NumGen==2 & PEMe2 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe2) 
            #   }
            #   if (PEMnumgens > 2) {
            #     PEMe3 <- Cond(NumGen==3 & PEMe3 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe3) 
            #   }
            #   if (PEMnumgens > 3) {
            #     PEMe4 <- Cond(NumGen==4 & PEMe4 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe4) 
            #   }
            # }
            #Calculate lifestage progression: for dd accum in correct lifestage, is it >= egg DD threshold?
            if (i == "OE") {
              #cat("date and doy, i =",d,doy,i,"\n")
              progressOW0 <- (DDaccum * LSOW0) >= OWeggDD
              #writeRaster(progressOW0,paste("ProgressOW0_",d,sep=""), format="GTiff",overwrite=TRUE)
              Lifestage <- Cond(LSOW0 == 1 & progressOW0 == 1, 1, Lifestage)
              #Reset the DDaccum cells to zero for cells that progressed to next lifestage
              DDaccum <- Cond(progressOW0 == 1,(DDaccum - OWeggDD) * LSOW0, DDaccum)
            } else if (i == "E") {
              #cat("date and doy, i =",d,doy,i,"\n")
              progress0 <- (DDaccum * LS0) >= eggDD
              #writeRaster(progress0,paste("Progress0_",d,sep=""), format="GTiff",overwrite=TRUE)
              Lifestage <- Cond(LS0 == 1 & progress0 == 1, 1, Lifestage)
              #Reset the DDaccum cells to zero for cells that progressed to next lifestage
              DDaccum <- Cond(progress0 == 1,(DDaccum - eggDD) * LS0, DDaccum)
            }
          }
          
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
            dd1tmp <- GDD[[paste("X", d, sep = "")]]
          }
          
          #else { # just use lifestage mask
          if (i == "OL") { dd1 <- dd1tmp * LSOW1 }
          else if (i == "L") { dd1 <- dd1tmp * LS1 }
          #}
          
          #Accumulate degree days, if dd1 > 0 otherwise exclusion masks get applied to this generation.
          dd.stat <- cellStats(dd1,stat='max',na.rm=TRUE)
          if (dd.stat > 0) {
            DDaccum <- DDaccum + dd1
            # if (pems & larvaeEventDD) {
            #   # set DOY in PEM if curr. gen, not set yet, DDs > eventDDs, Stage is larvae, else keep same value
            #   if (PEMnumgens > 0) {
            #     PEMl1 <- Cond(NumGen==1 & PEMl1 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl1) 
            #   }
            #   if (PEMnumgens > 1) {
            #     PEMl2 <- Cond(NumGen==2 & PEMl2 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl2) 
            #   }
            #   if (PEMnumgens > 2) {
            #     PEMl3 <- Cond(NumGen==3 & PEMl3 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl3) 
            #   }
            #   if (PEMnumgens > 3) {
            #     PEMl4 <- Cond(NumGen==4 & PEMl4 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl4) 
            #   }
          }
          #Calculate lifestage progression: for dd accum in correct lifestage, is it >= egg DD threshold?
          if (i == "OL") {
            progressOW1 <- (DDaccum * LSOW1) >= OWlarvaeDD
            #writeRaster(progressOW1,paste("ProgressOW1_",d,sep=""), format="GTiff",overwrite=TRUE)
            Lifestage <- Cond(LSOW1 == 1 & progressOW1 == 1, 2, Lifestage)
            #Reset the DDaccum cells to zero for cells that progressed to next lifestage
            DDaccum <- Cond(progressOW1 == 1,(DDaccum - OWlarvaeDD) * LSOW1, DDaccum)
          } else if (i == "L") {
            progress1 <- (DDaccum * LS1) >= larvaeDD
            #writeRaster(progress1,paste("Progress1_",d,sep=""), format="GTiff",overwrite=TRUE)
            Lifestage <- Cond(LS1 == 1 & progress1 == 1, 2, Lifestage)
            DDaccum <- Cond(progress1 == 1,(DDaccum - larvaeDD) * LS1, DDaccum)
          }
          
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
            dd2tmp <- GDD[[paste("X", d, sep = "")]]
          }
          
          if (i == "OP") { dd2 <- dd2tmp * LSOW2
          } else if (i == "P") { dd2 <- dd2tmp * LS2
          }
          
          #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
          dd.stat <- cellStats(dd2,stat='max',na.rm=TRUE)
          if (dd.stat > 0) {
            DDaccum <- DDaccum + dd2
            # if (pems & pupaeEventDD) {
            #   # set DOY in PEM if curr. gen, not set yet, DDs > eventDDs, Stage is pupa, else keep same value
            #   if (PEMnumgens > 0) {
            #     PEMp1 <- Cond(NumGen==1 & PEMp1 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp1) 
            #   }
            #   if (PEMnumgens > 1) {
            #     PEMp2 <- Cond(NumGen==2 & PEMp2 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp2) 
            #   }
            #   if (PEMnumgens > 2) {
            #     PEMp3 <- Cond(NumGen==3 & PEMp3 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp3) 
            #   }
            #   if (PEMnumgens > 3) {
            #     PEMp4 <- Cond(NumGen==4 & PEMp4 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp4) 
            #   }
          }
          #Calculate lifestage progression: for dd accum in correct lifestage, is it >= pupae DD threshold?
          if (i == "OP") {
            progressOW2 <- (DDaccum * LSOW2) >= OWpupDD
            #writeRaster(progressOW2,paste("ProgressOW2_",d,sep=""), format="GTiff",overwrite=TRUE)
            Lifestage <- Cond(LSOW2 == 1 & progressOW2 == 1, 3, Lifestage)
            #Reset the DDaccum cells to zero for cells that progressed to next lifestage
            DDaccum <- Cond(progressOW2 == 1,(DDaccum - OWpupDD) * LSOW2, DDaccum)
          } else if (i == "P") {
            progress2 <- (DDaccum * LS2) >= pupDD
            #writeRaster(progress2,paste("Progress2_",d,sep=""), format="GTiff",overwrite=TRUE)
            Lifestage <- Cond(LS2 == 1 & progress2 == 1, 3, Lifestage)
            DDaccum <- Cond(progress2 == 1,(DDaccum - pupDD) * LS2, DDaccum)
          }
          
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
            dd3tmp <- GDD[[paste("X", d, sep = "")]]
          }
          
          #else { # just use lifestage mask
          if (i == "OA") { dd3 <- dd3tmp * LSOW3 
          } else if (i == "A") { dd3 <- dd3tmp * LS3 }
          #}
          #  if (i == "OA") { dd3 <- dd3tmp * adultmin * LSOW3 }
          #  else if (i == "A") { dd3 <- dd3tmp * adultmin * LS3 }
          # }
          # else { # just use lifestage mask
          #  if (i == "OA") { dd3 <- dd3tmp * LSOW3 }
          #  else if (i == "A") { dd3 <- dd3tmp * LS3 }
          # }
          
          #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
          dd.stat <- cellStats(dd3,stat='max',na.rm=TRUE)
          if (dd.stat > 0) {
            DDaccum <- DDaccum + dd3
            # if(pems & adultEventDD) { # record DOYs in PEMS for adult events
            #   if (PEMnumgens > 0) {
            #     PEMa1 <- Cond(NumGen==1 & PEMa1 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa1) 
            #   }
            #   if (PEMnumgens > 1) {
            #     PEMa2 <- Cond(NumGen==2 & PEMa2 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa2) 
            #   }
            #   if (PEMnumgens > 2) {
            #     PEMa3 <- Cond(NumGen==3 & PEMa3 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa3) 
            #   }
            #   if (PEMnumgens > 3) {
            #     PEMa4 <- Cond(NumGen==4 & PEMa4 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa4) 
            #   }
          }
          
          #Calculate lifestage progression: for dd accum in correct lifestage, is it >= adult DD threshold?
          if (i == "OA") {
            progressOW3 <- (DDaccum * LSOW3) >= OWadultDD
            #writeRaster(progressOW3,paste("ProgressOW3_",d,sep=""), format="GTiff",overwrite=TRUE)
            DDaccum <- Cond(progressOW3 == 1,(DDaccum - OWadultDD) * LSOW3, DDaccum)
            progressOW3[is.na(progressOW3)] <- template[is.na(progressOW3)]
            Lifestage <- Cond(LSOW3 == 1 & progressOW3 == 1, 0, Lifestage)
            Lifestage[is.na(Lifestage)] <- template[is.na(Lifestage)]
            #Increment NumGen + 1
            NumGen <- NumGen + progressOW3
            #writeRaster(NumGen,paste("NumGen",d,sep=""), format="GTiff",overwrite=TRUE)
          } else if (i == "A") {
            progress3 <- (DDaccum * LS3) >= adultDD
            #writeRaster(progress3,paste("Progress3_",d,sep=""), format="GTiff",overwrite=TRUE)
            #Reset the DDaccum cells to zero for cells that progressed to next lifestage
            DDaccum <- Cond(progress3 == 1,(DDaccum - adultDD) * LS3, DDaccum)
            #Remove masking effect from progress counter so it doesn't perpetuate through NumGen raster
            progress3[is.na(progress3)] <- template[is.na(progress3)]
            Lifestage <- Cond(LS3 == 1 & progress3 == 1,0, Lifestage)
            Lifestage[is.na(Lifestage)] <- template[is.na(Lifestage)]
            #Increment NumGen + 1
            NumGen <- NumGen + progress3
            #writeRaster(NumGen,paste("NumGen",d,sep=""), format="GTiff",overwrite=TRUE)
          }
          
          
          ####  MAIN STEPS FOR END OF EACH DAY ####
        } else if (i == "F") { # end of the day placeholder 
          
          if (model_CDL == 1){
            # diapause decision
            # add one more integer to Lifestage, 4 == diapause
            sens_mask <- Cond(Lifestage %in% photo_sens, 1, 0)
            diap_mask <- Cond(photo < CDL, sens_mask, 0)
            Lifestage <- Cond(diap_mask == 1, 4, Lifestage)
            LS4 <- Lifestage == 4
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
            LS0stack <- stack(LS0)
            LS1stack <- stack(LS1)
            LS2stack <- stack(LS2)
            LS3stack <- stack(LS3)
            LS4stack <- stack(LS4)
            
          } else {
            # too slow
           # test <- system.time({
           #  LS0stack <- stack(LS0stack, LS0)
           #  LS1stack <- stack(LS1stack, LS1)
           #  LS2stack <- stack(LS2stack, LS2)
           #  LS3stack <- stack(LS3stack, LS3)
           #  LS4stack <- stack(LS4stack, LS4) 
           #  })
           # print(test)
            # faster but dangerous, doesn't check extent in each file
             LS0stack@layers[[index]] <- LS0
             LS1stack@layers[[index]] <- LS1
             LS2stack@layers[[index]] <- LS2
             LS3stack@layers[[index]] <- LS3
             LS4stack@layers[[index]] <- LS4
          }
          
          if (!exists("NumGenstack")){
            NumGenstack <- stack(NumGen)
          } else {
            # NumGenstack <- stack(NumGenstack, NumGen)
            NumGenstack@layers[[index]] <- NumGen
          }
          # # new Nov 2015 tally stages plus gens
          # StageCount <- Lifestage + (NumGen * 4)
          
        }  # lifestage F or 5 (end of day calcs)
      }  # lifestage for loop
    }  # daily loop
    
    # 
    # # this accumulates simulation values in one stack, first attempt at 
    # # variation in physiological params
    # if (sim == 1){
    #   AggNumGen <- NumGenstack
    #   AggLS0 <- LS0stack
    #   AggLS1 <- LS1stack
    #   AggLS2 <- LS2stack
    #   AggLS3 <- LS3stack
    #   AggLS4 <- LS4stack
    # }else{
    #   AggNumGen <- overlay(AggNumGen, NumGenstack, fun=function(x,y) x + y)
    #   AggLS0 <- overlay(AggLS0, LS0stack, fun=function(x,y) x + y)
    #   AggLS1 <- overlay(AggLS1, LS1stack, fun=function(x,y) x + y)
    #   AggLS2 <- overlay(AggLS2, LS2stack, fun=function(x,y) x + y)
    #   AggLS3 <- overlay(AggLS3, LS3stack, fun=function(x,y) x + y)
    #   AggLS4 <- overlay(AggLS4, LS4stack, fun=function(x,y) x + y)
    # }
    
    #saving files by layer makes a lot of files!
    
    # # if parallel by map chunk
    # #each band is a day, first index is map chunk index
    # mapcode <- formatC(map, width = 3, format = "d", flag = "0")
    # 
    # NumGenFile <- writeRaster(NumGenstack, 
    #                           filename = paste(newname, "/NumGen_", mapcode, "_sim", sim, sep = ""),
    #                           overwrite = TRUE)
    # LS0File <- writeRaster(LS0stack, filename = paste(newname, "/LS0_", mapcode, "_sim", sim, sep = ""),
    #                        overwrite = TRUE)
    # LS1File <- writeRaster(LS1stack, filename = paste(newname, "/LS1_", mapcode, "_sim", sim, sep = ""),
    #                        overwrite = TRUE)
    # LS2File <- writeRaster(LS2stack, filename = paste(newname, "/LS2_", mapcode, "_sim", sim, sep = ""),
    #                        overwrite = TRUE)
    # LS3File <- writeRaster(LS3stack, filename = paste(newname, "/LS3_", mapcode, "_sim", sim, sep = ""),
    #                        overwrite = TRUE)
    # LS4File <- writeRaster(LS4stack, filename = paste(newname, "/LS4_", mapcode, "_sim", sim, sep = ""),
    #                        overwrite = TRUE)
    
    
    NumGenFile <- writeRaster(NumGenstack, 
                              filename = paste(newname, "/NumGen_", "sim", sim, sep = ""),
                              overwrite = TRUE)
    LS0File <- writeRaster(LS0stack, filename = paste(newname, "/LS0_", "sim", sim, sep = ""),
                           overwrite = TRUE)
    LS1File <- writeRaster(LS1stack, filename = paste(newname, "/LS1_", "sim", sim, sep = ""),
                           overwrite = TRUE)
    LS2File <- writeRaster(LS2stack, filename = paste(newname, "/LS2_", "sim", sim, sep = ""),
                           overwrite = TRUE)
    LS3File <- writeRaster(LS3stack, filename = paste(newname, "/LS3_", "sim", sim, sep = ""),
                           overwrite = TRUE)
    LS4File <- writeRaster(LS4stack, filename = paste(newname, "/LS4_", "sim", sim, sep = ""),
                           overwrite = TRUE)
    
    rm(NumGenstack, LS0stack, LS1stack, LS2stack, LS3stack, LS4stack)
    
  # } #sim loop
  
  # #saving files by layer makes a lot of files!
  # #each band is a day, first index is map chunk index
  # mapcode <- formatC(map, width = 3, format = "d", flag = "0")
  # 
  # NumGenFile <- writeRaster(AggNumGen, filename = paste(newname, "/NumGen_", mapcode, sep = ""),
  #                            overwrite = TRUE)
  # LS0File <- writeRaster(AggLS0, filename = paste(newname, "/LS0_", mapcode, sep = ""),
  #                         overwrite = TRUE)
  # LS1File <- writeRaster(AggLS1, filename = paste(newname, "/LS1_", mapcode, sep = ""),
  #                         overwrite = TRUE)
  # LS2File <- writeRaster(AggLS2, filename = paste(newname, "/LS2_", mapcode, sep = ""),
  #                         overwrite = TRUE)
  # LS3File <- writeRaster(AggLS3, filename = paste(newname, "/LS3_", mapcode, sep = ""),
  #                         overwrite = TRUE)
  # LS4File <- writeRaster(AggLS4, filename = paste(newname, "/LS4_", mapcode, sep = ""),
  #                         overwrite = TRUE)
  # 
  # 
} #foreach map chunk loop
})
  
  
  
# TODO: rewrite for no map chunking
  # just skip below, files already saved appropriately
#   
# # then have to retrieve files and mosaic back together
# returnwd <- getwd()
# setwd(newname)
# 
# f <-list.files()
# rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]
# 
# # # for all sims summed together
# # ls <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,1])
# # maps <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,2])
# # maps <- gsub(pattern = ".grd", replacement = "", x = maps)
# 
# # for each sim with unique information to save
# ls <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,1])
# maps <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,2])
# sims <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,3])
# sims <- gsub(pattern = ".grd", replacement = "", x = sims)
# 
# # # for all sims summed together
# # for (i in ls){
# #   fs <- rasfiles[grep(pattern = i, x = rasfiles, fixed = TRUE)]
# #   bricklist <- list()
# #   for (m in 1:length(fs)){
# #     bricklist[[m]] <- brick(fs[m])
# #   }
# #   
# #   bricklist$filename <- paste(i, "all", sep = "_")
# #   bricklist$overwrite <- TRUE
# #   test <- do.call(merge, bricklist)
# # }
# 
# 
# # for each sim with unique information to save
# for (i in ls){
#   fs <- rasfiles[grep(pattern = i, x = rasfiles, fixed = TRUE)]
#   for (j in sims){
#     fs2 <- fs[grep(pattern = j, x = fs, fixed = TRUE)]
#     
#     bricklist <- list()
#     for (m in 1:length(fs2)){
#       bricklist[[m]] <- brick(fs2[m])
#     }
#     
#     bricklist$filename <- paste(i, j, "all", sep = "_")
#     bricklist$overwrite <- TRUE
#     test <- do.call(merge, bricklist)
#   }
# }
# 
# cleanup <- list.files()
# cleanup <- cleanup[-grep("all", x = cleanup)]
# lapply(cleanup, FUN = file.remove) # CAREFUL HERE!
# 
# setwd(returnwd)

# })


# 56 hours for 100 sims of CONUS with 16 cores


saveRDS(params, file = paste(newname, "/params.rds", sep = ""))

# map of Number of Generations
# decimals indicate some simulations reached next integer
res <- brick(paste(newname, "/", "NumGen_sim7.grd", sep = ""))
# res <- res[[365]]/100
plot(res[[seq(10, 360, 50)]])

# map of diapause
res <- brick(paste(newname, "/", "LS1_sim1_all.grd", sep = ""))
plot(res[[seq(100, 250, 30)]])

# map of lifestages
res <- brick(paste(newname, "/", "LS4_sim1_all.grd", sep = ""))
plot(res[[seq(10, 100, 10)]])


#TODO pretty raster facets with common scale, use ggplot or rastervis
#include zoom of NW

library(ggplot2)
library(viridis)

df = as.data.frame(res, xy=TRUE)
names(df)[3] <- "Voltinism"
test <- ggplot(df, aes(x, y, fill = Voltinism)) +
  geom_raster() +
  scale_fill_viridis(na.value = "white") + 
  theme_bw() +
  ggtitle("Voltinism at end of year")
test
ggsave(filename = "voltinism.png", plot = test, device = "png", 
       width = 6, height = 4)





# TODO plotting substages considering their weights

returnwd <- getwd()
setwd(newname)

f <-list.files()
rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]

# for each sim with unique information to save
ls <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,1])
# maps <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,2])
# sims <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,3])
# sims <- gsub(pattern = ".grd", replacement = "", x = sims)

# this loop takes a lot of time
for (i in ls){
  fs <- sort(rasfiles[grep(pattern = i, x = rasfiles, fixed = TRUE)])
  blank <- brick(fs[1]) * 0
  for (j in 1:length(fs)){
    ras_weighted <- brick(fs[j]) * substages[j, 2] # weights for each substage size
    blank <- overlay(blank, ras_weighted, fun=function(x,y) x + y)
  }
  outras <- writeRaster(blank, filename = paste(i, "weighted", sep = "_"),
                         overwrite = TRUE)
}



f <-list.files()
rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]
ras <- brick(rasfiles[47])
plot(ras[[seq(10, 310, 50)]])
plot(ras[[365]])

setwd(returnwd)


library(ggplot2)
library(viridis)

df = as.data.frame(ras[[80]], xy=TRUE)
names(df)[3] <- "Prop"
test <- ggplot(df, aes(x, y, fill = Prop)) +
  geom_raster() +
  scale_fill_viridis(na.value = "white") + 
  theme_bw() +
  ggtitle("Voltinism at end of year")
test


# loop to grab all lifestage results on a certain day to plot together
# save in one giant data.frame to use for plots
returnwd <- getwd()
setwd(newname)

f <-list.files()
# days <- seq(35, 365, 30)
days <- seq(10, 50, 10)
rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]
dflist <- list()
for (i in rasfiles){
  res <- brick(i)
  ls <- unique(stringr::str_split_fixed(i, pattern = "_", 2)[,1])
  for (j in days){
    df <- as.data.frame(res[[j]], xy=TRUE)
    names(df)[3] <- "Percent_of_simulations"
    df$doy <- j
    df$lifestage <- ls
    dflist[[length(dflist)+1]] <- df
  }
}
resdf <- dplyr::bind_rows(dflist)
names(resdf)[3] <- "Proportion"
setwd(returnwd)

nsims <- 100
library(ggplot2)
library(viridis)
library(gridExtra)
library(grid)
library(dplyr)
pltlist <- list()
lifestages <- unique(resdf$lifestage)
ls_labels <- c("egg", "larva", "pupa", "adult", "diapause", "voltinism")
for (d in days){
  for (p in 1:length(lifestages)){
    pltdf <- resdf %>% 
      filter(lifestage == lifestages[p], 
             doy == d)
    if (lifestages[p] == "NumGen"){
      pltdf <- pltdf %>% 
        mutate(Voltinism = Proportion / nsims)
      pltlist[[length(pltlist)+1]] <- ggplot(pltdf, aes(x, y, fill = Voltinism)) +
        geom_raster() +
        scale_fill_viridis(na.value = "white") + 
        theme_bw() +
        ggtitle(ls_labels[p])
    }else{
      pltdf <- pltdf %>% 
        mutate(Proportion = Proportion / nsims)
      pltlist[[length(pltlist)+1]] <- ggplot(pltdf, aes(x, y, fill = Proportion)) +
        geom_raster() +
        scale_fill_viridis(na.value = "white") + 
        theme_bw() +
        ggtitle(ls_labels[p])
    }  
  }
}

for (i in 1:length(days)){
  index <- (i*length(lifestages) - (length(lifestages) - 1)):(i*length(lifestages))
  plt <- grid.arrange(grobs = pltlist[index],
                      ncol = 2,
                      top = paste("DOY", days[i], sep = " "))
  ggsave(paste("DOY", days [i], ".png", sep = ""),
         plot = plt, device = "png", width = 12, height = 12, units = "in")
}





