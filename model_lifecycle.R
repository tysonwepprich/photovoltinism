# Condensed DDRP model
# Universal code for 3 biocontrol species
# Any number of lifestages allowed

# 1. Setup -------
# packages, options, functions loaded
# TODO: find ways to take out dplyr, purrr, mixsmsn functions?
pkgs <- c("sp", "rgdal", "raster", "lubridate", "mixsmsn", "dplyr",
          "stringr", "purrr", "prism", "foreach", "doParallel", "daymetr", "ncdf4")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them

# adjust raster options for your computer, likely this is fine
rasterOptions(overwrite = FALSE, 
              chunksize = 1e+07,
              maxmemory = 1e+08)

# load collection of functions for this model
source('CDL_funcs.R')
source('species_params.R')


# 2. User input -----
# directory with daily tmax and tmin raster files
# PRISM: need to download CONUS which takes a lot of space
# Daymet: can download what you need with bounding box
weather_path <- "prismDL/2017"
download_daily_weather <- 0 # 1 if you need to download PRISM/Daymet data first (20 minutes)
weather_data_source <-"prism" #  "daymet" # or 'prism'


# directory to hold temporary raster results
myrastertmp <- "~/REPO/photovoltinism/rastertmp/"
# run simulations with parallel processing
runparallel <- 1 # 1 for yes, 0 for no

# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model scope
yr           <- 2017
start_doy    <- 1
end_doy      <- 365
region_param <- "TEST" # TEST/WEST/EAST/CONUS/SOUTHWEST/NORTHWEST
species      <- "GCA" # GCA/APHA/DCA
biotype      <- "N" # TODO: add options for each species, N or S for APHA and GCA

# introducing individual variation, tracked with simulation for each substage
# assign to 1 to match previous model versions
nsim <- 7 # number of substages/cohorts to approximate emergence distribution

# photoperiod decision inclusion
# 2 for logistic, 1 for single value CDL, 0 for none
model_CDL  <- 2 


# Weather download ----
# using ropensci 'prism' package to access webservice, or 'daymetr' package
# downloads entire CONUS, so files are large
if (download_daily_weather == 1){
  if (weather_data_source == "prism"){
  startdate <- as.Date(start_doy - 1, origin = paste0(yr, "-01-01"))
  enddate <- as.Date(end_doy - 1, origin = paste0(yr, "-01-01"))
  
  # need to set download directory
  # options(prism.path = paste("prismDL", yr, sep = "/"))
  options(prism.path = prism_path)
  
  get_prism_dailys("tmin", minDate = startdate, 
                   maxDate = enddate, keepZip = FALSE)
  get_prism_dailys("tmax", minDate = startdate, 
                   maxDate = enddate, keepZip = FALSE)
  }
  if (weather_data_source == "daymet"){
    bbox <- assign_extent(region_param)
    # download ncss doesn't work with windows because no OpenDap functioning
    download_daymet_ncss(location = c(bbox@xmin, bbox@ymax, bbox@xmax, bbox@ymin),
                         start = yr,
                         end = yr,
                         frequency = "daily",
                         param = c("tmin","tmax"),
                         path = weather_path,
                         silent = FALSE)
    # this works with windows, IDK why
    download_daymet_tiles(tiles = 11207,
                          start = 2012,
                          end = 2012,
                          param = "tmin",
                          path = "daymet")
  }
}



# Derived parameters -----
REGION <- assign_extent(region_param = region_param)
params <- species_params(species, biotype, nsim, model_CDL)

# If GDD not pre-calculated, load list of PRISM files
# Add option for GDD pre-calculated (if all stages have same traits)
# send file names to for loop
pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
tminfiles <- list.files(path = weather_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = TRUE)
tminfiles <- ExtractBestPRISM(tminfiles, yr, leap = "keep")[start_doy:end_doy]

pattern = paste("(PRISM_tmax_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
tmaxfiles <- list.files(path = weather_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = TRUE)
tmaxfiles <- ExtractBestPRISM(tmaxfiles, yr, leap = "keep")[start_doy:end_doy]

# make template for all subsequent rasters to ensure same extent
# TODO: test that tmin and tmax have same result

template <- crop(raster(tminfiles[1]), REGION)
# template <- crop(aggregate(raster(files[1]), fact = 2), REGION)
template[!is.na(template)] <- 0

# splits template into list of 4 smaller map rasters to run in parallel
# only doing this for CONUS, because benefit lost for smaller map sizes
# TODO: alternative for speed would be increasing cell size by aggregation
if (region_param == "CONUS"){
  SplitMap <- SplitRas(template, ppside = 2, TRUE, FALSE)
}else{
  SplitMap <- list(template)
}

# create directory for model output
# named by date/time of model start
newname <- paste("run", gsub("-", "", 
                             gsub(" ", "_", 
                                  gsub(":", "", Sys.time(), fixed = TRUE),
                                  fixed = TRUE),
                             fixed = TRUE), sep = "")
dir.create(newname)

if ( runparallel == 1){
  # run simulations in parallel
  # make this option in parameter file?
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
}


# 3. Run model -----

# if errors with some sims/maps, rerun foreach loop below again for them
# use folders left in raster tmp directory to know which didn't work

# system.time({
outfiles <- foreach(sim = 1:nsim,
                    # outfiles <- foreach(sim = 5, # if some runs don't work, rerun individually
                    .packages= "raster",
                    .inorder = FALSE) %:% 
  foreach(maps = 1:length(SplitMap),
          # foreach(maps = 3, # if some runs don't work, rerun individually
          .packages = "raster",
          .inorder = FALSE) %dopar%{
            # .inorder = FALSE) %do%{
            
            # profvis::profvis({
            # creates unique filepath for directory to store results
            tmppath <- paste0(myrastertmp, "run", maps, sim)
            dir.create(path = tmppath, showWarnings = FALSE)
            # sets temp directory for raster files
            rasterOptions(tmpdir=file.path(tmppath)) 
            
            # Varying traits from parameter file
            stgorder   <- params$stgorder
            relpopsize <- params$relpopsize[sim]
            stage_dd   <- params$stage_dd[sim, ]
            stage_ldt  <- params$stage_ldt
            stage_udt  <- params$stage_udt
            photo_sens <- params$photo_sens
            CDL        <- params$CDL[1]
            cdl_b0     <- params$CDL[2]
            cdl_b1     <- params$CDL[3]
            
            # Initiate rasters that all match template so they have same extent
            template <- SplitMap[[maps]]
            # Track degree-day accumulation per cell per day
            dd_accum <- template
            # Overwintering as stage 1, all cells start here
            # Track lifestage for each cell per day
            lifestage <- template + 1 
            # Track voltinism and diapause per cell per day, starting at 0
            numgen <- template
            diap <- template
            
            # Loop over days ----
            for (index in seq_along(tminfiles)) {
              
              # get daily temperature, crop to REGION
              tmin <- crop(raster(tminfiles[index]), template)
              tmax <- crop(raster(tmaxfiles[index]), template)
              
              # photoperiod for this day across raster
              # TODO: could be done outside of loop for speed up
              if (model_CDL != 0){
                doy <- start_doy + index - 1
                photo <- RasterPhoto(template, doy, perc_twilight = 25)
              }

              # Assign lifestage raster -----
              # Each raster cell is assigned to a lifestage
              # These three rasters assign physiological parameters by cell
              ls_ldt <- setValues(lifestage, stage_ldt[getValues(lifestage)])
              ls_udt <- setValues(lifestage, stage_udt[getValues(lifestage)])
              ls_dd  <- setValues(lifestage, stage_dd[getValues(lifestage)])
              
              # Calculate stage-specific degree-days for each cell per day
              dd_tmp <- overlay(tmax, tmin, ls_ldt, ls_udt, fun = TriDD)
              rm(ls_ldt, ls_udt, tmax, tmin)
              # gc()
              
              # TODO: wrap inside function for memory efficiency
              # Accumulate degree days ----
              dd_accum <- dd_accum + dd_tmp
              # Calculate lifestage progression: Is accumulation > lifestage requirement
              progress <- dd_accum >= ls_dd
              lifestage <- lifestage + progress
              # Reset the DDaccum cells to zero for cells that progressed to next lifestage
              dd_accum <- dd_accum - (progress * ls_dd)
              # If adult stage progressed, that cell has oviposition
              ovip <- lifestage == (which(stgorder == "A") + 1)
              numgen <- numgen + ovip
              # Reassign cells with oviposition to egg stage
              lifestage <- Cond(ovip == 1, which(stgorder == "E"), lifestage)
              rm(dd_tmp, ovip, progress, ls_dd)
              # gc()
              
              # Diapause ----
              # keep running lifecycle as if all directly developing
              # but add new raster to track percent of population actually
              # in diapause following photoperiod response
              if (model_CDL == 1){
                sens_mask <- Cond(lifestage %in% photo_sens, 1, 0)
                prop_diap <- photo < CDL
                tmpdiap <- Cond(sens_mask == 1, prop_diap, diap)
                diap <- Cond(prop_diap < diap, diap, tmpdiap)
              }
              
              if (model_CDL == 2){
                # TODO: wrap inside function for memory efficiency
                # add logistic regression variation in CDL response
                sens_mask <- Cond(lifestage %in% photo_sens, 1, 0)
                # prop_diap <- exp(cdl_b0 + cdl_b1 * photo) /
                #   (1 + exp(cdl_b0 + cdl_b1 * photo))
                prop_diap <- calc(x = photo, function(x){
                  exp(cdl_b0 + cdl_b1 * x) /
                    (1 + exp(cdl_b0 + cdl_b1 * x))
                })
                tmpdiap <- Cond(sens_mask == 1, prop_diap, diap)
                # need to account for prop_diap declining up until solstice
                # only let proportion diapausing increase over season
                diap <- Cond(prop_diap < diap, diap, tmpdiap)
              }
              
              # memory management
              rm(sens_mask, prop_diap, tmpdiap)
              # gc()
              
              # TODO:
              # Add section to combine simulations into weighted values before stacking
              # Get one raster per day per 3 outputs first, fewer files to deal with later

              # stack each days raster
              if (!exists("LS_stack")){
                LS_stack <- stack(lifestage)
              } else {
                LS_stack <- addLayer(LS_stack, lifestage)
              }
              
              if (model_CDL != 0){
                diap1 <- round(diap * 1000)
                if (!exists("diap_stack")){
                  diap_stack <- stack(diap1)
                } else {
                  diap_stack <- addLayer(diap_stack, diap1)
                }    
              }
              
              if (!exists("numgen_stack")){
                numgen_stack <- stack(numgen)
              } else {
                numgen_stack <- addLayer(numgen_stack, numgen)
              }
              
              rm(diap1)
            }  # daily loop
            
            # Write results ----
            # if parallel by map chunk
            # each band is a day, first index is map chunk index
            mapcode <- formatC(maps, width = 3, format = "d", flag = "0")
            
            writeRaster(numgen_stack,
                        filename = paste(newname, "/NumGen_", mapcode, "_sim", sim, sep = ""),
                        overwrite = TRUE, datatype = "INT1U")
            writeRaster(LS_stack, filename = paste(newname, "/LS_", mapcode, "_sim", sim, sep = ""),
                        overwrite = TRUE, datatype = "INT1U")
            if (model_CDL != 0){
              writeRaster(diap_stack, filename = paste(newname, "/diap_", mapcode, "_sim", sim, sep = ""),
                          overwrite = TRUE, datatype = "INT2U")
            }
            # memory management
            rm(numgen_stack, LS_stack, diap_stack)
            # delete entire raster temp directory without affecting other running processes
            unlink(tmppath, recursive = TRUE)
            gc()
            
            # }) # profvis
          } # close foreach loop
# }) #system.time

if(exists("cl")){
  stopCluster(cl)
}


# TODO: rename results directory by species/yr/biotype/region

# 4. Process results ----- 
# Simulations combined to approximate continuous distributions
# Works slowly now (30 minutes), will speed up later
# Split combined lifestage raster into one/lifestage
# Mosaic CONUS if needed
# Remove unused files

# library(abind)

# use template to make sure NA values in right places
# will need to adjust results to fit datatype with integers
template <- crop(raster(tminfiles[1]), REGION)
template[!is.na(template)] <- 0
dataType(template) <- "INT2U"

# Input directories with results
newdirs <- c("DCA_2017_LL", "DCA_2017_TM")
for (newname in newdirs){
  # Weighted results by substage sizes
  # Not needed for older model with only one parameter per stage
  returnwd <- getwd()
  setwd(newname)
  
  f <-list.files()
  rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]
  rasfiles <- rasfiles[grep(pattern = "sim", x = rasfiles, fixed = TRUE)]
  
  # for each sim with unique information to save
  ls <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,1])
  maps <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,2])
  sims <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,3])
  sims <- gsub(pattern = ".grd", replacement = "", x = sims)
  nsim <- length(sims)
  # parallel backend for foreach loop
  if ( runparallel == 1){
    # run simulations in parallel
    # make this option in parameter file?
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
  }
  
  # run inside foreach now
  # this loop takes a lot of time
  foreach(i = ls,
          .packages= "raster") %dopar%{
            # foreach(m = maps,
            #         .packages = "raster") %dopar%{
            m <- maps[1]
            
            tmppath <- paste0("~/REPO/photovoltinism/rastertmp/", i, "run", m)
            dir.create(path = tmppath, showWarnings = FALSE)
            #sets temp directory
            rasterOptions(tmpdir=file.path(tmppath)) 
            fs <- sort(rasfiles[grep(pattern = i, x = rasfiles, fixed = TRUE)])
            fs <- sort(fs[grep(pattern = m, x = fs, fixed = TRUE)])
            tmp <- brick(fs[1])[[1]]
            tmp <- tmp + template # NA values added
            tmp[!is.na(tmp)] <- 0
            
            if (i == "diap"){
              ll <- replicate(nlayers(brick(fs[1])), tmp)
              blank <- brick(ll)
              for (j in 1:length(fs)){
                # Diapause, total / 1000
                ras_weighted <- round(brick(fs[j]) * params$relpopsize[j])
                blank <- overlay(blank, ras_weighted, fun=function(x,y) x + y)
              }
              outras <- writeRaster(blank, filename = paste(i, m, "weighted", sep = "_"),
                                    overwrite = TRUE, datatype = "INT2U")
            }
            
            # Lifestage
            if (i == "LS"){
              for (stg in seq_along(params$stgorder)){
                ll <- replicate(nlayers(brick(fs[1])), tmp)
                blank <- brick(ll)
                for (j in 1:length(fs)){
                  ras_weighted <- round(1000 * (brick(fs[j]) == stg) * params$relpopsize[j])
                  blank <- overlay(blank, ras_weighted, fun=function(x,y) x + y)
                }
                outras <- writeRaster(blank, filename = paste(params$stgorder[stg], m, "weighted", sep = "_"),
                                      overwrite = TRUE, datatype = "INT2U")
              }
            }
            
            
            if (i == "NumGen"){
              maxvolt <- max(unlist(lapply(fs,  function(x){
                max(getValues(brick(x)), na.rm = TRUE)
              })))
              for (stg in seq_len(maxvolt)){
                ll <- replicate(nlayers(brick(fs[1])), tmp)
                blank <- brick(ll)
                for (j in 1:length(fs)){
                  ras_weighted <- round(1000 * (brick(fs[j]) == stg) * params$relpopsize[j])
                  blank <- overlay(blank, ras_weighted, fun=function(x,y) x + y)
                }
                voltname <- paste0("NumGen", stg)
                outras <- writeRaster(blank, filename = paste(voltname, m, "weighted", sep = "_"),
                                      overwrite = TRUE, datatype = "INT2U")
              }
            }
            
            
            # might be faster, test later
            # test <- BigOverlayWeightedMean(fs, params$relpopsize)
            #   test <- test + template
            #   BigOverlayWeightedMean <- function(fs, simsize) {
            #     require(abind)
            #     require(raster)
            #     L <- as.list(fs)
            #     B <- lapply(seq_along(L), 
            #                  FUN = function(x, ras, pop){
            #                   brick(ras[x]) * pop[x]}, ras = L, pop = simsize)
            #     A <- do.call(abind, c(lapply(B, as.matrix), along=3))
            #     A <- aperm(A, c(3, 1, 2))
            #     z <- colSums(A)
            #     # z <- apply(A, c(1:2), FUN = weighted.mean, w = simsize)
            #     b <- setValues(brick(fs[1]), z)
            #   }
            
            removeTmpFiles(h = 0)
            unlink(tmppath, recursive = TRUE)
            
          } # end foreach
  
  setwd(returnwd)
} # end newdirs

stopCluster(cl) #WINDOWS
# setwd(returnwd)

# # quick plot of a few days from weighted raster to check
# test <- brick(paste(newname, "NumGen5_001_weighted.grd", sep = "/"))
# test <- brick(paste(newname, "L_001_weighted.grd", sep = "/"))
# test <- brick(paste(newname, "diap_001_weighted.grd", sep = "/"))
# plot(test[[c(150, 200, 250, 300)]])

# # mosaic maps together if CONUS
# # maybe hold off on this since its files so large
# # # LINUX
# # ncores <- length(ls)
# # registerDoMC(cores = ncores)
# # WINDOWS
# ncores <- 4
# cl <- makeCluster(ncores)
# registerDoSNOW(cl)
# f <-list.files()
# rasfiles <- f[grep(pattern = "weighted", x = f, fixed = TRUE)]
# lstage <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,1])
# 
# system.time({
#   foreach(i = lstage,
#           .packages = "raster") %dopar%{
#             fs <- rasfiles[grep(pattern = i, x = rasfiles, fixed = TRUE)]
#             
#             bricklist <- list()
#             for (m in 1:length(fs)){
#               bricklist[[m]] <- brick(fs[m])
#             }
#             
#             bricklist$filename <- paste(i, "all", sep = "_")
#             bricklist$overwrite <- TRUE
#             test <- do.call(merge, bricklist)
#           }
#   
#   cleanup <- list.files()
#   cleanup <- cleanup[-grep("all", x = cleanup)]
#   lapply(cleanup, FUN = file.remove) # CAREFUL HERE!
# })
# stopCluster(cl) #WINDOWS
# setwd(returnwd)


# 5. Clean up ------
# remove non-weighted results to save disk space
# make sure you have the correct directory!
cleanup <- list.files(newdirs[2], full.names = TRUE)
cleanup <- cleanup[-grep("weighted", x = cleanup)]
# lapply(cleanup, FUN = file.remove) # CAREFUL HERE!


