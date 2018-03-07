# Condensed DDRP model
# Universal code for 3 biocontrol species
# Any number of lifestages allowed

#####
# packages, options, functions loaded
#####
# TODO: find ways to take out dplyr, purrr, mixsmsn functions?
pkgs <- c("sp", "rgdal", "raster", "lubridate", "mixsmsn", "dplyr",
          "stringr", "purrr", "prism", "foreach", "doParallel")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them

# adjust raster options for your computer, likely this is fine
rasterOptions(overwrite = FALSE, 
              chunksize = 1e+07,
              maxmemory = 1e+08)

# load collection of functions for this model
source('CDL_funcs.R')
source('species_params.R')

##### 
#User input
#####
# directory with a year of PRISM tmax and tmin files
prism_path <- "prismDL/2017"
# prism_path <- "/data/PRISM/2014"
download_prism <- 0 # 1 if you need to download PRISM data first (20 minutes)
# directory to hold temporary raster results
myrastertmp <- "~/REPO/photovoltinism/rastertmp/"
# run simulations with parallel processing
runparallel <- 1 # 1 for yes, 0 for no

# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model scope
yr           <- 2017
start_doy    <- 200
end_doy      <- 250
region_param <- "WEST"
species      <- "DCA" # GCA/APHA/DCA
biotype      <- "Lovelock" # TODO: add options for each species

# introducing individual variation, tracked with simulation for each substage
# assign to 1 to match previous model versions
nsim <- 7 # number of substages/cohorts to approximate emergence distribution

# photoperiod decision inclusion
# 2 for logistic, 1 for single value CDL, 0 for none
model_CDL  <- 2 


# PRISM download if needed
# using ropensci 'prism' package to access webservice
# downloads entire CONUS, so files are large
if (download_prism == 1){
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



#####
# derived parameters
#####
REGION <- assign_extent(region_param = region_param)
params <- species_params(species, biotype, nsim, model_CDL)

# If GDD not pre-calculated, load list of PRISM files
# Add option for GDD pre-calculated (if all stages have same traits)
# send file names to for loop
pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
tminfiles <- list.files(path = prism_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = TRUE)
tminfiles <- ExtractBestPRISM(tminfiles, yr, leap = "keep")[start_doy:end_doy]

pattern = paste("(PRISM_tmax_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
tmaxfiles <- list.files(path = prism_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = TRUE)
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


# if errors with some sims/maps, rerun foreach loop below again for them
# use folders left in raster tmp directory to know which didn't work

#######
#Run model
#######

# system.time({
outfiles <- foreach(sim = 1:nsim,
                    # outfiles <- foreach(sim = 5, # if some runs don't work, rerun individually
                    .packages= "raster",
                    # .export = c("SplitMap", "tmaxfiles", "tminfiles",
                    #             "newname", "params"),
                    .inorder = FALSE) %:% 
  foreach(maps = 1:length(SplitMap),
          # foreach(maps = 3, # if some runs don't work, rerun individually
          .packages = "raster",
          # .export = c("SplitMap", "tmaxfiles", "tminfiles",
          #             "newname", "params"),
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
            
            # loop over days
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
              
              # Each raster cell is assigned to a lifestage
              # These three rasters assign physiological parameters by cell
              ls_ldt <- setValues(lifestage, stage_ldt[getValues(lifestage)])
              ls_udt <- setValues(lifestage, stage_udt[getValues(lifestage)])
              ls_dd  <- setValues(lifestage, stage_dd[getValues(lifestage)])
              
              # Calculate stage-specific degree-days for each cell per day
              dd_tmp <- overlay(tmax, tmin, ls_ldt, ls_udt, fun = TriDD)
              rm(ls_ldt, ls_udt, tmax, tmin)
              # gc()
              
              # Accumulate degree days
              # TODO: wrap inside function for memory efficiency
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
              
              # Diapause
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
                prop_diap <- calc(x = photo, fun = function(x){
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
              # But how to parallelize?
              
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
                                     overwrite = TRUE, datatype = "INT1U")
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
