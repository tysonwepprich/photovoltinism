# Speeding up raster-based model
# Universal code for 3 biocontrol species
# Any number of lifestages allowed

# Changes:
# Arrays instead of rasters
# Will try two-stage model: 1. Degree-days and lifestage 2. Photoperiod/diapause 
# Hope to make it compatible with individual-based models, simulate many individuals in place of cohorts

# 1. Setup -------
# packages, options, functions loaded
# TODO: find ways to take out dplyr, purrr, mixsmsn functions?
pkgs <- c("sp", "rgdal", "rgeos", "raster", "lubridate", "mixsmsn", "dplyr",
          "stringr", "prism", "purrr", "foreach", "doParallel", "abind")
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
weather_path <- "prism"
download_daily_weather <- 0 # 1 if you need to download PRISM/Daymet data first (20 minutes)
weather_data_source <- "prism"

# directory to hold temporary raster results
myrastertmp <- "~/REPO/photovoltinism/rastertmp/"
# run simulations with parallel processing
runparallel <- 1 # 1 for yes, 0 for no
ncores <- 15 # choose number of cores for parallel

# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model scope
yr           <- 2017
start_doy    <- 1
end_doy      <- 365
region_param <- "NORTHWEST" # TEST/WEST/EAST/CONUS/SOUTHWEST/NORTHWEST
species      <- "APHA" # GCA/APHA/DCA
biotype      <- "S" # TODO: add options for each species, N or S for APHA and GCA


# introducing individual variation, tracked with simulation for each substage
# assign to 1 to match previous model versions
ncohort <- 3 # number of substages/cohorts to approximate emergence distribution

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
    options(prism.path = weather_path)
    
    get_prism_dailys("tmin", minDate = startdate, 
                     maxDate = enddate, keepZip = FALSE)
    get_prism_dailys("tmax", minDate = startdate, 
                     maxDate = enddate, keepZip = FALSE)
  }
}



# Derived parameters -----
REGION <- assign_extent(region_param = region_param)
params <- species_params(species, biotype, ncohort, model_CDL)

if (weather_data_source == "prism"){
  
  # If GDD not pre-calculated, load list of PRISM files
  # Add option for GDD pre-calculated (if all stages have same traits)
  # send file names to for loop
  pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
  tminfiles <- list.files(path = weather_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = TRUE)
  tminfiles <- ExtractBestPRISM(tminfiles, yr, leap = "keep")[start_doy:end_doy]
  
  pattern = paste("(PRISM_tmax_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
  tmaxfiles <- list.files(path = weather_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = TRUE)
  tmaxfiles <- ExtractBestPRISM(tmaxfiles, yr, leap = "keep")[start_doy:end_doy]
}

if(weather_data_source == "prism"){
  template <- crop(raster(tminfiles[1]), REGION)
  # template <- crop(aggregate(raster(files[1]), fact = 2), REGION)
  template[!is.na(template)] <- 0
}


if ( runparallel == 1){
  # run simulations in parallel
  
  # avoid memory issues on laptop
  if(ncores > (parallel::detectCores() - 1)){
    ncores <- parallel::detectCores() - 1
  }
  
  # cores will be allocated by cohort and chunk of raster pixels
  nchunk <- floor(ncores / ncohort)
  
  cl <- makePSOCKcluster(ncores)
  registerDoParallel(cl)
}

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

chunk_list <- chunk2(getValues(template), nchunk)

# function needed to combine array results with foreach
acomb4 <- function(...) abind::abind(..., along = 4)
acomb1 <- function(...) abind::abind(..., along = 1)


# 3. Run lifestage model -----

test <- system.time({
  
  outlist <-  foreach(cohort = 1:ncohort,
                      .packages = "raster",
                      .inorder = TRUE,
                      .combine = 'acomb4',
                      .multicombine = TRUE) %:% 
    foreach(chunk = 1:nchunk,
            .packages = "raster",
            .inorder = TRUE,
            .combine = 'acomb1',
            .multicombine = TRUE) %dopar%{
              
              chunk_template <- chunk_list[[chunk]]
              
              # create array to hold results
              # output variables: lifestage, numgen, daily degree-days
              ls_array <- array(data = NA, dim = c(length(chunk_template),
                                                   length(start_doy:end_doy), 3))
              
              # Varying traits from parameter file
              stgorder   <- params$stgorder
              relpopsize <- params$relpopsize[cohort]
              stage_dd   <- params$stage_dd[cohort, ]
              stage_ldt  <- params$stage_ldt
              stage_udt  <- params$stage_udt
              photo_sens <- params$photo_sens
              CDL        <- params$CDL[1]
              cdl_b0     <- params$CDL[2]
              cdl_b1     <- params$CDL[3]
              
              # Initiate vectors (not rasters) that all match template
              # If you want to see them as rasters: plot(setValues(template, .vector.))
              template_vec <- chunk_template
              # Track degree-day accumulation per cell per day
              dd_accum <- template_vec
              # Overwintering as stage 1, all cells start here
              # Track lifestage for each cell per day
              lifestage <- template_vec + 1
              # Track voltinism and diapause per cell per day, starting at 0
              numgen <- template_vec
              diap <- template_vec
              
              # Loop over days ----
              for (index in start_doy:end_doy) {
                
                # # get daily temperature, crop to REGION
                tmin <- chunk2(getValues(crop(raster(tminfiles[index]), template)), nchunk)[[chunk]]
                tmax <- chunk2(getValues(crop(raster(tmaxfiles[index]), template)), nchunk)[[chunk]]
                
                # Assign lifestage raster -----
                # Each raster cell is assigned to a lifestage
                # These three rasters assign physiological parameters by cell
                ls_ldt <- stage_ldt[lifestage]
                ls_udt <- stage_udt[lifestage]
                ls_dd  <- stage_dd[lifestage]
                
                # Calculate stage-specific degree-days for each cell per day
                dd_tmp <- TriDD(tmax, tmin, ls_ldt, ls_udt)
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
                
                
                # fill in array for each day's results
                ls_array[, index, 1] <- lifestage 
                ls_array[, index, 2] <- numgen
                ls_array[, index, 3] <- dd_tmp
                
                
              }  # daily loop
              return(ls_array)
              
            } # close foreach loop
}) # system.time
stopCluster(cl)

# check results by plotting values assigned to raster
# array dimensions: [pixel, day, results, cohort]
# results are 1. lifestage 2. Number of Generation 3. Daily degree-days
plot(setValues(template, outlist[, 360, 1, 1]))

# plot time series of lifestage at a single pixel for cohort 1
plot(outlist[50000, , 1, 1])

# weight the cohorts for each lifestage
cl <- makePSOCKcluster(15)
registerDoParallel(cl)
test <- system.time({
  eggs <- Cond(outlist[, , 1 , ] == 2, 1, 0)
  w_array <- parApply(cl = cl, X = eggs, MARGIN = c(1:2), FUN = stats::weighted.mean, w = params$relpopsize, na.rm = TRUE)
  stopCluster(cl)
})


# TODO: 
# -arrays that write to disk so R doesn't run out of memory (ff/bigmemory pkg)
# -weighted averages by dimension in the array for cohorts
# -not saving every single day in the results to save disk space





# 4. Diapause ----
# Take lifestage result from degree-day model and estimate diapause over the year
# Need to know photo-sensitive lifestage


# photoperiod for this day across raster
# this still needs raster to extract coordinates for latitude
if (model_CDL != 0){
  doy <- start_doy + index - 1
  photo <- getValues(RasterPhoto(template, doy, perc_twilight = 25))
}

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
  prop_diap <- exp(cdl_b0 + cdl_b1 * photo) /
    (1 + exp(cdl_b0 + cdl_b1 * photo))
  tmpdiap <- Cond(sens_mask == 1, prop_diap, diap)
  # need to account for prop_diap declining up until solstice
  # only let proportion diapausing increase over season
  diap <- Cond(prop_diap < diap, diap, tmpdiap)
}

res_array[, index, 3] <- diap 
res_array[, index, 4] <- dd_tmp 

