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
weather_path <- "/data/PRISM"
download_daily_weather <- 0 # 1 if you need to download PRISM/Daymet data first (20 minutes)
weather_data_source <- "prism"

# directory to hold temporary raster results
myrastertmp <- "~/REPO/photovoltinism/rastertmp/"
# run simulations with parallel processing
runparallel <- 1 # 1 for yes, 0 for no
ncores <- 5 # choose number of cores for parallel

# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model scope
yr           <- 2015
start_doy    <- 1
end_doy      <- 365
region_param <- "CONUS" # TEST/WEST/EAST/CONUS/SOUTHWEST/NORTHWEST
species      <- "APHA" # GCA/APHA/DCA
biotype      <- "S" # TODO: add options for each species, N or S for APHA and GCA


# introducing individual variation, tracked with simulation for each substage
# assign to 1 to match previous model versions
ncohort <- 7 # number of substages/cohorts to approximate emergence distribution

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
  yr_path <- paste(weather_path, yr, sep = "/")
  # If GDD not pre-calculated, load list of PRISM files
  # Add option for GDD pre-calculated (if all stages have same traits)
  # send file names to for loop
  pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
  tminfiles <- list.files(path = yr_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = FALSE)
  tminfiles <- ExtractBestPRISM(tminfiles, yr, leap = "keep")[start_doy:end_doy]
  
  pattern = paste("(PRISM_tmax_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
  tmaxfiles <- list.files(path = yr_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = FALSE)
  tmaxfiles <- ExtractBestPRISM(tmaxfiles, yr, leap = "keep")[start_doy:end_doy]
}

if(weather_data_source == "prism"){
  geo_template <- crop(raster(tminfiles[1]), REGION)
  # template <- crop(aggregate(raster(files[1]), fact = 2), REGION)
  geo_template[!is.na(geo_template)] <- 0
}


if ( runparallel == 1){
  # run simulations in parallel
  
  # avoid memory issues on laptop
  if(ncores > (parallel::detectCores() / 2)){
    ncores <- parallel::detectCores() / 2
  }
  
  cl <- makePSOCKcluster(ncores)
  registerDoParallel(cl)
}

chunk2 <- function(x,n){
  if(n > 1){
    split(x, cut(seq_along(x), n, labels = FALSE))
  }else{
    list(x)
  }
}
chunk_list <- chunk2(getValues(geo_template), ncores)


# function needed to combine array results with foreach
acomb4 <- function(...) abind::abind(..., along = 4)
acomb1 <- function(...) abind::abind(..., along = 1)


# 3. Run lifestage model -----

test <- system.time({
  
  outlist <- foreach(chunk = 1:ncores,
                     .packages = "raster",
                     .inorder = TRUE,
                     .combine = 'acomb1',
                     .multicombine = TRUE) %dopar%{
                       
                       
                       template <- chunk_list[[chunk]]
                       
                       template_mat <- matrix(template, 
                                              nrow = ncohort,
                                              ncol = length(template), 
                                              byrow=TRUE)
                       
                       
                       
                       # create array to hold results
                       # only track weekly results to save disk space/memory
                       # output variables: lifestage, numgen, daily degree-days
                       ls_array <- array(data = NA, dim = c(length(template), 
                                                            length(seq(start_doy, end_doy, by = 7)), 
                                                            length(params$stgorder) + 3)) 
                       
                       # Varying traits from parameter file
                       stgorder   <- params$stgorder
                       relpopsize <- params$relpopsize
                       stage_dd   <- params$stage_dd
                       stage_ldt  <- params$stage_ldt
                       stage_udt  <- params$stage_udt
                       photo_sens <- params$photo_sens
                       CDL        <- params$CDL[1]
                       cdl_b0     <- params$CDL[2]
                       cdl_b1     <- params$CDL[3]
                       
                       # Track degree-day accumulation per cell per day
                       dd_accum <- template_mat
                       # dd_total <- template_mat
                       # Overwintering as stage 1, all cells start here
                       # Track lifestage for each cell per day
                       lifestage <- template_mat + 1
                       # Track voltinism and diapause per cell per day, starting at 0
                       numgen <- template_mat
                       fullgen <- template_mat
                       diap <- template_mat
                       
                       # Loop over days ----
                       for (index in 1:length(start_doy:end_doy)) {
                         doy <- start_doy + index - 1
                         # # get daily temperature, crop to REGION
                         tmin <- chunk2(getValues(crop(raster(tminfiles[index]), geo_template)), ncores)[[chunk]]
                         tmax <- chunk2(getValues(crop(raster(tmaxfiles[index]), geo_template)), ncores)[[chunk]]
                         lats <- chunk2(getValues(init(geo_template, 'y')), ncores)[[chunk]]
                         
                         tmax <- matrix(tmax, 
                                        nrow = ncohort,
                                        ncol = length(template), 
                                        byrow=TRUE)
                         tmin <- matrix(tmin, 
                                        nrow = ncohort,
                                        ncol = length(template), 
                                        byrow=TRUE)
                         
                         # Assign lifestage raster -----
                         # Each raster cell is assigned to a lifestage
                         # These three rasters assign physiological parameters by cell
                         ls_ldt <- matrix(stage_ldt[as.vector(lifestage)], 
                                          nrow = ncohort,
                                          ncol = length(template), 
                                          byrow=FALSE)
                         ls_udt <- matrix(stage_udt[as.vector(lifestage)], 
                                          nrow = ncohort,
                                          ncol = length(template), 
                                          byrow=FALSE)
                         ls_dd  <- matrix(stage_dd[as.vector(lifestage)],
                                          nrow = ncohort,
                                          ncol = length(template), 
                                          byrow=FALSE)
                         
                         # Calculate stage-specific degree-days for each cell per day
                         dd_tmp <- TriDD(tmax, tmin, ls_ldt, ls_udt)
                         # TODO: wrap inside function for memory efficiency
                         # Accumulate degree days ----
                         dd_accum <- dd_accum + dd_tmp # resets with lifestage change
                         # dd_total <- dd_total + dd_tmp # tracks annual accumulation
                         
                         # Calculate lifestage progression: Is accumulation > lifestage requirement
                         progress <- dd_accum >= ls_dd
                         # If reproductive adult stage progressed, that cell has oviposition and generation count increases
                         numgen <- numgen + (progress == 1 & lifestage == which(stgorder == "A")) 
                         fullgen <- fullgen + (progress == 1 & lifestage == (length(stgorder) - 1)) # reaches OW stage
                         # Reset the DDaccum cells to zero for cells that progressed to next lifestage
                         lifestage <- lifestage + progress
                         dd_accum <- dd_accum - (progress * ls_dd)
                         
                         # Reassign cells that progressed past end of stgorder to first non-overwintering stage
                         lifestage <- Cond(lifestage == (length(stgorder) + 1), 2, lifestage)
                         
                         # 4. Diapause ----
                         # Take lifestage result from degree-day model and estimate diapause over the year
                         # Need to know photo-sensitive lifestage
                         
                         # photoperiod for this day across vector
                         if (model_CDL != 0){
                           photo <- photoperiod(lats, doy, p = 1.5)
                         }
                         
                         if (model_CDL == 1){
                           sens_mask <- matrix(Cond(as.vector(lifestage) %in% photo_sens, 1, 0),
                                               nrow = ncohort,
                                               ncol = length(template), 
                                               byrow = FALSE)
                           prop_diap <- photo < CDL
                           tmpdiap <- Cond(sens_mask == 1, prop_diap, diap)
                           diap <- Cond(prop_diap < diap, diap, tmpdiap)
                         }
                         
                         if (model_CDL == 2){
                           # TODO: wrap inside function for memory efficiency
                           # add logistic regression variation in CDL response
                           sens_mask <- matrix(Cond(as.vector(lifestage) %in% photo_sens, 1, 0),
                                               nrow = ncohort,
                                               ncol = length(template), 
                                               byrow= FALSE)
                           prop_diap <- exp(cdl_b0 + cdl_b1 * photo) /
                             (1 + exp(cdl_b0 + cdl_b1 * photo))
                           tmpdiap <- Cond(sens_mask == 1, prop_diap, diap)
                           # need to account for prop_diap declining up until solstice
                           # only let proportion diapausing increase over season
                           diap <- Cond(prop_diap < diap, diap, tmpdiap)
                         }
                         
                         
                         # weights each cohort now with relpopsize and colSums
                         
                         if (index %in% seq(start_doy, end_doy, by = 7)){
                           week <- which(seq(start_doy, end_doy, by = 7) == index)
                           for (i in 1:length(stgorder)){
                             ls_array[, week, i] <- colSums(relpopsize * (lifestage == i), na.rm = FALSE)
                           }
                           ls_array[, week, length(stgorder) + 1] <-  colSums(relpopsize * numgen, na.rm = FALSE) # generations by oviposition
                           ls_array[, week, length(stgorder) + 2] <-  colSums(relpopsize * fullgen, na.rm = FALSE) # max generations potentially completed
                           ls_array[, week, length(stgorder) + 3] <-  colSums(relpopsize * diap, na.rm = FALSE)
                         }
                       }  # daily loop 
                       
                       
                       return(ls_array)
                       
                     } # close foreach loop
}) # system.time

# try to stop cluster, otherwise old rstudio sessions stay on server indefinitely
stopCluster(cl)


# check results by plotting values assigned to raster
# array dimensions: [pixel, day, results]
plot(setValues(geo_template, outlist[, 53, 3]))

gdd <- setValues(template, outlist[, 365, 4, 1])
outlist <- outlist[,,-4,]
# plot time series of lifestage at a single pixel for cohort 1
plot(outlist[50000, , 1, 1])

# weight the cohorts for each lifestage
cl <- makePSOCKcluster(15)
registerDoParallel(cl)
test <- system.time({
  eggs <- Cond(outlist[, , 1 , ] == 2, 1, 0)
  
  eggs <- apply(eggs, MARGIN = c(1, 2), FUN = function(x) x + params$relpopsize)
  
  # w_array <- parApply(cl = cl, X = eggs, MARGIN = c(1:2), FUN = stats::weighted.mean, w = params$relpopsize, na.rm = TRUE)
  w_array <- parApply(cl = cl, X = eggs, MARGIN = 2, FUN = rowSums, na.rm = TRUE)
  stopCluster(cl)
})


# TODO: 
# -arrays that write to disk so R doesn't run out of memory (ff/bigmemory pkg)
# -weighted averages by dimension in the array for cohorts
# -not saving every single day in the results to save disk space










res_array[, index, 4] <- dd_tmp 

