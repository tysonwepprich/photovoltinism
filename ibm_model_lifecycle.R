# Speeding up raster-based model with arrays (cohort x pixel)
# Universal code for 3 biocontrol species
# Any number of lifestages allowed
# IBM version


# 1. Setup -------
# packages, options, functions loaded
# TODO: find ways to take out dplyr, purrr, mixsmsn functions?
pkgs <- c("sp", "rgdal", "rgeos", "raster", "lubridate", "mixsmsn", "dplyr",
          "stringr", "prism", "purrr", "foreach", "doParallel", "abind", "ncdf4")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them

# load collection of functions for this model
source('CDL_funcs.R')
source('species_params.R')

# 2. User input -----
# directory with daily tmax and tmin raster files
# weather_path <- "prism"
# weather_path <- "data/maca"
weather_path <- "prism" # PRISM data on grub server (needs to have stable files downloaded)
download_daily_weather <- 0 # 1 if you need to download PRISM/Daymet data first (20 minutes)
# weather_data_source <- "macav2" # could also have daymet or macav2
weather_data_source <- "prism"

# use parallel processing
# splits raster into chunks for running through simulation
runparallel <- 0 # 1 for yes, 0 for no
ncores <- 1 # choose number of cores for parallel\
time_interval <- 1 # number of days between results that are saved (7 = weekly), saves disk space and memory if not daily

# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model scope
yr           <- 2015
start_doy    <- 1
end_doy      <- 365
region_param <- "NORTHWEST" # TEST/WEST/EAST/CONUS/SOUTHWEST/NORTHWEST
species      <- "GCA" # GCA/APHA/DCA
biotype      <- "S" # TODO: add options for each species, N or S for APHA and GCA


# individual variation
ncohort <- 10 # number of cohorts to approximate emergence distribution

# photoperiod decision inclusion
# 2 for logistic, 1 for single value CDL, 0 for none
model_CDL  <- 0


# Weather download ----
# using ropensci 'prism' package to access webservice, or 'daymetr' package
# downloads entire CONUS, so files are large
if (download_daily_weather == 1){
  if (weather_data_source == "prism"){
    startdate <- as.Date(start_doy - 1, origin = paste0(yr, "-01-01"))
    enddate <- as.Date(end_doy - 1, origin = paste0(yr, "-01-01"))
    
    # need to set download directory
    yr_path <- paste(weather_path, yr, sep = "/")
    
    if(!exists(yr_path)){
      dir.create(yr_path)
    }
    
    # options(prism.path = paste("prismDL", yr, sep = "/"))
    options(prism.path = yr_path)
    
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
  tminfiles <- list.files(path = yr_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = TRUE)
  tminfiles <- ExtractBestPRISM(tminfiles, yr, leap = "keep")[start_doy:end_doy]
  
  pattern = paste("(PRISM_tmax_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
  tmaxfiles <- list.files(path = yr_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = TRUE)
  tmaxfiles <- ExtractBestPRISM(tmaxfiles, yr, leap = "keep")[start_doy:end_doy]
  
  geo_template <- crop(raster(tminfiles[1]), REGION)
  # template <- crop(aggregate(raster(files[1]), fact = 2), REGION)
  geo_template[!is.na(geo_template)] <- 0
}

if (weather_data_source == "macav2"){
  yr_path <- paste(weather_path, yr, sep = "/")
  
  tmin <- brick(paste(yr_path, list.files(yr_path)[grep(x = list.files(yr_path), pattern = "tasmin", fixed = TRUE)], sep = "/"), 
                varname = "air_temperature", lvar = 3, level = 4)
  tmax <- brick(paste(yr_path, list.files(yr_path)[grep(x = list.files(yr_path), pattern = "tasmax", fixed = TRUE)], sep = "/"),
                varname = "air_temperature", lvar = 3, level = 4)
  tminfiles <- shift(tmin[[which(lubridate::year(as.Date(getZ(tmin))) == yr)]], x = -360)
  tmaxfiles <- shift(tmax[[which(lubridate::year(as.Date(getZ(tmax))) == yr)]], x = -360)
  
  geo_template <- crop(tminfiles[[1]], REGION)
  # template <- crop(aggregate(raster(files[1]), fact = 2), REGION)
  geo_template[!is.na(geo_template)] <- 0
}


if ( runparallel == 1){
  # avoid memory issues on laptop
  if(ncores > (parallel::detectCores() / 2)){
    ncores <- parallel::detectCores() / 2
  }
  
  cl <- makePSOCKcluster(ncores)
  registerDoParallel(cl)
}
if ( runparallel == 0){
  ncores <- 1
}


chunk2 <- function(x,n){
  if(n > 1){
    split(x, cut(seq_along(x), n, labels = FALSE))
  }else{
    list(x)
  }
}

# testing one pixel at a time
ncores <- length(geo_template)
chunk_list <- chunk2(getValues(geo_template), ncores)
# choose one pixel to run
chunk <- which(round(getValues(init(geo_template, "x")), 1) == -119.5 & round(getValues(init(geo_template, "y")), 1) == 46.5)[1]

# function needed to combine array results with foreach
acomb4 <- function(...) abind::abind(..., along = 4)
acomb1 <- function(...) abind::abind(..., along = 1)


# 3. Run lifestage model -----

test <- system.time({
  
  # outlist <- foreach(chunk = 1:ncores,
  #                    .packages = "raster",
  #                    .inorder = TRUE,
  #                    .combine = 'acomb1',
  #                    .multicombine = TRUE) %dopar%{
  #                      
                       
                       template <- chunk_list[[chunk]]
                       
                       template_mat <- rep(template, times = ncohort)
                       
                       # latitudes for photoperiod model
                       lats <- chunk2(getValues(init(geo_template, 'y')), ncores)[[chunk]]
                       
                       
                       # create array to hold results
                       # output variables: lifestage, numgen, daily degree-days
                       ls_array <- array(data = NA, dim = c(
                         length(seq(start_doy, end_doy, by = time_interval)), 
                         length(params$stgorder) + 6)) 
                       
                       # Varying traits from parameter file
                       stgorder   <- params$stgorder
                       # relpopsize <- params$relpopsize
                       stage_dd   <- params$stage_dd
                       stage_ldt  <- params$stage_ldt
                       stage_udt  <- params$stage_udt
                       photo_sens <- params$photo_sens
                       CDL        <- params$CDL[1]
                       cdl_b0     <- params$CDL[2]
                       cdl_b1     <- params$CDL[3]
                       
                       # vector of which individuals to simulate from data.frame
                       ind_df <- data.frame(stage_dd)
                       names(ind_df) <- stgorder
                       
                       ls_complete <- data.frame(array(data = 0, dim = dim(stage_dd)))
                       names(ls_complete) <- paste(stgorder, "DOY", sep = "_")
                       
                       ind_df <- cbind(ind_df, ls_complete)
                       ind_df$dd_accum <- 0
                       ind_df$dd_total <- 0
                       ind_df$first_frost <- 0
                       ind_df$last_frost <- 0
                       ind_df$diapause <- -1
                       ind_df$numgen <- 0
                       ind_df$parent <- NA
                       ind_df$lifestage <- 1
                       ind_df$indiv <- 1:nrow(ind_df)
                       ind_df$active <- 1    # this will indicate rows to simulate each day
                       
                       # Loop over days ----
                       for (index in 1:length(start_doy:end_doy)) {
                         
                         ncohort <- sum(ind_df$active)
                         if (ncohort == 0) break
                         
                         # data.frames to track individuals
                         newdf <- ind_df[which(ind_df$active == 1), ]
                         olddf <- ind_df[-which(ind_df$active == 1), ]
                         
                         lifestage <- newdf$lifestage
                         dd_accum <- newdf$dd_accum
                         dd_total <- newdf$dd_total
                         diap <- newdf$diapause
                         stage_dd <- as.matrix(newdf[, stgorder])
                         
                         doy <- start_doy + index - 1
                         # # get daily temperature, crop to REGION
                         if (weather_data_source == "prism"){
                           tmin <- chunk2(getValues(crop(raster(tminfiles[index]), geo_template)), ncores)[[chunk]]
                           tmax <- chunk2(getValues(crop(raster(tmaxfiles[index]), geo_template)), ncores)[[chunk]]
                         }
                         if (weather_data_source == "macav2"){
                           tmin <- chunk2(getValues(crop(tminfiles[[index]], geo_template)), ncores)[[chunk]]
                           tmax <- chunk2(getValues(crop(tmaxfiles[[index]], geo_template)), ncores)[[chunk]]
                         }
                         # frost boundaries to growing season
                         if(doy <= 182){
                           newdf$last_frost <- Cond(tmin <= -2, doy, newdf$last_frost)
                         }
                         if(doy >= 183){
                           newdf$first_frost <- Cond(newdf$first_frost == 0, (tmin <= -2) * doy, newdf$first_frost)
                         }
                         
                         tmax <- rep(tmax, times = ncohort) 
                         tmin <- rep(tmin, times = ncohort)
                         
                         
                         # Assign lifestage raster -----
                         # Each raster cell is assigned to a lifestage
                         # These three rasters assign physiological parameters by cell
                         ls_ldt <- stage_ldt[as.vector(lifestage)] 
                         ls_udt <- stage_udt[as.vector(lifestage)] 
                         ls_dd  <- sapply(1:ncohort, FUN = function(x) stage_dd[x, ][lifestage[x]])
                         
                         # Calculate stage-specific degree-days for each cell per day
                         dd_tmp <- as.vector(TriDD(tmax, tmin, ls_ldt, ls_udt))
                         # Accumulate degree days ----
                         dd_accum <- dd_accum + dd_tmp # resets with lifestage change
                         newdf$dd_total <- dd_total + dd_tmp[1] # tracks annual accumulation
                         
                         # Calculate lifestage progression: Is accumulation > lifestage requirement
                         progress <- dd_accum >= ls_dd
                         # DOY of progress tracked
                         for (ls in 1:length(stgorder)){
                          newdf[, paste(stgorder[ls], "DOY", sep = "_")] <- Cond(progress == 1 & lifestage == ls, 
                                                                                         index, 
                                                                                         newdf[, paste(stgorder[ls], "DOY", sep = "_")])
                         }
                         
                         
                         
                         # 4. Diapause ----
                         # Take lifestage result from degree-day model and estimate diapause over the year
                         # Need to know photo-sensitive lifestage
                         
                         # Only lets individual make diapause decision once, even if in stage for longer time
                         
                         # photoperiod for this day across vector
                         if (model_CDL != 0){
                           photo <- photoperiod(lats, doy, p = 1.5)
                         }
                         
                         if (model_CDL == 1){
                           sens_mask <- Cond(lifestage %in% photo_sens & diap == -1, 1, 0)
                           prop_diap <- photo < CDL
                           diap <- Cond(sens_mask == 1, rbinom(ncohort, size = 1, prob = prop_diap), diap)
                         }
                         
                         if (model_CDL == 2){
                           sens_mask <- Cond(lifestage %in% photo_sens & diap == -1, 1, 0)
                           prop_diap <- exp(cdl_b0 + cdl_b1 * photo) /
                             (1 + exp(cdl_b0 + cdl_b1 * photo))
                           diap <- Cond(sens_mask == 1, rbinom(ncohort, size = 1, prob = prop_diap), diap)
                         }
                         newdf$diapause <- diap
                         
                         # track individual progress, diapause decisions, and oviposition/mortality
                         ovip <- Cond(progress == 1 & 
                                        lifestage %in% which(stgorder %in% c("A", "OA")) &
                                        diap < 1,
                                      newdf$indiv,
                                      0)
                         lifestage <- Cond(progress == 1 & ovip == 0, lifestage + progress, lifestage)
                         # Reassign cells that progressed past end of stgorder to first non-overwintering stage
                         lifestage <- Cond(lifestage == (length(stgorder) + 1), 2, lifestage)
                         newdf$lifestage <- lifestage
                         # Inactive if decided to diapause and reaches overwinter stage
                         newdf$active[which(newdf$lifestage == (length(stgorder)) & newdf$diapause == 1)] <- 0
                         
                         dd_accum <- dd_accum - (progress * ls_dd) # should this reset to zero? (assuming no accumulation if progress)
                         newdf$dd_accum <- dd_accum
                         
                         if(max(ovip) > 0){
                           newrows <- newdf[which(newdf$indiv %in% ovip), ]
                           if(nrow(newrows) > 0){
                             
                             # TODO:
                             # Draw degree-day stage requirements here for variation 
                             
                             newrows[, paste(stgorder, "DOY", sep = "_")] <- 0
                             newrows$diapause <- -1
                             newrows$numgen <- newrows$numgen + 1
                             newrows$parent <- ovip[which(ovip > 0)]
                             newrows$lifestage <- which(stgorder == "E")
                             newrows$dd_accum <- 0
                             newrows$indiv <- (max(ind_df$indiv) + 1):(max(ind_df$indiv) + nrow(newrows))
                             newrows$active <- 1
                           }
                           
                           newdf <- rbind(newdf, newrows)
                           newdf$active[which(newdf$indiv %in% ovip)] <- 0
                         }
                         
                         ind_df <- rbind(olddf, newdf)
                         
                         
                       #   # weights each cohort now with relpopsize and colSums
                       #   # time interval called week, but could be any number of days for saving results
                       #   if (index %in% seq(start_doy, end_doy, by = time_interval)){
                       #     week <- which(seq(start_doy, end_doy, by = time_interval) == index)
                       #     for (i in 1:length(stgorder)){
                       #       ls_array[, week, i] <- colSums(relpopsize * (lifestage == i), na.rm = FALSE)
                       #     }
                       #     ls_array[, week, length(stgorder) + 1] <-  colSums(relpopsize * numgen, na.rm = FALSE) # generations by oviposition
                       #     ls_array[, week, length(stgorder) + 2] <-  colSums(relpopsize * fullgen, na.rm = FALSE) # max generations potentially completed
                       #     ls_array[, week, length(stgorder) + 3] <-  colSums(relpopsize * diap, na.rm = FALSE)
                       #     # weather data, doesn't really need to be weekly, but easiest to include in same array
                       #     ls_array[, week, length(stgorder) + 4] <- dd_total
                       #     ls_array[, week, length(stgorder) + 5] <- last_frost
                       #     ls_array[, week, length(stgorder) + 6] <- first_frost
                       #   }
                         
                       }  # daily loop
                       
                       # return(ind_df)
                       
                     # } # close foreach loop
}) # system.time

# try to stop cluster, otherwise old rstudio sessions stay on server indefinitely
stopCluster(cl)

fullgen <- fullgen + (progress == 1 & lifestage == (length(stgorder) - 1)) # reaches OW stage



# check results by plotting values assigned to raster
# array dimensions: [pixel, week, results]
plot(setValues(geo_template, outlist[, 53, 5]))

# plot time series of lifestage at a single pixel for weighted cohorts
plot(outlist[50000, , 1])



