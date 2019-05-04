# Speeding up raster-based model with arrays (cohort x pixel)
# Universal code for 3 biocontrol species
# Any number of lifestages allowed
# Multiple years run at once, 5 yr climate forecasts for example

# TODO:
# 1. Download stable PRISM files to GRUB, many months/years not updated
# 2. MACAV2 option for forecasts
# 3. Redo mapping functions to use new results format
# 4. Add CDL values to foreach loop with goal to find optimal photoperiod response for each pixel

# Run Time
# 40 minutes for 5 years of Northwest region with 7 cohorts, weekly output save 2.3 GB


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
weather_path <- "data/maca"
# weather_path <- "/data/PRISM" # PRISM data on grub server (needs to have stable files downloaded)
weather_data_source <- "macav2" # could also have daymet or macav2
# weather_data_source <- "prism"

if(weather_data_source == "macav2"){
  # in Kelvin, adjust thresholds/frost or add 273.15 to giant rasters?
  weather_files <- list.files(path = weather_path, full.names = TRUE, recursive = TRUE)
}

# use parallel processing
# splits raster into chunks for running through simulation
runparallel <- 1 # 1 for yes, 0 for no
ncores <- 10 # choose number of cores for parallel

# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model scope
yr           <- 2016 # if multi-year, will be replaced
start_doy    <- 1
end_doy      <- 365
region_param <- "NW_SMALL" # TEST/WEST/EAST/CONUS/SOUTHWEST/NORTHWEST
# if CONUS, wouldn't need to crop PRISM or MACAv2
species      <- "APHA" # GCA/APHA/DCA
biotype      <- "S" # TODO: add options for each species, N or S for APHA and GCA


# introducing cohort variation
# assign to 1 to match previous model versions
ncohort <- 7 # number of cohorts to approximate emergence distribution

# photoperiod decision inclusion
# 2 for logistic, 1 for single value CDL, 0 for none
model_CDL  <- 2 

# Derived parameters -----
REGION <- assign_extent(region_param = region_param)
params <- species_params(species, biotype, ncohort, model_CDL)

# PRISM data is in one file per day, one directory per year
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

# MACAv2 files have 5 years in one netcdf file
if (weather_data_source == "macav2"){
  yr_path <- paste(weather_path, yr, sep = "/")
  tminbrick <- brick(paste(yr_path, list.files(yr_path)[grep(x = list.files(yr_path), pattern = "tasmin", fixed = TRUE)], sep = "/"), 
                varname = "air_temperature", lvar = 3, level = 4)
  tmaxbrick <- brick(paste(yr_path, list.files(yr_path)[grep(x = list.files(yr_path), pattern = "tasmax", fixed = TRUE)], sep = "/"),
                varname = "air_temperature", lvar = 3, level = 4)
  
  macafile <- basename(tmaxbrick@file@name)
  yrs <- gregexpr(pattern = "[0-9]{4}", text = macafile)
  yr1 <- as.numeric(substr(macafile, start = yrs[[1]][1], stop = yrs[[1]][1] + 3))
  yr2 <- as.numeric(substr(macafile, start = yrs[[1]][2], stop = yrs[[1]][2] + 3))
  yr <- yr1:yr2
  
  tminbrick <- shift(tminbrick, x = -360)
  tmaxbrick <- shift(tmaxbrick, x = -360)
  
  geo_template <- crop(tminbrick[[1]], REGION)
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
chunk_list <- chunk2(getValues(geo_template), ncores)


# function needed to combine array results with foreach
acomb4 <- function(...) abind::abind(..., along = 4)
acomb1 <- function(...) abind::abind(..., along = 1)


# 3. Run lifestage model -----

test <- system.time({
  
  # parallel by yr if MACA 5-year dataset
  outlist <- foreach(y = 1:length(yr),
                     .packages = c("raster", "lubridate"), 
                     .inorder = TRUE,
                     .combine = 'acomb4',
                     .multicombine = TRUE) %:%
    foreach(chunk = 1:ncores,
                     .packages = c("raster", "lubridate"),
                     .inorder = TRUE,
                     .combine = 'acomb1',
                     .multicombine = TRUE) %dopar%{
                       
                       modyear <- yr[y]
                       
                       if (weather_data_source == "macav2"){
                         tminfiles <- tminbrick[[which(lubridate::year(as.Date(getZ(tminbrick))) == modyear)]]
                         tmaxfiles <- tmaxbrick[[which(lubridate::year(as.Date(getZ(tmaxbrick))) == modyear)]]
                       }

                       template <- chunk_list[[chunk]]
                       
                       template_mat <- matrix(template, 
                                              nrow = ncohort,
                                              ncol = length(template), 
                                              byrow=TRUE)
                       
                       # latitudes for photoperiod model
                       lats <- chunk2(getValues(init(geo_template, 'y')), ncores)[[chunk]]
                       
                       
                       # create array to hold results
                       # only track weekly results to save disk space/memory
                       # output variables: lifestage, numgen, daily degree-days
                       ls_array <- array(data = NA, dim = c(length(template), 
                                                            length(seq(start_doy, end_doy, by = 7)), 
                                                            length(params$stgorder) + 6)) 
                       
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
                       # Weather variables to track (no cohort dimension needed)
                       dd_total <- template
                       first_frost <- template
                       last_frost <- template
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
                         if (weather_data_source == "prism"){
                           tmin <- chunk2(getValues(crop(raster(tminfiles[index]), geo_template)), ncores)[[chunk]]
                           tmax <- chunk2(getValues(crop(raster(tmaxfiles[index]), geo_template)), ncores)[[chunk]]
                         }
                         if (weather_data_source == "macav2"){
                           tmin <- -273.15 + chunk2(getValues(crop(tminfiles[[index]], geo_template)), ncores)[[chunk]]
                           tmax <- -273.15 + chunk2(getValues(crop(tmaxfiles[[index]], geo_template)), ncores)[[chunk]]
                         }
                         # frost boundaries to growing season
                         if(doy <= 182){
                           last_frost <- Cond(tmin <= -2, doy, last_frost)
                         }
                         if(doy >= 183){
                           first_frost <- Cond(first_frost == 0, (tmin <= -2) * doy, first_frost)
                         }
                         
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
                                          ncol = length(template) , 
                                          byrow = FALSE)
                         ls_udt <- matrix(stage_udt[as.vector(lifestage)], 
                                          nrow = ncohort,
                                          ncol = length(template), 
                                          byrow = FALSE)
                         ls_dd  <- t(sapply(1:ncohort, FUN = function(x) stage_dd[x, ][lifestage[x, ]]))
                         
                         # Calculate stage-specific degree-days for each cell per day
                         dd_tmp <- TriDD(tmax, tmin, ls_ldt, ls_udt)
                         # Accumulate degree days ----
                         dd_accum <- dd_accum + dd_tmp # resets with lifestage change
                         dd_total <- dd_total + dd_tmp[1, ] # tracks annual accumulation
                         
                         # Calculate lifestage progression: Is accumulation > lifestage requirement
                         progress <- dd_accum >= ls_dd
                         # If reproductive adult stage progressed, that cell has oviposition and generation count increases
                         numgen <- numgen + (progress == 1 & lifestage %in% which(stgorder %in% c("A", "OA"))) 
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
                           # weather data, doesn't really need to be weekly, but easiest to include in same array
                           ls_array[, week, length(stgorder) + 4] <- dd_total
                           ls_array[, week, length(stgorder) + 5] <- last_frost
                           ls_array[, week, length(stgorder) + 6] <- first_frost
                         }
                       }  # daily loop 
                       
                       
                       return(ls_array)
                       
                     } # close foreach loop
}) # system.time

# try to stop cluster, otherwise old rstudio sessions stay on server indefinitely
stopCluster(cl)

# 
# # check results by plotting values assigned to raster
# # array dimensions: [pixel, week, results]
plot(setValues(geo_template, outlist[, 53, 5, 1]))
# 
# # plot time series of lifestage at a single pixel for weighted cohorts
# plot(outlist[50000, , 1])
# 
# 


sites <- data.frame(ID = c("Corvallis, OR", "JB Lewis-McChord, WA", 
                            "Yakima Training Center, WA", "Camp Rilea, OR",
                            "S Portland, OR",
                           "Sutherlin, OR", "Bellingham, WA"),
                    x = c(-123.263, -122.53, -120.461073,
                          -123.934759, -122.658887,
                          -123.315854, -122.479482),
                    y = c(44.564, 47.112, 46.680138,
                          46.122867, 45.470532,
                          43.387721, 48.756105))
coordinates(sites) <- ~x+y

# get cells from geo_template to extract
# We store the cell number
id.cell <- extract(geo_template, sites, cellnumbers=TRUE)[, 1]

# site results from giant array
site_array <- outlist[id.cell,,,]
saveRDS(site_array, "site_array_apha.rds")
