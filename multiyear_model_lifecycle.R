# Speeding up raster-based model with arrays (cohort x pixel)
# Universal code for 3 biocontrol species
# Any number of lifestages allowed
# Multiple years run at once, 5 yr climate forecasts for example

# TODO:
# 1. Download stable PRISM files to GRUB, many months/years not updated
# 2. MACAV2 option for forecasts
# 3. Redo mapping functions to use new results format
# 4. Add CDL values to foreach loop with goal to find optimal photoperiod response for each pixel
# 5. Crop netcdf on command line with CDO: cdo timmean -sellonlatbox,lon1,lon2,lat1,lat2 -seldate,date1,date2 in.nc out.nc



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
# weather_path <- "data/maca/2016"
weather_path <- "data/PRISM" # PRISM data on grub server (needs to have stable files downloaded)
download_daily_weather <- 0 # 1 if you need to download PRISM/Daymet data first (20 minutes)
# weather_data_source <- "macav2" # could also have daymet or macav2
weather_data_source <- "prism"

# use parallel processing
# splits raster into chunks for running through simulation
runparallel <- 1 # 1 for yes, 0 for no
ncores <- 3 # choose number of cores for parallel
nchunk <- 1 # not used now, for splitting map if too large for memory

# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model scope
yr           <- 2016 # if multi-year, this is start year of 5-year period
start_doy    <- 1
end_doy      <- 365
region_param <- "NW_SMALL" # TEST/WEST/EAST/CONUS/SOUTHWEST/NORTHWEST
# if CONUS, wouldn't need to crop PRISM or MACAv2
species      <- "APHA" # GCA/APHA/DCA
biotype      <- "N" # TODO: add options for each species, N or S for APHA and GCA

if(weather_data_source == "macav2"){
  # in Kelvin, adjust thresholds/frost or add 273.15 to giant rasters?
  weather_files <- list.files(path = weather_path, full.names = TRUE, recursive = TRUE)
  weather_files <- weather_files[grep(pattern = yr, x = weather_files, fixed = TRUE)]
  weather_files <- weather_files[grep(pattern = "tas", x = weather_files, fixed = TRUE)]
}

# introducing cohort variation
# assign to 1 to match previous model versions
ncohort <- 11 # number of cohorts to approximate emergence distribution

# photoperiod decision inclusion
# 2 for logistic, 1 for single value CDL, 0 for none
model_CDL  <- 2 

# Derived parameters -----
REGION <- assign_extent(region_param = region_param)
params <- species_params(mod_type = "cohort", species, biotype, nsim = ncohort, model_CDL, dd_sd = 0)

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
  # yr_path <- paste(weather_path, yr, sep = "/")
  # tmin <- brick(paste(yr_path, list.files(yr_path)[grep(x = list.files(yr_path), pattern = "tasmin", fixed = TRUE)], sep = "/"), 
  #               varname = "air_temperature", lvar = 3, level = 4)
  # tmax <- brick(paste(yr_path, list.files(yr_path)[grep(x = list.files(yr_path), pattern = "tasmax", fixed = TRUE)], sep = "/"),
  #               varname = "air_temperature", lvar = 3, level = 4)
  
  tmax <- brick(weather_files[grep(x = weather_files, pattern = "tasmax", fixed = TRUE)],
                varname = "air_temperature", lvar = 3, level = 4)
  tmin <- brick(weather_files[grep(x = weather_files, pattern = "tasmin", fixed = TRUE)],
                varname = "air_temperature", lvar = 3, level = 4)
  
  macafile <- basename(tmax@file@name)
  yrs <- gregexpr(pattern = "[0-9]{4}", text = macafile)
  yr1 <- as.numeric(substr(macafile, start = yrs[[1]][1], stop = yrs[[1]][1] + 3))
  yr2 <- as.numeric(substr(macafile, start = yrs[[1]][2], stop = yrs[[1]][2] + 3))
  yr <- yr1:yr2
  
  tmin_all <- shift(tmin, x = -360)
  tmax_all <- shift(tmax, x = -360)
  
  geo_template <- crop(tmin_all[[1]], REGION)
  # template <- crop(aggregate(raster(files[1]), fact = 2), REGION)
  geo_template[!is.na(geo_template)] <- 0
}


if ( runparallel == 1){
  # # avoid memory issues on laptop
  # if(ncores > (parallel::detectCores() / 2)){
  #   ncores <- parallel::detectCores() / 2
  # }
  
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
chunk_list <- chunk2(getValues(geo_template), nchunk) # replace with ncores if using
chunk <- 1

# function needed to combine array results with foreach
acomb4 <- function(...) abind::abind(..., along = 4)
acomb1 <- function(...) abind::abind(..., along = 1)


# # FOR PRISM
# # Crop rasters and convert them to a matrix, using multiple cores 
# # .inorder must be set to TRUE so that output files are in correct order!
# tmax_list <- foreach(t=tmaxfiles, .packages = "raster", .inorder = TRUE) %dopar% {
#   m <- as.matrix(crop(raster(t),template))
# }
# 
# tmin_list <- foreach(t=tminfiles, .packages = "raster", .inorder = TRUE) %dopar% {
#   m <- as.matrix(crop(raster(t),template))
# }




# 3. Run lifestage model -----

yr <- c(2016:2018)
test <- system.time({
  
  # parallel by yr if MACA 5-year dataset
  outlist <- foreach(numyrs = 1:length(yr),
                     .packages = c("raster", "lubridate", "stringr", "purrr", "dplyr"), 
                     .inorder = FALSE,
                     .combine = 'acomb4',
                     .multicombine = TRUE) %dopar%{
                       # foreach(chunk = 1:nchunk,
                       # .packages = c("raster", "lubridate"),
                       #         .inorder = TRUE,
                       #         .combine = 'acomb1',
                       #         .multicombine = TRUE) %dopar%{
                       
                       
                       
                       if (weather_data_source == "macav2"){
                         modyear <- yr[numyrs]
                         tminfiles <- tmin_all[[which(lubridate::year(as.Date(getZ(tmin_all))) == modyear)]]
                         tmaxfiles <- tmax_all[[which(lubridate::year(as.Date(getZ(tmax_all))) == modyear)]]
                       }
                       
                       if (weather_data_source == "prism"){
                         modyear <- yr[numyrs]
                         yr_path <- paste(weather_path, modyear, sep = "/")
                         # If GDD not pre-calculated, load list of PRISM files
                         # Add option for GDD pre-calculated (if all stages have same traits)
                         # send file names to for loop
                         pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
                         tminfiles <- list.files(path = yr_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = TRUE)
                         tminfiles <- ExtractBestPRISM(tminfiles, modyear, leap = "keep")[start_doy:end_doy]
                         
                         pattern = paste("(PRISM_tmax_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
                         tmaxfiles <- list.files(path = yr_path, pattern=pattern, all.files=FALSE, full.names=TRUE, recursive = TRUE)
                         tmaxfiles <- ExtractBestPRISM(tmaxfiles, modyear, leap = "keep")[start_doy:end_doy]
                       }
                       
                       template <- chunk_list[[chunk]]
                       
                       template_mat <- matrix(template, 
                                              nrow = ncohort,
                                              ncol = length(template), 
                                              byrow=TRUE)
                       
                       # latitudes for photoperiod model
                       lats <- chunk2(getValues(init(geo_template, 'y')), nchunk)[[chunk]]
                       
                       
                       # create array to hold results
                       # only track weekly results to save disk space/memory
                       # output variables: lifestage, numgen, daily degree-days
                       ls_array <- array(data = NA, dim = c(length(template), 
                                                            length(seq(start_doy, end_doy, by = 7)), 
                                                            length(params$stgorder) + 7)) 
                       
                       # Varying traits from parameter file
                       stgorder   <- params$stgorder
                       relpopsize <- params$relpopsize / sum(params$relpopsize) # so that sum = 1
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
                       voltinism <- template_mat
                       diap <- template_mat
                       newdiap <- template_mat
                       diap_sens <- template_mat + 1
                       
                       # Loop over days ----
                       for (index in 1:length(start_doy:end_doy)) {
                         doy <- start_doy + index - 1
                         # # get daily temperature, crop to REGION
                         if (weather_data_source == "prism"){
                           tmin <- chunk2(getValues(crop(raster(tminfiles[index]), geo_template)), nchunk)[[chunk]]
                           tmax <- chunk2(getValues(crop(raster(tmaxfiles[index]), geo_template)), nchunk)[[chunk]]
                         }
                         if (weather_data_source == "macav2"){
                           tmin <- -273.15 + chunk2(getValues(crop(tminfiles[[index]], geo_template)), nchunk)[[chunk]]
                           tmax <- -273.15 + chunk2(getValues(crop(tmaxfiles[[index]], geo_template)), nchunk)[[chunk]]
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
                                          ncol = length(template), 
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
                         # Problem with voltinism, may decide to diapause before completing generation
                         # need to track who's newly in diapause, but also remove them when complete
                         # don't happen in the same day/week
                         voltinism <- voltinism + newdiap * fullgen * (progress == 1 & lifestage == (length(stgorder) - 1))
                         newdiap <- newdiap - newdiap * (progress == 1 & lifestage == (length(stgorder) - 1))

                        # Take lifestage result from degree-day model and estimate diapause over the year
                         # Need to know photo-sensitive lifestage
                         
                         # photoperiod for this day across vector
                         if (model_CDL != 0){
                           photo <- photoperiod(lats, doy, p = 1.5)
                           photo_mat <- matrix(photo,
                                               nrow = ncohort,
                                               ncol = length(template),
                                               byrow=TRUE)
                           
                         }
                         
                         if (model_CDL == 1){
                           sens_mask <- matrix(Cond(as.vector(lifestage) %in% photo_sens &
                                                      diap_sens == 1, 1, 0),
                                               ncol = length(template), 
                                               byrow = FALSE)
                           prop_diap <- photo < CDL
                           newdiap <- newdiap + sens_mask * prop_diap * (1 - diap)
                           diap <- diap + sens_mask * prop_diap * (1 - diap)
                           diap_sens <- Cond(sens_mask == 1, 0, diap_sens)
                         }
                         
                         if (model_CDL == 2){
                           # TODO: wrap inside function for memory efficiency
                           # add logistic regression variation in CDL response
                           sens_mask <- matrix(Cond(as.vector(lifestage) %in% photo_sens &
                                                      diap_sens == 1, 1, 0),
                                               nrow = ncohort,
                                               ncol = length(template), 
                                               byrow= FALSE) # because as.vector goes by column
                           prop_diap <- exp(cdl_b0 + cdl_b1 * photo_mat) /
                             (1 + exp(cdl_b0 + cdl_b1 * photo_mat))
                           # need to account for prop_diap declining up until solstice
                           # only let proportion diapausing increase over season
                           newdiap <- newdiap + sens_mask * prop_diap * (1 - diap)
                           diap <- diap + sens_mask * prop_diap * (1 - diap)
                           diap_sens <- Cond(sens_mask == 1, 0, diap_sens)
                         }
                         
                         # reset if sensitive to diapause decision when reaches next generation
                         diap_sens <- Cond(progress == 1 & lifestage %in% which(stgorder %in% c("A", "OA")),
                                           1, diap_sens) 
                         
                         # Reset the DDaccum cells to zero (not actually zero!!!) for cells that progressed to next lifestage
                         lifestage <- lifestage + progress
                         dd_accum <- dd_accum - (progress * ls_dd)
                    
                         # Reassign cells that progressed past end of stgorder to first non-overwintering stage
                         lifestage <- Cond(lifestage == (length(stgorder) + 1), 2, lifestage)
                         
                         # weights each cohort now with relpopsize and colSums
                         if (index %in% seq(start_doy, end_doy, by = 7)){
                           week <- which(seq(start_doy, end_doy, by = 7) == index)
                           for (i in 1:length(stgorder)){
                             ls_array[, week, i] <- colSums(relpopsize * (lifestage == i), na.rm = FALSE)
                           }
                           ls_array[, week, length(stgorder) + 1] <-  colSums(relpopsize * numgen, na.rm = FALSE) # generations by oviposition
                           ls_array[, week, length(stgorder) + 2] <-  colSums(relpopsize * fullgen, na.rm = FALSE) # max generations potentially completed
                           ls_array[, week, length(stgorder) + 3] <-  colSums(relpopsize * voltinism, na.rm = FALSE) # max generations potentially completed
                           ls_array[, week, length(stgorder) + 4] <-  colSums(relpopsize * diap, na.rm = FALSE)
                           # weather data, doesn't really need to be weekly, but easiest to include in same array
                           ls_array[, week, length(stgorder) + 5] <- dd_total
                           ls_array[, week, length(stgorder) + 6] <- last_frost
                           ls_array[, week, length(stgorder) + 7] <- first_frost
                           
                           # print(round(ls_array[10000, week, 6:9], 2))
                         }
                       }  # daily loop 
                       
                       
                       return(ls_array)
                       
                     } # close foreach loop
}) # system.time

# try to stop cluster, otherwise old rstudio sessions stay on server indefinitely
stopCluster(cl)

# extract site data for testing against IBM

sites <- data.frame(ID = c("Corvallis, OR", "JB Lewis-McChord, WA", 
                           "Yakima Training Center, WA", "Camp Rilea, OR",
                           "S Portland, OR",
                           "Sutherlin, OR", "Bellingham, WA",
                           "McArthur, CA", "Palermo, CA", "Baskett Slough, OR"),
                    x = c(-123.263, -122.53, -120.461073,
                          -123.934759, -122.658887,
                          -123.315854, -122.479482,
                          -121.41, -121.58, -123.27),
                    y = c(44.564, 47.112, 46.680138,
                          46.122867, 45.470532,
                          43.387721, 48.756105,
                          41.10, 39.41, 44.98))
coordinates(sites) <- ~x+y

# get cells from geo_template to extract
# We store the cell number
id.cell <- extract(geo_template, sites, cellnumbers=TRUE)[, 1]

# site results from giant array
site_array <- outlist[id.cell,,,]
saveRDS(site_array, "site_array_gca_17coh.rds")




# check results by plotting values assigned to raster
# array dimensions: [pixel, week, results]
plot(setValues(geo_template, outlist[, 50, 8, 1]))
plot(setValues(geo_template, outlist[, 50, 9, 1] - outlist[,50,8,1]))

plot(setValues(geo_template, ls_array[, 25, 8]))

# plot time series of lifestage at a single pixel for weighted cohorts
plot(outlist[50000, , 1])

saveRDS(outlist, "APHA_output_2016to2018_NWSMALL_N.rds")

# trying out voltinism quantification

da <- readRDS("test_MACA_output_2016.rds")

tsdat <- outlist[5000, ,  ,1]
tsdat <- ls_array[5000,,]
tsdat <- data.frame(tsdat)
names(tsdat) <- c(params$stgorder, "attempt", "compl", "volt", "diap", "gdd", "lf", "ff")
df <- tsdat %>% 
  dplyr::select(attempt, compl, volt, diap) %>% 
  mutate(week = 1:nrow(.)) %>% 
  # group_by(week) %>% 
  arrange(week) %>% 
  # mutate(compl = compl / sum(params$relpopsize),
  #        diap = diap / sum(params$relpopsize)) %>% 
  mutate(diapnewgen = cumsum(c(0, diff(attempt)) * (1-diap))) %>% 
  tidyr::gather(key = "var", value = "val", -week)




ggplot(df, aes(x = week, y = val, group = var, color = var)) + 
  geom_line() +
  facet_wrap(~var, scales = "free_y", ncol = 1)



br <- brick(replicate(53, geo_template))

attem <- setValues(br, ls_array[, , 6])
poten <- setValues(br, ls_array[, , 7])
weivolt <- setValues(br, ls_array[, , 8])
diap <- setValues(br, ls_array[,, 9])
VoltDiap <- function(x, y){
  ss <- (subset(x, 2:nlayers(x)) - subset(x, 1:(nlayers(x)-1))) * (1 - y)[[2:nlayers(x)]]
  tmp <- calc(ss, fun = cumsum)
  return(tmp)
}
voltdiap <- VoltDiap(attem, diap)
plot(voltdiap[[52]])
# plot(setValues(geo_template, voltdiap[[52]]))
plot(poten[[53]])
plot(weivolt[[53]])
plot(voltdiap[[52]] - poten[[53]])
which.max(poten[[53]])
