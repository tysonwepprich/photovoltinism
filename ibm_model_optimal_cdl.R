# Speeding up raster-based model with arrays (cohort x pixel)
# Universal code for 3 biocontrol species
# Any number of lifestages allowed
# IBM version


# 1. Setup -------
# packages, options, functions loaded
# TODO: find ways to take out dplyr, purrr, mixsmsn functions?
pkgs <- c("sp", "rgdal", "rgeos", "raster", "lubridate", "mixsmsn", "dplyr", "daymetr",
          "stringr", "prism", "purrr", "foreach", "doParallel", "abind", "ncdf4", "ggplot2")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them

# load collection of functions for this model
source('CDL_funcs.R')
source('species_params.R')

# 2. User input -----

# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model scope
yr           <- 2016
start_doy    <- 1
end_doy      <- 365
region_param <- "NW_SMALL" # TEST/WEST/EAST/CONUS/SOUTHWEST/NORTHWEST
species      <- "GCA" # GCA/APHA/DCA
biotype      <- "S" # TODO: add options for each species, N or S for APHA and GCA

nsim <- 500
# lambda <- 1.1
lat <- 43.39
lon <- -123.35
startyear <- 1993
endyear <- 2018
site <- "Sutherlin"

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
# coordinates(sites) <- ~x+y
REGION <- assign_extent(region_param = region_param)

# Sutherlin: 43.39, -123.35
# Yakima: 46.75, -119.98
# Bellingham: 48.72, -122.48
# Palermo: 39.41, -121.58
# Baskett: 44.98, -123.27
# McArthur: 41.10, -121.41
# Montesano: 46.95, -123.64

# photoperiod decision inclusion
# 2 for logistic, 1 for single value CDL, 0 for none
model_CDL  <- 1


# Derived parameters -----
params <- species_params(mod_type = "ibm", species, biotype, nsim, model_CDL, dd_sd = 5)

ldt <- params$stage_ldt[1]
udt <- params$stage_udt[1]

# Download weather data ----
# # janky because two sites in water without daymet data
# gddsites <- readRDS("lythrum.rds")
# # siteslist <- list()
# for(i in 202:nrow(gddsites)){
#   temp <- download_daymet(site = as.character(i), lat = gddsites$Latitude[i], lon = gddsites$Longitude[i],
#                           start = startyear, end = endyear, internal = TRUE,
#                           silent = TRUE, force = FALSE)
#   outdf <- temp$data %>% 
#     mutate(elev = temp$altitude)
#   siteslist[[i-2]] <- outdf
# }
# gdd <- bind_rows(siteslist) %>% 
#   mutate(Site = rep(1:746, each=9490),
#          Latitude = gddsites$Latitude[-c(181,201)][Site],
#          Longitude = gddsites$Longitude[-c(181,201)][Site])
# sites <- 1:746

# for dataframe of sites
siteslist <- list()
for(i in 1:nrow(sites)){
  temp <- download_daymet(site = sites$ID[i], lat = sites$y[i], lon = sites$x[i],
                          start = startyear, end = endyear, internal = TRUE,
                          silent = TRUE, force = FALSE)
  outdf <- temp$data %>%
    mutate(elev = temp$altitude,
           Site = sites$ID[i],
           Latitude = sites$y[i],
           Longitude = sites$x[i])
  siteslist[[i]] <- outdf
}
gdd <- bind_rows(siteslist)


# saveRDS(gdd, "lythrum_gdd.rds")
# Calculate degree days and photoperiod ----
g1 <- readRDS("ngermany.rds")
g2 <- readRDS("sgermany.rds")
gdd_all <- bind_rows(g1,g2) %>% 
  rename(Year = YEAR, yday = YDAY)


gdd <- readRDS("lythrum_gdd.rds")

gdd_all <- gdd %>%
  filter(yday <= end_doy & yday >= start_doy) %>%
  rowwise() %>%
  mutate(degday = ifelse(yday < 210, TriDD(tmax..deg.c., tmin..deg.c., ldt, udt), 0)) %>%
  ungroup() %>%
  group_by(Site, year) %>%
  arrange(yday) %>%
  mutate(accumdegday = cumsum(degday),
         daylength = photoperiod(Latitude, yday),
         daylength_daymet = dayl..s. / 60 / 60,
         lightfrost = tmin..deg.c. <= 0,
         hardfrost = tmin..deg.c. <= -2)

# naive max voltinism by gdd alone
ow_dd <- params$stage_dd[,1]
gen_dd <- rowSums(params$stage_dd[,-1])
maxvolt <- gdd_all %>% 
  filter(yday == 365) %>% 
  summarise(maxgen = mean(floor((accumdegday - ow_dd) / gen_dd)))



# # MACA
# # MACAv2 files have 5 years in one netcdf file
# weather_path <- "data/maca"
# 
# # if (weather_data_source == "macav2"){
#   yr_path <- paste(weather_path, yr, sep = "/")
#   tminbrick <- brick(paste(yr_path, list.files(yr_path)[grep(x = list.files(yr_path), pattern = "tasmin", fixed = TRUE)], sep = "/"), 
#                      varname = "air_temperature", lvar = 3, level = 4)
#   tmaxbrick <- brick(paste(yr_path, list.files(yr_path)[grep(x = list.files(yr_path), pattern = "tasmax", fixed = TRUE)], sep = "/"),
#                      varname = "air_temperature", lvar = 3, level = 4)
#   
#   macafile <- basename(tmaxbrick@file@name)
#   yrs <- gregexpr(pattern = "[0-9]{4}", text = macafile)
#   yr1 <- as.numeric(substr(macafile, start = yrs[[1]][1], stop = yrs[[1]][1] + 3))
#   yr2 <- as.numeric(substr(macafile, start = yrs[[1]][2], stop = yrs[[1]][2] + 3))
#   yr <- yr1:yr2
#   
#   tminbrick <- shift(tminbrick, x = -360)
#   tmaxbrick <- shift(tmaxbrick, x = -360)
#   
#   # geo_template <- crop(tminbrick[[1]], REGION)
#   # template <- crop(aggregate(raster(files[1]), fact = 2), REGION)
#   # geo_template[!is.na(geo_template)] <- 0
# # }
# 
#   geo_template <- crop(tminbrick[[1]], REGION)
#   tminbrick <- crop(tminbrick, geo_template)
#   tmaxbrick <- crop(tmaxbrick, geo_template)
#   
#   daylist <- list()
#   for (modyear in yr){
#     for(yd in 1:365){
#     tminfiles <- tminbrick[[which(lubridate::year(as.Date(getZ(tminbrick))) == modyear)]][[yd]]
#     tmaxfiles <- tmaxbrick[[which(lubridate::year(as.Date(getZ(tmaxbrick))) == modyear)]][[yd]]
#     
#     tmin <- -273.15 + raster::extract(tminfiles, sites)
#     tmax <- -273.15 + raster::extract(tmaxfiles, sites)
#     
#     outdf <- data.frame(SiteID = sites$ID, lat = sites@coords[,2], lon = sites@coords[,1],
#                         day = yd, year = modyear, tmin = tmin, tmax = tmax)
#     
#     daylist[[length(daylist) + 1]] <- outdf
#     }
#   }
#     
#     gdd_all <- bind_rows(daylist) %>% 
#       mutate(degday = TriDD(tmax, tmin, ldt, udt)) %>% 
#       group_by(year) %>% 
#       arrange(day) %>% 
#       mutate(accumdegday = cumsum(degday),
#              daylength = photoperiod(lat, day),
#              lightfrost = tmin <= 0,
#              hardfrost = tmin <= -2) %>% 
#       mutate(yday = day)
# 
#   saveRDS(gdd_all, "gdd_all_maca.rds")
# gdd_all <- readRDS("gdd_all_maca.rds")

# 3. Run lifestage model -----

# just far western sites
# took 10.8 hours with 40 cores
# gdd_west <- gdd_all %>%
#   filter(Longitude < -120)
gdd_west <- gdd_all %>%
  filter(Longitude > -85 & Longitude < -70)

sites <- unique(gdd_west$Site)
cdls <- expand.grid(cdl = seq(12, 17.5, by = .25),
                    cdl_sd = 0.33, # seq(0, .66, length.out = 3),
                    lambda = 1.5, # seq(0.5, 2, by = 0.5),
                    site = sites)

# cdls <- expand.grid(cdl = NA,
#                     cdl_sd = NA,
#                     lambda = 1)


# startyear <- 2001
# endyear <- 2005
startyear <- min(gdd_all$year)
endyear <- max(gdd_all$year)

ncores <- 40
cl <- makePSOCKcluster(ncores)
registerDoParallel(cl)

test <- system.time({
  
  outlist <- foreach(ncdl = 1:nrow(cdls),
                     .packages = c("raster", "dplyr"),
                     .inorder = FALSE) %:%
    # outlist <- foreach(site = 1:length(sites),
    #                    .packages = "dplyr",
    #                    .inorder = TRUE) %:%
    foreach(y = 1:length(startyear:endyear),
            .packages = c("raster", "dplyr"),
            .inorder = FALSE) %dopar%{
              
              
              
              allyrs <- c(startyear:endyear)
              thisyr <- allyrs[y]
              gdd <- gdd_west %>% 
                filter(year == thisyr & Site == cdls[ncdl, "site"])
              # filter(year == thisyr & SiteID == sites$ID[site])
              
              set.seed(thisyr)
              # Varying traits from parameter file
              stgorder   <- params$stgorder
              stage_dd   <- params$stage_dd
              stage_ldt  <- params$stage_ldt
              stage_udt  <- params$stage_udt
              photo_sens <- params$photo_sens
              # CDL        <- params$CDL[1]
              # cdl_b0     <- params$CDL[2]
              # cdl_b1     <- params$CDL[3]
              
              cdl_mu <- cdls[ncdl, "cdl"]
              cdl_sd <- cdls[ncdl, "cdl_sd"]
              lambda <- cdls[ncdl, "lambda"]
              # lambda <- 1
              
              # vector of which individuals to simulate from data.frame
              ind_df <- data.frame(stage_dd)
              names(ind_df) <- stgorder
              
              ls_complete <- data.frame(array(data = 0, dim = dim(stage_dd)))
              names(ls_complete) <- paste(stgorder, "DOY", sep = "_")
              
              ind_df <- cbind(ind_df, ls_complete)
              ind_df$dd_accum <- 0
              ind_df$diapause <- -1
              ind_df$numgen <- 0
              ind_df$parent <- NA
              ind_df$lifestage <- 1
              # ind_df$cdl <- 0
              ind_df$photo_used <- 0
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
                diap <- newdf$diapause
                photo_used <- newdf$photo_used
                stage_dd <- as.matrix(newdf[, stgorder])
                
                doy <- start_doy + index - 1
                
                # Assign lifestage raster -----
                # Each raster cell is assigned to a lifestage
                # These three rasters assign physiological parameters by cell
                ls_dd  <- sapply(1:ncohort, FUN = function(x) stage_dd[x, ][lifestage[x]])
                
                # TODO: Calculate stage-specific degree-days for each cell per day
                #  Now just has gdd calculated ahead of time with same thresholds
                
                dd_tmp <- rep(gdd$degday[which(gdd$yday == doy)], times = ncohort)
                # Accumulate degree days ----
                dd_accum <- dd_accum + dd_tmp # resets with lifestage change
                
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
                  # photo <- photoperiod(lats, doy, p = 1.5)
                  photo <- gdd$daylength[which(gdd$yday == doy)]
                }
                
                if (model_CDL == 1){
                  ## if CDL is mu and sd
                  CDL <- rnorm(n = ncohort, mean = cdl_mu, sd = cdl_sd)
                  sens_mask <- Cond(lifestage %in% photo_sens & diap == -1, 1, 0)
                  prop_diap <- photo < CDL
                  diap <- Cond(sens_mask == 1, rbinom(ncohort, size = 1, prob = prop_diap), diap)
                  
                  ## if CDL is one number
                  # sens_mask <- Cond(lifestage %in% photo_sens & diap == -1, 1, 0)
                  # prop_diap <- photo < CDL
                  # diap <- Cond(sens_mask == 1, rbinom(ncohort, size = 1, prob = prop_diap), diap)
                }
                
                if (model_CDL == 2){
                  sens_mask <- Cond(lifestage %in% photo_sens & diap == -1, 1, 0)
                  prop_diap <- exp(cdl_b0 + cdl_b1 * photo) /
                    (1 + exp(cdl_b0 + cdl_b1 * photo))
                  diap <- Cond(sens_mask == 1, rbinom(ncohort, size = 1, prob = prop_diap), diap)
                }
                
                newdf$photo_used <- Cond(diap != newdf$diapause, photo, photo_used)
                # newdf$cdl <- Cond(diap != newdf$diapause, 0, CDL)
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
                  # simulate oviposition variation
                  egg_lam <- ifelse(newrows$numgen == 0, 1, lambda)
                  negg <- rpois(n = nrow(newrows), lambda = egg_lam)
                  newrows <- newrows[rep(row.names(newrows), negg), ]
                  
                  if(nrow(newrows) > 0){
                    
                    
                    # Draw degree-day stage requirements here for variation (just from existing params in ind_df for simplicity)
                    newrows[, stgorder] <- apply(X = ind_df[, stgorder], MARGIN = 2,
                                                 FUN = function(x) base::sample(x, size = nrow(newrows), replace = TRUE))
                    
                    newrows[, paste(stgorder, "DOY", sep = "_")] <- 0
                    newrows$dd_accum <- 0
                    newrows$diapause <- -1
                    newrows$numgen <- newrows$numgen + 1
                    newrows$parent <- newrows$indiv
                    newrows$lifestage <- which(stgorder == "E")
                    newrows$photo_used <- 0
                    newrows$indiv <- (max(ind_df$indiv) + 1):(max(ind_df$indiv) + nrow(newrows))
                    newrows$active <- 1
                  }
                  
                  newdf <- rbind(newdf, newrows)
                  newdf$active[which(newdf$indiv %in% ovip)] <- 0
                }
                
                ind_df <- rbind(olddf, newdf)
              }  # daily loop
              
              # Summarize results -----
              voltinism <- ind_df %>% 
                filter(diapause == 1 & active == 0) %>% 
                group_by(numgen) %>% 
                tally() %>% 
                stats::weighted.mean(x = .$numgen, w = .$n)
              
              if(is.nan(voltinism)){ 
                voltinism <- data.frame(voltinism = NA)
              }
              
              attvoltinism <- ind_df %>% 
                filter(diapause == 1 | active == 1) %>% 
                # filter(numgen > 0) %>% 
                group_by(numgen) %>% 
                tally() %>% 
                stats::weighted.mean(x = .$numgen, w = .$n)
              
              volt_comp <- ind_df %>%
                filter(diapause == 1 & active == 0) %>% 
                group_by(numgen) %>% 
                tally() %>% 
                tidyr::spread(numgen, n)
              
              
              if(nrow(volt_comp) == 0){
                volt_comp <- data.frame(comp_1 = NA)
              }else{
                names(volt_comp) <- paste("comp", names(volt_comp), sep = "_")
              }
              
              volt_att <- ind_df %>%
                filter(numgen > 0) %>% 
                group_by(numgen) %>% 
                tally() %>% 
                tidyr::spread(numgen, n)
              
              names(volt_att) <- paste("att", names(volt_att), sep = "_")
              
              
              lost <- ind_df %>% 
                filter(active == 1) %>% 
                group_by(numgen) %>% 
                nrow()
              
              ann_lam <- (ind_df %>% 
                            filter(diapause == 1 & active == 0) %>% 
                            nrow()) / nsim
              
              results <- data.frame(SiteID = cdls[ncdl, "site"], Year = thisyr, cdl_mu, cdl_sd, 
                                    lambda, voltinism, attvoltinism, lost, volt_att, volt_comp, ann_lam)
              
              # results <- data.frame(SiteID = sites$ID[site],
              # Year = thisyr, lambda, voltinism, lost, volt_att, volt_comp, ann_lam)
              
              return(results)
              
            } # close foreach loop
}) # system.time

stopCluster(cl)


# quirks: indiv can decide to diapause at beginning of TA, but still be active because is doesn't make it to complete TA stage before end of year

res <- bind_rows(flatten(outlist)) %>% 
  arrange(Year, SiteID)
saveRDS(res, "ibm_results_nesites.rds")


# Compare to cohort model results

# UNKNOWN: Does voltinism from IBM == voltinism from cohort model?
# Cohort has weighted completed x weighted diapause, is this the same as
# weighted mean of [completed(1) x diapause(1), completed(2) x diapause(2), ...]?
# Go to cohort model and test by putting voltinism calc (1st for each cohort)
# in matrix of weighted results 

coh <- readRDS("site_array_gca_17coh.rds")

# extract gdd, voltinism, lost, 
cohlist <- list()
for (w in 1:53){
  
  coh_att <- data.frame(coh[, w, 7,])
  names(coh_att) <- paste0("X", 2016:2020)
  coh_att <- coh_att %>% 
    mutate(ID = sites$ID) %>% 
    tidyr::gather(year, attempted, X2016:X2020) %>% 
    mutate(year = gsub(pattern = "X", replacement = "", x = year, fixed = TRUE),
           week = w)
  
  
  coh_volt <- data.frame(coh[, w, 8,])
  names(coh_volt) <- paste0("X", 2016:2020)
  coh_volt <- coh_volt %>% 
    mutate(ID = sites$ID) %>% 
    tidyr::gather(year, completed, X2016:X2020) %>% 
    mutate(year = gsub(pattern = "X", replacement = "", x = year, fixed = TRUE),
           week = w)
  
  coh_wvolt <- data.frame(coh[, w, 9,])
  names(coh_wvolt) <- paste0("X", 2016:2020)
  coh_wvolt <- coh_wvolt %>% 
    mutate(ID = sites$ID) %>% 
    tidyr::gather(year, wvolt, X2016:X2020) %>% 
    mutate(year = gsub(pattern = "X", replacement = "", x = year, fixed = TRUE),
           week = w)
  
  coh_diap <- data.frame(coh[, w, 10,])
  names(coh_diap) <- paste0("X", 2016:2020)
  coh_diap <- coh_diap %>% 
    mutate(ID = sites$ID) %>% 
    tidyr::gather(year, diap, X2016:X2020) %>% 
    mutate(year = gsub(pattern = "X", replacement = "", x = year, fixed = TRUE),
           week = w)
  
  coh_all <-  coh_att %>% 
    left_join(coh_volt) %>% 
    left_join(coh_wvolt) %>% 
    left_join(coh_diap)
  cohlist[[w]] <- coh_all
}
cohplot <- bind_rows(cohlist) %>% 
  tidyr::gather(var, val, attempted, completed, wvolt, diap)

# plot helps show that diapause * completed doesn't make sense, if 

ggplot(cohplot, aes(x = week, y = val, group = ID, color = ID)) +
  geom_line() +
  facet_grid(var ~ year, scales = "free") +
  theme_bw(base_size = 14)

# best diapause cutoff to make IBM??
# cohplot <- bind_rows(cohlist)
# ggplot(cohplot, aes(x = week, y = ifelse(diap > .75, NA, 1) * completed, group = ID, color = ID)) +
#   geom_line(aes(size = 1-diap)) +
#   facet_grid(. ~ year) +
#   theme_bw(base_size = 14)
# 
# coh_volt <- cohplot %>% 
#   mutate(d25 = ifelse(diap > .25, NA, 1) * completed,
#          d50 = ifelse(diap > .5, NA, 1) * completed,
#          d75 = ifelse(diap > .75, NA, 1) * completed,
#          d90 = ifelse(diap > .9, NA, 1) * completed,
#          d95 = ifelse(diap > .95, NA, 1) * completed,
#          d98 = ifelse(diap > .98, NA, 1) * completed) %>% 
#   group_by(ID, year) %>% 
#   summarise_all(max, na.rm = TRUE) %>% 
#   arrange(year, ID)
# 
# # good correlation between the cohort and IBM voltinism estimates now
# # but what does it really mean? 
# # is this just weighted voltinism? Why does it work for the cohort model?
# coh_volt <- coh_volt %>% ungroup() %>%  mutate(res$voltinism)
# plot(coh_volt$`res$voltinism`, coh_volt$d95)
# abline(0, 1)

# this matches IBM well if attempted used instead of completed, why?
cohdf <- bind_rows(cohlist) %>% 
  group_by(ID, year) %>% 
  arrange(week) %>% 
  # mutate(diapnewgen = cumsum(c(0, diff(attempted)) * (1-diap))) %>% 
  filter(week == 53) %>% 
  ungroup() %>% 
  rename(SiteID = ID, Year = year) %>% 
  mutate(Year = as.numeric(Year),
         mismatch = wvolt - completed)

ind_diap <- res %>% 
  rowwise() %>% 
  mutate(ind_comp = sum(c(comp_1, comp_2, comp_3, comp_4), na.rm = TRUE),
         ind_diap = ind_comp / (ind_comp + lost)) %>% 
  left_join(cohdf)

plot(ind_diap$voltinism, ind_diap$wvolt)
abline(0, 1)
plot(ind_diap$ind_diap, ind_diap$diap)
abline(0, 1)


# Optimal CDL plots
library(ggplot2)
library(viridis)
###############
res <- readRDS("ibm_results_nwsites.rds") %>% 
  filter(cdl_sd == .33, lambda == 1.5, SiteID == "Sutherlin, OR", Year > 2009)
# filter(cdl_sd == .33, lambda == 1.5, SiteID %in% unique(SiteID)[c(1, 2, 6, 7, 9, 10)])

plt <- ggplot(res, aes(x = cdl_mu, y = attvoltinism)) +
  geom_point(aes(color = ann_lam), size = 3.5) +
  scale_x_reverse() +
  scale_color_viridis(begin = 0, end = 1) +
  facet_wrap(~Year, ncol = 3) +
  # theme_bw() +
  ylab("Voltinism: weighted average of generations attempted") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated voltinism and projected growth rate by critical photoperiod")
plt

ggsave(filename = paste0("GCA_voltinism_", site, ".png"), plot = plt, device = "png", width = 9, height = 6, units = "in")


# heatmap
plt <- ggplot(res, aes(x = cdl_mu, y = lambda, fill = log(ann_lam))) +
  geom_raster() +
  scale_x_reverse() +
  scale_fill_viridis(begin = 1, end = 0) +
  facet_wrap(~Year, ncol = 3) +
  theme_bw() +
  ylab("Between generation lambda") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated log(Annual Population Growth Rate)", subtitle = "Accounting for lost generations")
plt


# lines
plt <- ggplot(res, aes(x = cdl_mu, y = log(ann_lam), group = lambda, color = lambda)) +
  geom_line() +
  scale_x_reverse() +
  scale_color_viridis(begin = 1, end = 0) +
  facet_wrap(~Year, ncol = 3) +
  theme_bw() +
  ylab("Log(annual lambda)") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated log(Annual Population Growth Rate)", subtitle = "Accounting for lost generations")
plt

# across years
res <- readRDS("ibm_results_west.rds") # %>%  filter(SiteID == 123)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geom_lam <- res %>% 
  # filter(Year > 2007) %>% 
  group_by(SiteID, cdl_mu, cdl_sd, lambda) %>% 
  # tidyr::complete(Year, cdl_mu, cdl_sd, lambda, fill = list(ann_lam = 0)) %>% 
  summarise(mean_annual_lambda = gm_mean(ann_lam + 0.001))

plt <- ggplot(geom_lam, aes(x = cdl_mu, y = lambda, fill = log(mean_annual_lambda))) +
  geom_raster() +
  scale_x_reverse() +
  scale_fill_viridis(begin = 1, end = 0) +
  theme_bw() +
  facet_wrap(~cdl_sd) +
  ylab("Between generation lambda") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated mean log(Annual Population Growth Rate)", subtitle = "1993-2017 geometric mean log(lambda)")
plt

# lines
plt <- ggplot(geom_lam, aes(x = cdl_mu, y = log(mean_annual_lambda), color = SiteID, group = SiteID)) +
  geom_line() +
  scale_x_reverse() +
  # scale_color_viridis(begin = 1, end = 0) +
  facet_wrap(~cdl_sd) +
  # theme_bw() +
  ylab("Log(annual lambda)") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated log(Annual Population Growth Rate)", subtitle = "Accounting for lost generations")
plt


# best cdl combo
res <- readRDS("ibm_results_midwest.rds")
geom_lam <- res %>% 
  left_join(maxvolt, by = c("Year" = "year", "SiteID" = "Site")) %>% 
  group_by(SiteID, cdl_mu, cdl_sd, lambda) %>% 
  summarise(mean_annual_lambda = gm_mean(ann_lam + 0.001),
            mean_lost = mean(lost),
            mean_wvolt = mean(ifelse(is.na(voltinism), 0, voltinism)),
            mean_attvolt = mean(attvoltinism),
            mean_mismatch = mean(attvoltinism - maxgen))

# ggplot(geom_lam, aes(x = mean_annual_lambda)) +
#   geom_density() +
#   facet_wrap(~lambda)

# maximize geometric mean of annual lambda (non including overwintering)
best_cdl <- geom_lam %>% 
  filter(lambda == 1.5) %>% 
  group_by(SiteID) %>% 
  mutate(diffbest = mean_annual_lambda - max(mean_annual_lambda)) %>% 
  filter(diffbest == 0) %>% 
  arrange(SiteID, cdl_mu) %>% 
  slice(1) #%>% 
  # filter(cdl_mu > 12.2)

best_wvolt <- geom_lam %>% 
  filter(lambda == 1.5) %>% 
  group_by(SiteID) %>% 
  mutate(diffbest = mean_wvolt - max(mean_wvolt)) %>% 
  filter(diffbest == 0) %>% 
  arrange(SiteID, cdl_mu) %>% 
  slice(1)

best_avolt <- geom_lam %>% 
  filter(lambda == 1.5) %>% 
  group_by(SiteID) %>% 
  mutate(diffbest = mean_attvolt - max(mean_attvolt)) %>% 
  filter(diffbest == 0) %>% 
  arrange(SiteID, cdl_mu) %>% 
  slice(1)

best_lost <- geom_lam %>% 
  filter(lambda == 1.5) %>% 
  group_by(SiteID) %>% 
  mutate(diffbest = mean_lost - min(mean_lost)) %>% 
  filter(diffbest == 0) %>% 
  arrange(SiteID, cdl_mu) %>% 
  slice(1)

best_mm <- geom_lam %>% 
  filter(lambda == 1.5) %>% 
  group_by(SiteID) %>% 
  mutate(diffbest = abs(mean_mismatch) - min(abs(mean_mismatch))) %>% 
  filter(diffbest == 0) %>% 
  arrange(SiteID, cdl_mu) %>% 
  slice(1)

plot(jitter(best_cdl$cdl_mu), jitter(best_lost$cdl_mu), ylab = "CP minimizing lost gen", xlab = "CP maximizing lambda")
abline(0,1)
plot(jitter(best_cdl$cdl_mu), jitter(best_mm$cdl_mu), ylab = "CP minimizing mismatch", xlab = "CP maximizing lambda")
abline(0,1)
# 
# # extract voltinism/lost/attempted info for the row with best cdl
# for (site in unique(res$SiteID)){
#   best <- best_cdl %>% filter(SiteID == site)
#   tmp <- res %>% 
#     filter(SiteID == site, cdl_mu == best$cdl_mu, cdl_sd == best$cdl_sd, lambda == best$lambda)
#   
#   attempts <- tmp %>% 
#     tidyr::gather(gen, num, starts_with("att"))
#   
#   
# }


# site variables from gdd data
site_gdd <- gdd_all %>% 
  group_by(Site) %>% 
  summarise(mean_accumdd = mean(accumdegday[yday == 365]),
            early_accumdd = mean(accumdegday[yday == 170]),
            late_accumdd = mean(accumdegday[yday == 365] - accumdegday[yday == 170]),
            percgdd_early = early_accumdd / late_accumdd,
            lat = Latitude[1],
            lon = Longitude[1])

# map
dat <- left_join(best_cdl, site_gdd, by = c("SiteID" = "Site"))

# # gdd better predictor of cdl than latitude
summary(lm(cdl_mu ~ scale(mean_accumdd), data = dat))
summary(lm(cdl_mu ~ scale(lat), data = dat))

region <- rgdal::readOGR("./src/ref/ne_50m_admin_1_states_provinces_lakes", 'ne_50m_admin_1_states_provinces_lakes', encoding='UTF-8')
# region <-spTransform(region, CRS(proj4string(template)))
reg.points = fortify(region, region="name_en")
reg.df = left_join(reg.points, region@data, by = c("id" = "name_en"))
theme_set(theme_bw(base_size = 20))

xs <- range(dat$lon) + c(-.5, .5)
ys <- range(dat$lat) + c(-.5, .5)

tmpplt <- ggplot() +
  geom_polygon(data = reg.df, aes(x = long, y = lat, group = group), fill = NA, color = "light gray", inherit.aes = FALSE, size = 1, alpha = .3) +
  # geom_point(data = dat, aes(x = lon, y = lat, color = as.factor(round(mean_attvolt))), alpha = 1, size = 3, inherit.aes = FALSE) +
  # scale_color_viridis(discrete = TRUE, option = "C") +
  geom_point(data = dat, aes(x = lon, y = lat, color = mean_attvolt), alpha = 1, size = 4, inherit.aes = FALSE) +
  scale_color_viridis(discrete = FALSE, name = "Attempted") +
  coord_fixed(1.3, xlim = xs, ylim = ys, expand = FALSE, clip = "on") +
  # coord_fixed(1.3, xlim = c(-126.5, -119.5), ylim = c(36, 52), expand = FALSE, clip = "on") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
tmpplt



# adjust for lambda
lamlist <- list()
lamrange <- seq(.9, 3, length.out = 25)
for (i in 1:nrow(res)){
  tmp <- res[i, ]
  
  comp <- tmp %>% 
    dplyr::select(contains("comp")) %>% 
    tidyr::gather(gen, num) %>% 
    mutate(gen = as.numeric(str_split_fixed(gen, pattern = coll("_"), 2)[,2]))
  
  for (j in 1:length(lamrange)){
    lam <- comp %>% 
      filter(complete.cases(.)) %>% 
      mutate(lambda = (num * (lamrange[j]) ^ gen))
    
    lamlist[[(length(lamlist) + 1)]] <- tmp %>% 
      dplyr::select(Year, cdl_mu, cdl_sd, voltinism, lost) %>% 
      mutate(gen_lambda = lamrange[j],
             annual_lambda = sum(lam$lambda) / 100)
  }
}


lams <- bind_rows(lamlist)

plt <- ggplot(lams, aes(x = CDL, y = gen_lambda, fill = log(annual_lambda))) +
  geom_raster() +
  scale_x_reverse() +
  scale_fill_viridis(begin = 1, end = 0) +
  facet_wrap(~Year, ncol = 5) +
  theme_bw() +
  ylab("Between generation lambda") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated log(Annual Population Growth Rate)", subtitle = "Accounting for lost generations")
plt
ggsave(filename = paste0("GCA_AnnualLambda_", site, ".png"), plot = plt, device = "png", width = 9, height = 6, units = "in")



gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geom_lam <- lams %>% 
  # filter(Year > 2007) %>% 
  group_by(CDL, gen_lambda) %>% 
  tidyr::complete(Year, CDL, gen_lambda, fill = list(annual_lambda = 0)) %>% 
  summarise(mean_annual_lambda = gm_mean(annual_lambda))

plt <- ggplot(geom_lam, aes(x = CDL, y = gen_lambda, fill = log(mean_annual_lambda))) +
  geom_raster() +
  scale_x_reverse() +
  scale_fill_viridis(begin = 1, end = 0) +
  theme_bw() +
  ylab("Between generation lambda") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated mean log(Annual Population Growth Rate)", subtitle = "1993-2017 geometric mean log(lambda)")
plt

ggsave(filename = paste0("GCA_meanlambda_", site, ".png"), plot = plt, device = "png", width = 9, height = 6, units = "in")





