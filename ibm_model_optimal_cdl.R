# Speeding up raster-based model with arrays (cohort x pixel)
# Universal code for 3 biocontrol species
# Any number of lifestages allowed
# IBM version


# 1. Setup -------
# packages, options, functions loaded
# TODO: find ways to take out dplyr, purrr, mixsmsn functions?
pkgs <- c("sp", "rgdal", "rgeos", "raster", "lubridate", "mixsmsn", "dplyr", "daymetr",
          "stringr", "prism", "purrr", "foreach", "doParallel", "abind", "ncdf4", "ggplot2", "tidyr")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them

# load collection of functions for this model
source('CDL_funcs.R')
source('species_params.R')

# 2. User input -----

# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model scope
yr           <- 2018
start_doy    <- 1
end_doy      <- 365
region_param <- "NW_SMALL" # TEST/WEST/EAST/CONUS/SOUTHWEST/NORTHWEST
species      <- "GCA" # GCA/APHA/DCA
biotype      <- "S" # TODO: add options for each species, N or S for APHA and GCA

nsim <- 1000
# lambda <- 1.1
lat <- 43.39
lon <- -123.35
startyear <- 1994
endyear <- 2018
# site <- "Sutherlin"
# sites <- read.csv("data/GCA_modeling_sites.csv", header = TRUE)[8:9, ]
# 
# sites <- data.frame(ID = c("Corvallis, OR", "JB Lewis-McChord, WA", 
#                            "Yakima Training Center, WA", "Camp Rilea, OR",
#                            "S Portland, OR",
#                            "Sutherlin, OR", "Bellingham, WA",
#                            "McArthur, CA", "Palermo, CA", "Baskett Slough, OR"),
#                     x = c(-123.263, -122.53, -120.461073,
#                           -123.934759, -122.658887,
#                           -123.315854, -122.479482,
#                           -121.41, -121.58, -123.27),
#                     y = c(44.564, 47.112, 46.680138,
#                           46.122867, 45.470532,
#                           43.387721, 48.756105,
#                           41.10, 39.41, 44.98))

sites <- data.frame(ID = c(
  "Ephrata, WA", "Yakima TC, WA",
  "Sutherlin, OR", "Bellingham, WA",
  "McArthur, CA", "Palermo, CA", "Rickreall, OR"),
  x = c(-119.655253, -119.986530,
           -123.315854, -122.479482,
           -121.41, -121.58, -123.27),
  y = c(47.161647, 46.756318,
          43.387721, 48.756105,
          41.10, 39.41, 44.98))

# sites <- sites[c(3, 6, 7, 8, 9, 10), ]
# sites <- sites[c(4, 6, 7, 8, 9, 22), ]

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
gddsites <- sites
siteslist <- list()
for(i in 1:nrow(gddsites)){
  temp <- download_daymet(site = gddsites$ID[i], lat = gddsites$y[i], lon = gddsites$x[i],
                          start = startyear, end = endyear, internal = TRUE,
                          silent = TRUE, force = FALSE)
  # temp <- download_daymet(site = as.character(i), lat = gddsites$Latitude[i], lon = gddsites$Longitude[i],
  #                         start = startyear, end = endyear, internal = TRUE,
  #                         silent = TRUE, force = FALSE)
  outdf <- temp$data %>%
    mutate(elev = temp$altitude,
           Site = gddsites$ID[i],
           Latitude = gddsites$y[i],
           Longitude = gddsites$x[i])
  siteslist[[i]] <- outdf
}

gdd <- bind_rows(siteslist)
# gdd <- bind_rows(siteslist) %>% 
#   mutate(Site = rep(1:746, each=9490),
#          Latitude = gddsites$Latitude[-c(181,201)][Site],
#          Longitude = gddsites$Longitude[-c(181,201)][Site])
# sites <- 1:746

# gdd <- bind_rows(siteslist) %>% 
#   mutate(Site = gddsites$ID,
#                   Latitude = gddsites$y,
#                   Longitude = gddsites$x)
# saveRDS(gdd, "lythrum_gdd.rds")

# summarise sites
sitesumm <- gdd %>% 
  group_by(Site) %>% 
  summarise(meanTmax = mean(tmax..deg.c.),
            meanTmin = mean(tmin..deg.c.),
            meanPrec = mean(prcp..mm.day.) * 365 / 10,
            elevation = elev[1],
            latitude = Latitude[1],
            longitude = Longitude[1]) %>% 
  arrange(-latitude)
write.csv(sitesumm, "site_summ.csv")

  
# Calculate degree days and photoperiod ----
# g1 <- readRDS("ngermany.rds") %>% 
#   mutate(latitude = 54.4)
# g2 <- readRDS("sgermany.rds") %>% 
#   mutate(latitude = 50.2)
# gdd_west <- bind_rows(g1,g2) %>% 
#   rename(year = YEAR, yday = YDAY) %>% 
#   tidyr::complete(Site, year, yday, fill = list(tmax = 0, tmin = 0, degday = 0)) %>% 
#   mutate(daylength = photoperiod(latitude, yday)) %>% 
#   group_by(Site, year) %>%
#   arrange(yday) %>%
#   mutate(accumdegday = cumsum(degday))
# 
# 
# gdd <- readRDS("lythrum_gdd.rds")

gdd_all <- gdd %>%
  filter(yday <= end_doy & yday >= start_doy) %>%
  rowwise() %>%
  mutate(degday = TriDD(tmax..deg.c., tmin..deg.c., ldt, udt)) %>%
  ungroup() %>%
  group_by(Site, year) %>%
  arrange(yday) %>%
  mutate(accumdegday = cumsum(degday),
         daylength = photoperiod(Latitude, yday, p = 1.5),
         daylength_daymet = dayl..s. / 60 / 60,
         lightfrost = tmin..deg.c. <= 0,
         hardfrost = tmin..deg.c. <= -2)

gdd_all$Site <- plyr::revalue(gdd_all$Site, c("Yakima Training Center, WA"="Vantage, WA", 
                                              "Baskett Slough, OR"="Rickreall, OR"))
gdd_all$Site <- factor(gdd_all$Site, levels = c("Bellingham, WA", "Vantage, WA", "Rickreall, OR",
                                                "Sutherlin, OR", "McArthur, CA", "Palermo, CA"))

# season length
gdd_seas <- gdd_all %>% 
  group_by(Site, year) %>% 
  summarise(dd = max(accumdegday),
            lastfdd = accumdegday[max(which(lightfrost == TRUE & yday < 170))],
            firstfdd = accumdegday[min(which(lightfrost == TRUE & yday > 170))],
            lastfday = yday[max(which(lightfrost == TRUE & yday < 170))],
            firstfday = yday[min(which(lightfrost == TRUE & yday > 170))],            
            dd100 = yday[which.min(abs(accumdegday - 100))],
            dd497 = yday[which.min(abs(accumdegday - 497))],
            dl100 = daylength[which.min(abs(accumdegday - 100))],
            dl497 = daylength[which.min(abs(accumdegday - 497))],
            seaslength = firstfday - lastfday,
            seaslengthdd = firstfdd - lastfdd) %>% 
  group_by(Site) %>% 
  summarise_all(.funs = list(mean, sd), na.rm = TRUE)
  
write.csv(gdd_seas, "gdd_seas.csv")  

tapply(gdd_all$daylength, gdd_all$Site, FUN = max)  


# naive max voltinism by gdd alone
ow_dd <- params$stage_dd[,1]
gen_dd <- rowSums(params$stage_dd[,-1])
# maxvolt <- gdd_west %>% 
  # filter(yday == 365) %>% 
  # summarise(maxgen = mean(floor((accumdegday - ow_dd) / gen_dd)))

maxvolt <- gdd_all %>% 
  filter_all(all_vars(!is.na(.)))  %>% 
  group_by(Site, year) %>% 
  arrange(yday) %>%
  mutate(nday = n(),
    accumdegday = cumsum(degday)) %>% 
  dplyr::filter(yday == max(yday),
         nday > 350) %>% 
  summarise(maxgen = floor((accumdegday - mean(ow_dd)) / mean(gen_dd)))

# Voltinism by degree-days (constant)
# example years 2011,13,15 for cool to warm
# G&C2015 uses 100dd OW ovip, 497 1st gen, 523 subsequent adult emergence

expand_df <- function(df){
  newdf <- data.frame(gen = 0:df$maxgen) %>% 
    mutate(ddreq = median(ow_dd) + gen * median(gen_dd)) %>% 
    rowwise() %>% 
    mutate(accumdegday = gdd_all %>% 
             filter(Site == df$Site[1] & year == df$year[1]) %>% 
             filter(accumdegday - ddreq[1] >= 0) %>% 
             arrange(yday) %>% 
             slice(1) %>% pull(accumdegday),
           yday = gdd_all %>% 
             filter(Site == df$Site[1] & year == df$year[1]) %>% 
             filter(accumdegday - ddreq[1] >= 0) %>% 
             arrange(yday) %>% 
             slice(1) %>% pull(yday))
  return(newdf)
}

milestones <- maxvolt %>% 
  filter(year %in% c(2011, 2013, 2015)) %>% 
  group_by(Site, year) %>% 
  do(expand_df(.))
milestones$year <- as.character(milestones$year)
milestones$date <- as.Date(milestones$yday, origin=as.Date("2015-12-31"))




dds <- gdd_all %>% filter(year %in% c(2011, 2013, 2015)) %>% 
                            mutate(date = as.Date(yday, origin=as.Date("2015-12-31")))
dds$year <- as.character(dds$year)

glab <- milestones %>% 
  group_by(Site, gen) %>% 
  summarise(accumdegday = mean(accumdegday),
            date = min(dds$date) + 10,
            lab = paste("Gen", gen, sep = " ")[1]) 

theme_set(theme_bw(base_size = 12) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())) 

plt <- ggplot(dds, aes(x = date, y = accumdegday, group = year, color = year)) +
  geom_line(size = 1.5) +
  scale_color_brewer(name = "Year", type = "qual") + 
  geom_point(data = milestones, aes(x = date, y = accumdegday), color = "black") +
  geom_hline(data = glab, aes(yintercept = accumdegday), linetype = "dotted") +
  geom_label(data = glab, aes(x = date, y = accumdegday, label = lab), fill = "white", inherit.aes = FALSE) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = c(as.Date("2016-01-01"), as.Date("2016-12-01"))) +
  facet_wrap(~Site, nrow = 3) +
  theme(legend.position = c(.25, .96),
        legend.direction = "horizontal") +
  xlab("Calendar date") + ylab("Accumulated degree-days")
plt
# saved at 850x700 pixels


# Voltinism by photothermograph (constant)
# add photosensitive time to milestones
cdl <- 15
expand_df <- function(df){
  cycle <- "reproduce"
  newdf <- data.frame(gen = 1, start = median(ow_dd), end = median(ow_dd) + 1 * median(gen_dd),
                      sense = median(ow_dd) + 1 * median(gen_dd) - median(rowSums(params$stage_dd[,c(5:6)]))) %>% 
    mutate(daylength = df %>% filter(accumdegday - sense >= 0) %>% 
             arrange(yday) %>% 
             slice(1) %>% pull(daylength),
           choice = ifelse(daylength < cdl, "diapause", "reproduce"))
  while(cycle == "reproduce"){
    if(newdf$choice[nrow(newdf)] == "reproduce"){
      newrow <- newdf[nrow(newdf),] %>% 
        mutate(gen = gen + 1,
               start = end,
               end = end + median(gen_dd),
               sense = sense + median(gen_dd),
               daylength = ifelse(sense < max(df$accumdegday),
                                  df %>% filter(accumdegday - sense >= 0) %>% 
                 arrange(yday) %>% 
                 slice(1) %>% pull(daylength), NA),
               choice = ifelse(daylength < cdl, "diapause", "reproduce"))
      newdf <- bind_rows(newdf, newrow)
    }else{
      cycle = "diapause"
    }
    if(newdf$end[nrow(newdf)] > max(df$accumdegday)){
      cycle = "lost"
      newdf$choice[nrow(newdf)] <- "lost"
    }
  }
  return(newdf)
}

events <- dds %>% 
  group_by(Site, year) %>% 
  do(expand_df(.)) %>% 
  mutate(ymin = 7 + (2015 - as.numeric(year))/3,
         ymax = ymin + .3,
         shading = ifelse(choice == "lost", .4, 1),
         sense = ifelse(choice == "lost", NA, sense),
         choice = as.character(ifelse(choice == "lost", NA, choice)))

sens_time <- events %>% 
  dplyr::select(Site,year,sense) %>% 
  distinct()

plt <- ggplot(dds, aes(x = accumdegday, y = daylength, group = year, color = year)) +
  geom_line(size = 1.5) +
  scale_color_brewer(name = "Year", type = "qual") + 
  geom_rect(data = events, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, 
                               fill = year, alpha = shading), color = "black", inherit.aes = FALSE) +
  scale_fill_brewer(name = "Year", type = "qual") + 
  geom_point(data = events, aes(x = sense, y = (ymin + ymax)/2, shape = choice), color = "black", size = 2, inherit.aes = FALSE) +
  scale_shape_manual(name = "Choice", values=c(4, 1), na.translate=FALSE)+
  # geom_text(data = glab, aes(x = date, y = accumdegday, label = lab), hjust = .1, inherit.aes = FALSE) +
  facet_wrap(~Site, nrow = 3) +
  geom_hline(aes(yintercept = cdl), linetype = "dashed") +
  geom_vline(data = sens_time, aes(xintercept = sense), linetype = "dotted") +
  # theme(legend.position = c(.25, .96),
        # legend.direction = "horizontal") +
  ylab("Daylength (hours)") + xlab("Accumulated degree-days")
plt





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
#       group_by(SiteID, year) %>% 
#       arrange(day) %>% 
#       mutate(accumdegday = cumsum(degday),
#              daylength = photoperiod(lat, day),
#              lightfrost = tmin <= 0,
#              hardfrost = tmin <= -2) %>% 
#       mutate(yday = day)
# 
#   saveRDS(gdd_all, "gdd_all_maca.rds")
# gdd_all <- readRDS("gdd_all_maca.rds")
# gdd_all <- gdd_all %>% 
#         mutate(degday = TriDD(tmax, tmin, ldt, udt)) %>%
#         group_by(SiteID, year) %>%
#         arrange(day) %>%
#         mutate(accumdegday = cumsum(degday),
#               yday = day)
# 
# ow_dd <- params$stage_dd[,1]
# gen_dd <- rowSums(params$stage_dd[,-1])
# maxvolt <- gdd_all %>% 
#   filter(day == 365) %>% 
#   group_by(SiteID, year) %>% 
#   summarise(maxgen = mean(floor((accumdegday - ow_dd) / gen_dd)))


# 3. Run lifestage model -----

# just far western sites
# took 10.8 hours with 40 cores
# gdd_west <- gdd_all %>%
  # filter(Longitude < -120)
# gdd_midwest <- gdd_all %>%
#   filter(Longitude > -95 & Longitude < -90)

sites <- unique(gdd_all$Site)
# sites <- unique(sites$ID)
cdls <- expand.grid(cdl = seq(14, 18, by = 0.25),
                    cdl_sd = 0.25, # c(0, 0.25, 0.5),
                    lambda = seq(1, 3, by = 1),
                    site = sites,
                    sim = 1)

# cdls <- expand.grid(cdl = NA,
#                     cdl_sd = NA,
#                     lambda = 1)

gdd_west <- gdd_all
startyear <- 2011
endyear <- 2015
# startyear <- min(gdd_west$year)
# endyear <- max(gdd_west$year)

ncores <- 7
cl <- makePSOCKcluster(ncores)
registerDoParallel(cl)

test <- system.time({
  # outlist2 <- list()
  outlist <- foreach(ncdl = 1:nrow(cdls),
                     .packages = c("raster", "dplyr"),
                     .inorder = FALSE) %:%
    # outlist <- foreach(site = 1:length(sites),
    #                    .packages = "dplyr",
    #                    .inorder = TRUE) %:%
    foreach(y = 1:length(startyear:endyear),
            .packages = c("raster", "dplyr"),
            .inorder = FALSE) %dopar%{
              
              print(ncdl)
              print(y)
              
              allyrs <- c(startyear:endyear)
              thisyr <- allyrs[y]
              gdd <- gdd_west %>% 
                filter(year == thisyr & Site == cdls[ncdl, "site"])
              # filter(year == thisyr & SiteID == sites$ID[site])
              # set.seed(cdls[ncdl, "site"] * thisyr)
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
              ind_df$diap_gdd <- 0
              ind_df$numgen <- 0
              ind_df$parent <- NA
              ind_df$lineage <- 1:nrow(ind_df)
              ind_df$lifestage <- 1
              ind_df$cdl <- 0
              ind_df$photo_used <- 0
              ind_df$indiv <- 1:nrow(ind_df)
              ind_df$active <- 1    # this will indicate rows to simulate each day
              
              pop_list <- list()
              
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
                cdl <- newdf$cdl
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
                
                #TODO: This seems like a bug
                olddf$diap_gdd <- olddf$diap_gdd + gdd$degday[which(gdd$yday == doy)]
                
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
                  # CDL <- rnorm(n = ncohort, mean = cdl_mu, sd = cdl_sd)
                  sens_mask <- Cond(lifestage %in% photo_sens & diap == -1, 1, 0)
                  prop_diap <- photo < cdl
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
                    newrows$lineage <- newrows$lineage
                    newrows$lifestage <- which(stgorder == "E")
                    newrows$cdl <- rnorm(n = nrow(newrows), mean = cdl_mu, sd = cdl_sd)
                    newrows$photo_used <- 0
                    newrows$indiv <- (max(ind_df$indiv) + 1):(max(ind_df$indiv) + nrow(newrows))
                    newrows$active <- 1
                  }
                  
                  newdf <- rbind(newdf, newrows)
                  newdf$active[which(newdf$indiv %in% ovip)] <- 0
                }
                
                ind_df <- rbind(olddf, newdf)
                
                # daily population counts to save
                pop_list[[index]] <- ind_df %>% 
                  group_by(numgen, lifestage, active, diapause) %>% 
                  tally() %>% 
                  mutate(doy = index)
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

              # redundant with diap below              
              # volt_comp <- ind_df %>%
              #   filter(diapause == 1 & active == 0) %>% 
              #   group_by(numgen) %>% 
              #   tally() %>% 
              #   tidyr::spread(numgen, n)
              # 
              # 
              # if(nrow(volt_comp) == 0){
              #   volt_comp <- data.frame(comp_1 = NA)
              # }else{
              #   names(volt_comp) <- paste("comp", names(volt_comp), sep = "_")
              # }
              
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
              
              ## Can show which OW adults had lost descendants
              # if(lost > 0){
              # lost_parents <- ind_df %>% 
              #   filter(active == 1) %>% 
              #   pull(unique(parent))
              # lost_parents2 <- ind_df %>% 
              #   filter(indiv %in% lost_parents) %>% 
              #   pull(unique(parent))
              # lost_parents3 <- ind_df %>% 
              #   filter(indiv %in% lost_parents2) %>% 
              #   pull(unique(parent))
              # lost_parents4 <- ind_df %>% 
              #   filter(indiv %in% lost_parents3) %>% 
              #   pull(unique(parent))
              # 
              # hist(ind_df$OA_DOY[which(ind_df$indiv %in% lost_parents4)])
              # hist(ind_df$OA_DOY[which(ind_df$indiv %!in% lost_parents4 & ind_df$OA_DOY > 0)])
              # 
              # }else{
              #   lost_parents = NA
              # }
              
              # diap_gdd <- ind_df %>% 
              #   filter(diapause == 1 & active == 0) %>% 
              #   group_by(numgen) %>% 
              #   summarise(meandiapgdd = mean(diap_gdd)) %>% 
              #   tidyr::spread(numgen, meandiapgdd)
              # names(diap_gdd) <- paste("diapgdd", names(diap_gdd), sep = "_")
              # 
              diap <- ind_df %>%
                filter(diapause == 1 & active == 0) %>% 
                group_by(numgen) %>% 
                tally() %>% 
                tidyr::spread(numgen, n)
              
              if(nrow(diap) == 0){
                diap <- data.frame(diap_1 = NA)
              }else{
                names(diap) <- paste("diap", names(diap), sep = "_")
              }
              
              
              repro <- ind_df %>%
                filter(diapause == 0) %>% 
                group_by(numgen) %>% 
                tally() %>% 
                tidyr::spread(numgen, n)
              
              if(nrow(repro) == 0){
                repro <- data.frame(repro_1 = NA)
              }else{
                names(repro) <- paste("repro", names(repro), sep = "_")
              }
              
              
              ann_lam <- (ind_df %>% 
                            filter(diapause == 1 & active == 0) %>% 
                            nrow()) / nsim
              
              # ann_lam_adj_5000 <- sum(ind_df %>% 
              #               filter(diapause == 1 & active == 0) %>% 
              #               mutate(gdd_penalty = 1 - diap_gdd / 5000) %>% 
              #               pull(gdd_penalty))/ nsim
              #   
              # ann_lam_adj_2500 <- sum(ind_df %>% 
              #                           filter(diapause == 1 & active == 0) %>% 
              #                           mutate(gdd_penalty = 1 - diap_gdd / 2500) %>% 
              #                           pull(gdd_penalty))/ nsim
              
              
              results <- data.frame(SiteID = cdls[ncdl, "site"], Year = thisyr, cdl_mu, cdl_sd, 
                                    lambda, voltinism, attvoltinism, lost, volt_att,  
                                    diap, repro, ann_lam)
              
              # thinking about which lineages are doing better or worse
              
              # oa <- ind_df %>%
              #   filter(numgen == 0) %>% 
              #   dplyr::select(lineage, OA)
              # cp <- ind_df %>% 
              #   filter(numgen > 0) %>% 
              #   group_by(lineage, numgen) %>% 
              #   summarise(meancp = mean(ind_df$cdl[ind_df$indiv %in% parent]))
              # lam <- ind_df %>% 
              #   group_by(lineage) %>% 
              #   summarise(lam = length(which(diapause == 1 & active == 0)))
              # linlost <- ind_df %>% 
              #   group_by(lineage, numgen) %>% 
              #   summarise(linlost = length(which(active == 1))) %>% 
              #   filter(linlost > 0)
              # if(length(unique(linlost$numgen)) == 1){
              #   cp <- cp %>% 
              #     filter(numgen == unique(linlost$numgen)[1] - 1)
              #   linlost <- linlost %>% 
              #     ungroup() %>% 
              #     dplyr::select(-numgen)
              # }
              # lineages <- right_join(oa, cp) %>% 
              #   left_join(lam) %>% 
              #   left_join(linlost)
              # lineages$linlost[is.na(lineages$linlost)] <- 0
                              
              # results <- data.frame(SiteID = sites$ID[site],
              # Year = thisyr, lambda, voltinism, lost, volt_att, volt_comp, ann_lam)
              # outlist2[[length(outlist2)+1]] <- results
              out <- list()
              out$indiv <- ind_df %>% 
                mutate(SiteID = cdls[ncdl, "site"], 
                       Year = thisyr, 
                       cdl_mu, cdl_sd,  lambda)
              out$results <- results
              out$pop <- bind_rows(pop_list) %>% 
                mutate(SiteID = cdls[ncdl, "site"], 
                       Year = thisyr, 
                       cdl_mu, cdl_sd,  lambda)
              return(out)
              
              # 
              # # plot proportion over season
              # ggplot(out$pop, aes(x = doy, y = n, fill = numgen)) +
              #   geom_bar(stat = "identity", width = 1)
              # ggplot(out$pop, aes(x = doy, y = n, fill = interaction(diapause, active))) +
              #   geom_bar(stat = "identity", width = 1) +
              #   facet_wrap(~numgen)
              # ggplot(out$pop, aes(x = doy, y = n, fill = lifestage)) +
              #   geom_bar(stat = "identity", width = 1) +
              #   facet_wrap(~numgen)
              
              
              
            } # close foreach loop
}) # system.time

stopCluster(cl)

saveRDS(outlist, "ibm_outlist_GCA_6sites_daymet.rds")

# quirks: indiv can decide to diapause at beginning of TA, but still be active because is doesn't make it to complete TA stage before end of year
# res <- bind_rows(outlist)
# res <- bind_rows(flatten(outlist)) 
# res <- bind_rows(flatten(outlist)) %>% # if >1 sim per combo
  # group_by(Year, SiteID, cdl_mu, cdl_sd, lambda) %>% 
  # summarise_at(vars(voltinism:diap_2), mean, na.rm = TRUE)
# saveRDS(res, "ibm_results_GCA_6sites_daymet.rds")

# res <- readRDS("ibm_outlist_GCA_6sites_daymet.rds")
inds <- outlist %>% 
  flatten() %>%
  purrr::map(1) %>% dplyr::bind_rows()

firsts <- outlist %>% 
  flatten() %>%
  purrr::map(2) %>% dplyr::bind_rows()
seconds <- outlist %>% 
  flatten() %>%
  purrr::map(3) %>% dplyr::bind_rows()

# firsts$SiteID <- plyr::revalue(firsts$SiteID, c("Yakima Training Center, WA"="Vantage, WA", 
#                                                       "Baskett Slough, OR"="Rickreall, OR"))
# firsts$SiteID <- factor(firsts$SiteID, levels = c("Bellingham, WA", "Vantage, WA", "Rickreall, OR",
#                                                         "Sutherlin, OR", "McArthur, CA", "Palermo, CA"))
# seconds$SiteID <- plyr::revalue(seconds$SiteID, c("Yakima Training Center, WA"="Vantage, WA", 
#                                                 "Baskett Slough, OR"="Rickreall, OR"))
# seconds$SiteID <- factor(seconds$SiteID, levels = c("Bellingham, WA", "Vantage, WA", "Rickreall, OR",
#                                                   "Sutherlin, OR", "McArthur, CA", "Palermo, CA"))


# Plot voltinism photothermograph with ridge plots
inds <- inds %>% 
  filter(Year %in% c(2011, 2013, 2015), 
               cdl_mu == 15,
               lambda == 1)
gdd_yr <- gdd_all %>% 
  filter(year %in% c(2011, 2013, 2015))


# Marginal histograms ----
mean_cdl <- mean(inds$cdl[inds$cdl > 0], na.rm = TRUE)
maxgdd <- max(gdd_yr$accumdegday, na.rm = TRUE)
maxphoto <- max(c(max(gdd_yr$daylength), max(inds$cdl, na.rm = TRUE)))
minphoto <- min(gdd_yr$daylength)

cdl_hist <- conseq %>% # conseq from below
  group_by(SiteID) %>% 
  filter(cdl > 0) %>% 
  mutate(breaks = cut(cdl, breaks=seq(min(cdl),max(cdl),length.out = 20), 
                      labels=seq(min(cdl),max(cdl), length.out = 20)[-1], 
                      include.lowest=TRUE)) %>% 
  mutate(daylength = as.numeric(as.character(breaks))) %>%
  group_by(SiteID, daylength) %>% 
  summarise(n = n()) %>% 
  mutate(accumdegday = -50 + 50 * n/max(n),
         Site = SiteID)

# oa <- inds %>% filter(OA_DOY > 0) %>% 
#   left_join(gdd_all, by = c("Year" = "year", "SiteID" = "Site", "OA_DOY" = "yday"))
# teneral <- inds %>% filter(TA_DOY > 0) %>% 
#   left_join(gdd_all, by = c("Year" = "year", "SiteID" = "Site", "TA_DOY" = "yday"))
# sens_stage <- bind_rows(oa, teneral)
sens_hist <- conseq %>% # conseq from below
  group_by(SiteID) %>% 
  mutate(breaks = cut(accumdegday, breaks = seq(min(accumdegday), max(accumdegday), length.out = 100),
                      labels = seq(min(accumdegday), max(accumdegday), length.out = 100)[-1],
                      include.lowest = TRUE)) %>% 
  mutate(accumdegday = as.numeric(as.character(breaks))) %>% 
  group_by(SiteID, accumdegday) %>% 
  summarise(n = n()) %>% 
  mutate(daylength = minphoto - .5 + .5 * n/max(n),
         Site = SiteID)
sens_lines <- sens_hist %>% filter(accumdegday > 400)

plt <- ggplot(gdd_yr, aes(x = accumdegday, y = daylength, group = as.factor(year), color = as.factor(year))) +
  geom_hline(data = cdl_hist, aes(yintercept = daylength, alpha = n/100000), color = "light grey", size = 2) +
  geom_vline(data = sens_lines, aes(xintercept = accumdegday, alpha = n/100000), color = "light grey", size = 2) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(-50, maxgdd + 100), ylim = c(minphoto - .5, 17.2), expand = FALSE) +
  # scale_color_viridis(discrete = TRUE, name = "Year", option = "D") +
  scale_color_brewer(name = "Year", type = "qual") + 
  scale_y_continuous(breaks = c(9:17)) +
  guides(color = guide_legend(reverse=FALSE)) +
  geom_segment(data=cdl_hist, size=1.5, show.legend=FALSE,
               aes(x=-50, xend=accumdegday, y=daylength, yend=daylength), inherit.aes = FALSE) +
  geom_segment(data=sens_hist, size=2, show.legend=FALSE,
               aes(x=accumdegday, xend=accumdegday, y=minphoto - .5, yend=daylength), inherit.aes = FALSE) +
  # geom_rug(data = results, aes(x = accumdegday), color = "gray", alpha = .1) +
  # geom_rug(data = results, aes(y = cdl), color = "gray", alpha = .1) +
  # ggtitle("Photothermographs with critical photoperiod\nand simulated sensitive stage emergence") +
  xlab("Accumulated degree-days") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Daylength (hours)") +
  facet_wrap(~Site, nrow = 3)
plt




# Consequences ridgelines (stacked bar doesn't look as good here, 
# but don't know how to combine these plots with photothermograph)

uniqsims <- distinct(firsts[,1:5])

outcomelist <- list()
for (nr in 1:nrow(uniqsims)){
  ids <- uniqsims[nr, ]
  temp <- inds %>% filter(Year == ids$Year, 
                          SiteID == ids$SiteID, 
                          cdl_mu == ids$cdl_mu,
                          cdl_sd == ids$cdl_sd,
                          lambda == ids$lambda)
  gdd <- gdd_all %>% filter(year == ids$Year, 
                            Site == ids$SiteID) %>% 
    mutate(doy = yday)
  
  temp_conseq <- temp %>%
    mutate(outcome = case_when(active == 1 ~ "lost",
                               diapause == 1 ~ "diapause",
                               diapause != 1 ~ "fits"))
  
  
  out_lost <- temp_conseq %>% 
    filter(indiv %in% parent[outcome == "lost"])
  oa <- out_lost %>% filter(OA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "OA_DOY" = "yday"))
  teneral <- out_lost %>% filter(TA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "TA_DOY" = "yday"))
  out_lost <- bind_rows(oa, teneral) %>% 
    mutate(offspring = "lost")
  
  out_diap <- temp_conseq %>% 
    filter(indiv %in% parent[outcome == "diapause"])
  oa <- out_diap %>% filter(OA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "OA_DOY" = "yday"))
  teneral <- out_diap %>% filter(TA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "TA_DOY" = "yday"))
  out_diap <- bind_rows(oa, teneral) %>% 
    mutate(offspring = "diapause")
  
  out_fits <- temp_conseq %>% 
    filter(indiv %in% parent[outcome == "fits"])
  oa <- out_fits %>% filter(OA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "OA_DOY" = "yday"))
  teneral <- out_fits %>% filter(TA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "TA_DOY" = "yday"))
  out_fits <- bind_rows(oa, teneral) %>% 
    mutate(offspring = "fits")
  
  outcomes <- bind_rows(out_lost, out_diap, out_fits)
  
  no_off <- anti_join(temp_conseq, outcomes, by = "indiv") %>% 
    filter(outcome == "diapause") %>% 
    filter(TA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "TA_DOY" = "yday")) %>% 
    mutate(offspring = "none, diapause")
  
  outcomes <- bind_rows(outcomes, no_off)
  outcomelist[[length(outcomelist) + 1]] <- outcomes
}

outcomes <- bind_rows(outcomelist)
conseq <- outcomes %>%
  mutate(Consequence = case_when(outcome == "diapause" ~ "Diapauses",
                                 offspring == "lost" ~ "Reproduces, offspring lost",
                                 offspring != "lost" ~ "Reproduces, offspring fits"),
         date = as.Date(doy, origin=as.Date("2015-12-31"))) %>% 
  filter(cdl_mu == 15, 
         lambda == 1,
         Year %in% c(2011, 2013, 2015))

conseq$Consequence <- factor(conseq$Consequence, levels = c("Diapauses", "Reproduces, offspring lost", "Reproduces, offspring fits"))
conseq$Year <- factor(as.character(conseq$Year), levels = c("2015", "2013", "2011"))

library(ggridges)
library(viridis)

plt2 <- ggplot(conseq, aes(x = accumdegday, y = Year, 
                           fill = Consequence, height = ..count..)) +
  # geom_histogram() +
  geom_density_ridges(scale = .95, rel_min_height=0.001, stat = "binline", bins = 100) +
  scale_fill_viridis(discrete = TRUE, alpha = .25, drop = FALSE) +
  # scale_x_continuous(limits = c(-50, maxgdd + 100), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0.1, .1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Year") + xlab("Accumulated degree-days (base 10C)") +
  facet_wrap(~SiteID, nrow = 3)
plt2







# TODO
# 1. Adult "counts" for each generation
# 2. Egg "counts" for each generation as alternative monitoring stages
# 3. Proportion diapause in the field from each generation ("typical" penultimate most important)
# 4. Adult counts to voltinism based on relative size to "typical" penultimate
# 5. Lineage analysis too much, what about comparing parents in penultimate generations and their CP/emergence time

census <- firsts %>% 
  group_by(SiteID, Year, cdl_mu, cdl_sd, lambda) %>% 
  pivot_longer(
    cols = one_of(names(firsts)[c(9:14, 16:26)]),
    names_to = c("choice", "numgen"),
    names_sep = "_",
    values_to = "count",
    values_drop_na = FALSE
  )

# voltinism by proportion diapause
# assumption: count the first generation where >.5 diapause
# problem: voltinism decimal never above x.5
census$count[is.na(census$count)] <- 0

diapvolt <- census %>% 
  # filter(cdl_mu %in% c(14, 15, 16)) %>% 
  group_by(SiteID, Year, cdl_mu, cdl_sd, lambda, numgen) %>% 
  summarise(propdiap = count[choice == "diap"] / count[choice == "att"],
            n = count[choice == "att"],
            voltdiap = as.numeric(numgen[1]) + (1 - propdiap),
            lost = lost[1]) %>% 
  ungroup() %>% 
  filter(complete.cases(.)) %>% 
  group_by(SiteID, Year, cdl_mu, cdl_sd, lambda) %>%  
  mutate(maxgen = max(numgen),
         keeprow = case_when(numgen == 1 ~ "yes",
                             propdiap > .95 ~ "no",
                             n < 200 ~ "no",
                             TRUE ~ "yes")) %>% 
  filter(keeprow != "no") %>% 
  arrange(-n) %>% 
  slice(1)
  


counts <- seconds %>% 
  group_by(SiteID, Year, cdl_mu, cdl_sd, lambda, numgen, lifestage) %>% 
  summarise(peakn = max(n), 
            peakdoy = doy[n == max(n)][1])

# how to remove weird years? what is typical voltinism
relvolt <- counts %>% 
  # filter(cdl_mu %in% c(14, 15, 16)) %>% 
  filter(lifestage == 5) %>%
  filter(lambda > 1) %>% 
  # filter(SiteID == "Vantage, WA", Year == 1994, cdl_mu <= 14, lifestage == 5) %>% 
  group_by(SiteID, Year, cdl_mu, cdl_sd, lambda) %>% 
  arrange(SiteID, Year, cdl_mu, cdl_sd, lambda, numgen) %>% 
  mutate(dpeak = c(0, diff(peakn))) %>% 
  summarise(maxgen = max(numgen),
            lastn = case_when(maxgen == 1 ~ peakn[1],
                              maxgen == 2 ~ peakn[2],
                              min(dpeak) < 0 ~ peakn[which(dpeak < 0)[1]],
                              TRUE ~ peakn[maxgen]),
            penultn = case_when(maxgen == 1 ~ peakn[1]-peakn[1],
                                maxgen == 2 ~ peakn[1],
                                maxgen > 2 & min(dpeak) < 0 ~ peakn[which(dpeak < 0)[1] - 1L],
                                # maxgen > 2 & min(dpeak) >= 0 ~ peakn[maxgen - 1],
                                TRUE ~ NA_integer_),
            penultgen = case_when(maxgen == 1 ~ as.integer(0L),
                                maxgen == 2 ~ as.integer(1L),
                                min(dpeak) < 0 ~ as.integer(numgen[which(dpeak < 0)[1]-1]),
                                TRUE ~ as.integer(maxgen - 1L)),
            lastprop = lastn / (lastn + penultn),
            voltcount = penultgen + lastprop)
                              

dat <- firsts %>% 
  mutate(voltmodel = attvoltinism) %>% 
  dplyr::select(SiteID, Year, cdl_mu, lambda, ann_lam, voltmodel) %>% 
  right_join(relvolt) %>% 
  dplyr::select(-maxgen) %>% 
  right_join(diapvolt) %>% 
  mutate(lost = log(lost + 1),
         ann_lam = log(ann_lam + 0.1))

pltdat <- dat %>% 
  group_by(SiteID, Year, cdl_mu, lambda) %>% 
  pivot_longer(
    cols = one_of(c("ann_lam", "voltmodel", "voltcount", "voltdiap", "lost")),
    names_to = "var",
    values_to = "value",
    values_drop_na = FALSE
  ) #%>% 
  # filter(lambda == 1.5)

pltdat$var <- as.factor(pltdat$var)
pltdat$var <- plyr::revalue(pltdat$var, c("ann_lam"="Log(Annual lambda)", 
                                                "voltmodel"="Modeled voltinism",
                                          "voltcount"="Counted voltinism",
                                          "voltdiap"="Diapause voltinism",
                                          "lost"="Log(Lost)"))
pltdat$var <- factor(pltdat$var, levels = c("Modeled voltinism", "Diapause voltinism", 
                                            "Counted voltinism", "Log(Lost)", "Log(Annual lambda)"))

ggplot(pltdat, aes(x = Year, y = value, group = cdl_mu, color = as.factor(cdl_mu))) +
  geom_line() +
  scale_color_brewer(type = "qual", palette = "Dark2", name = "Critical\nphotoperiod") +
  facet_grid(var~SiteID, scales = "free") +
  ylab("Annual value") +
  theme_bw()

cor.test(pltdat$value[which(pltdat$var == "Log(Annual lambda)")],
         pltdat$value[which(pltdat$var == "Modeled voltinism")], method = "kendall")
cor.test(pltdat$value[which(pltdat$var == "Log(Annual lambda)")],
         pltdat$value[which(pltdat$var == "Counted voltinism")], method = "kendall")
cor.test(pltdat$value[which(pltdat$var == "Log(Annual lambda)")],
         pltdat$value[which(pltdat$var == "Diapause voltinism")], method = "kendall")

cordat <- pltdat %>% 
  pivot_wider(names_from = var, values_from = value)

ggplot(cordat, aes(x = `Diapause voltinism`, y = `Log(Annual lambda)`, color = `Log(Lost)`)) +
  geom_point() +
  theme_bw()

ggplot(relvolt, aes(x = peakn, group = numgen, color = numgen)) +
  geom_density() +
  facet_grid(SiteID ~ cdl_mu, scales = "free")

ggplot(firsts %>% filter(lambda == 2, cdl_mu %in% c(14, 15, 16)), aes(x = Year, y = attvoltinism, group = cdl_mu, color = as.factor(cdl_mu))) +
  geom_line() +
  facet_wrap(~SiteID) +
  theme_bw()
ggplot(firsts %>% filter(lambda == 2, cdl_mu %in% c(14, 15, 16)), aes(x = Year, y = log(ann_lam), group = cdl_mu, color = as.factor(cdl_mu))) +
  geom_line() +
  facet_wrap(~SiteID) +
  theme_bw()
ggplot(firsts %>% filter(cdl_mu == 15, lambda %in% c(1, 2, 3)), aes(x = Year, y = attvoltinism, group = lambda, color = as.factor(lambda))) +
  geom_line() +
  facet_wrap(~SiteID) +
  theme_bw()
ggplot(firsts %>% filter(cdl_mu == 15, lambda %in% c(1, 2, 3)), aes(x = Year, y = log(ann_lam), group = lambda, color = as.factor(lambda))) +
  geom_line() +
  facet_wrap(~SiteID) +
  theme_bw()

library(ggpubr) # for combining plots

pltdat2 <- seconds %>% 
  # filter(cdl_mu == 13.5, lambda == 1.5, SiteID == "Sutherlin, OR", Year == 2015, lifestage == 5)
  filter(cdl_mu == 14, lambda == 1.5, SiteID == "Palermo, CA", Year == 2015, lifestage == 5)
  # filter(cdl_mu == 16, lambda == 1, SiteID == "Bellingham, WA", Year == 2014, lifestage == 5)
  # filter(cdl_mu == 16, lambda == 2, SiteID == "Bellingham, WA", Year == 2014, lifestage == 5)
  # filter(cdl_mu == 15.5, lambda == 1, SiteID == "Bellingham, WA", Year == 2014, lifestage == 5)
  # filter(cdl_mu == 16, lambda == 1, SiteID == "Bellingham, WA", Year == 2015, lifestage == 5)

pltdat2$diapause <- as.factor(pltdat2$diapause)
pltdat2$diapause <- plyr::revalue(pltdat2$diapause, c("-1"="Emergence", 
                                          "0"="Reproduction",
                                          "1"="Diapause"))
pltdat2$diapause <- factor(pltdat2$diapause, levels = c("Diapause", "Reproduction", "Emergence"))

b <- ggplot(pltdat2,
       aes(x = doy, y = n, fill = diapause)) +
  scale_fill_brewer(type = "qual", palette = "Dark2", name = "Adult decision") +
  # scale_x_continuous(limits = c(170, 350)) +
  # scale_y_continuous(limits = c(0, 1000)) +
  geom_histogram(stat = "identity", width = 1) +
  theme_bw(base_size = 16) +
  theme(legend.position = c(.8, .8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Count") +
  xlab("Day of year")
b

e <- ggarrange(a,b,c,d, 
                labels = c("A", "B", "C", "D"),
                ncol = 2, nrow = 2, align = "v", common.legend = TRUE)
e

e <- ggarrange(a,b, 
               labels = c("A", "B"),
               ncol = 2, nrow = 1, align = "h", common.legend = TRUE)
e

pltdat %>% 
  # filter(cdl_mu == 14, lambda == 1.5, SiteID == "Palermo, CA", Year == 2015) %>% 
  filter(cdl_mu == 13.5, lambda == 1.5, SiteID == "Sutherlin, OR", Year == 2015) %>%
# filter(cdl_mu == 16, lambda == 1, SiteID == "Bellingham, WA", Year == 2014) %>% 
# filter(cdl_mu == 16, lambda == 2, SiteID == "Bellingham, WA", Year == 2014) %>%
# filter(cdl_mu == 15.5, lambda == 1, SiteID == "Bellingham, WA", Year == 2014) %>%
# filter(cdl_mu == 16, lambda == 1, SiteID == "Bellingham, WA", Year == 2015) %>%
  data.frame()



# Quanitfying voltinism
library(viridis)
plt <- ggplot(res %>% filter(lambda == 1.5, Year <= 2018, cdl_sd == 0, SiteID == "Vantage, WA"), aes(x = cdl_mu, y = attvoltinism, color = log(ann_lam_adj_2500))) +
  geom_point() +
  scale_color_viridis(name = "Log(lambda)") +
  facet_wrap(~Year, ncol = 5) +
  theme_bw(base_size = 14) +
  xlab("Simulated critical photoperiod") +
  ylab("Mean number attempted") +
  theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plt





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
  res1 <- res %>% 
  filter(SiteID == "Vantage, WA", cdl_sd == .25, lambda == 1.5, Year > 1993)

plt <- ggplot(res1, aes(x = cdl_mu, y = attvoltinism)) +
  geom_point(aes(color = log(ann_lam + 0.01))) +
  scale_x_reverse() +
  scale_color_viridis(begin = 0, end = 1) +
  facet_wrap(~Year, ncol = 5) +
  theme_bw() +
  ylab("Voltinism: weighted average of generations completed") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated voltinism and proportion lost by critical photoperiod")
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
plt <- ggplot(res1, aes(x = cdl_mu, y = log(ann_lam_adj_2500), group = lambda, color = lambda)) +
  geom_line() +
  # scale_x_reverse() +
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

geom_lam <- firsts %>% 
  # filter(Year > 2007) %>% 
  group_by(SiteID, cdl_mu, cdl_sd, lambda) %>% 
  # tidyr::complete(Year, cdl_mu, cdl_sd, lambda, fill = list(ann_lam = 0)) %>% 
  summarise(mean_annual_lambda_0 = gm_mean(ann_lam + 0.001))
            # mean_annual_lambda_1 = gm_mean(ann_lam_adj_5000 + 0.001),
            # mean_annual_lambda_2 = gm_mean(ann_lam_adj_2500 + 0.001))

plt <- ggplot(geom_lam, aes(x = cdl_mu, y = lambda, fill = log(mean_annual_lambda_0))) +
  geom_raster() +
  # scale_x_reverse() +
  scale_fill_viridis(begin = 0, end = 1) +
  theme_bw() +
  facet_grid(~SiteID) +
  ylab("Between generation lambda") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated critical photoperiod and mean growth rates", subtitle = "1994-2018 geometric mean log(lambda)")
plt

# lines
plt <- ggplot(geom_lam, aes(x = cdl_mu, y = log(mean_annual_lambda_0), group = lambda, color = lambda)) +
  geom_path() +
  # scale_x_reverse() +
  scale_color_viridis(name = "Between generation\ngrowth rate", begin = .9, end = .1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(cdl_sd~SiteID) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  ylab("Log(geometric mean annual growth rate)") +
  xlab("Critical photoperiod")
  # ggtitle("Simulated log(Annual Population Growth Rate)", subtitle = "Accounting for lost generations")
plt

# Geometric mean vs annual
# lines
res$SiteID <- factor(res$SiteID, levels = levels(res$SiteID)[c(1,6,4,5,2,3)])


plt <- ggplot(res %>% filter(lambda == 2, cdl_sd == .25), aes(x = cdl_mu, y = log(ann_lam), group = Year, color = Year)) +
  geom_line(alpha = .5) +
  # scale_color_viridis(begin = .9, end = 0) +
  geom_line(data = geom_lam %>% filter(lambda == 2, cdl_sd == .25), aes(x = cdl_mu, y = log(mean_arith_lambda)), size = 1, inherit.aes = FALSE) +
  facet_wrap(~SiteID, ncol = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Log(annual lambda)") +
  xlab("Mean critical photoperiod")
  # ggtitle("Simulated log(Annual Population Growth Rate)", subtitle = "Accounting for lost generations")
plt




# best cdl combo
res <- readRDS("ibm_results_west.rds")
geom_lam <- res %>% 
  # left_join(maxvolt, by = c("Year" = "year")) %>% 
  left_join(maxvolt, by = c("Year" = "year", "SiteID" = "Site")) %>%
  group_by(SiteID, cdl_mu, cdl_sd, lambda) %>% 
  summarise(mean_annual_lambda = gm_mean(ann_lam + 0.001),
            mean_arith_lambda = mean(ann_lam + 0.001),
            mean_loglam = log(mean_annual_lambda),
            mean_lost = mean(lost),
            mean_wvolt = mean(ifelse(is.na(voltinism), 0, voltinism)),
            mean_attvolt = mean(attvoltinism),
            mean_mismatch = mean(attvoltinism - maxgen),
            sd_annual_lambda = sd(ann_lam + 0.001),
            sd_loglambda = sd(log(ann_lam + 0.001)),
            sd_lost = sd(lost),
            sd_wvolt = sd(ifelse(is.na(voltinism), 0, voltinism)),
            sd_attvolt = sd(attvoltinism),
            sd_mismatch = sd(attvoltinism - maxgen))
geom_lam$SiteID <- factor(geom_lam$SiteID, levels = levels(geom_lam$SiteID)[c(1,6,4,5,2,3)])

# ggplot(geom_lam, aes(x = mean_annual_lambda)) +
#   geom_density() +
#   facet_wrap(~lambda)

year_lam <- res %>% 
  left_join(maxvolt, by = c("Year" = "year", "SiteID" = "Site")) %>% 
  group_by(Year, SiteID, cdl_mu, cdl_sd, lambda) %>% 
  mutate(mean_mismatch = mean(attvoltinism - maxgen))
best_yr <- year_lam %>% 
  filter(lambda == 2, cdl_sd == .25) %>%
  group_by(SiteID, Year) %>%
  mutate(diffbest = ann_lam - max(ann_lam)) %>% 
  filter(diffbest == 0) %>% 
  arrange(SiteID, -cdl_mu) %>% 
  slice(1) %>% 
  mutate(type = "yr")
opt_site <- bestcdltab %>% 
  ungroup() %>% 
  dplyr::select(SiteID, cdl_mu) %>% 
  rename(bestcdl = cdl_mu)
opt_yr <- res %>% 
  left_join(opt_site) %>% 
  filter(cdl_mu == bestcdl) %>% 
  filter(lambda == 2, cdl_sd == .25) %>% 
  mutate(type = "avg",
         opt_volt = voltinism, 
         opt_ann_lam = ann_lam) %>% 
  dplyr::select(SiteID, Year, bestcdl, opt_volt, opt_ann_lam)

# mismatch <- bind_rows(best_yr, opt_yr) %>% 
#   group_by(SiteID, Year) %>% 
#   arrange(SiteID, Year, type) %>% 
#   summarise(cdl_diff = diff(cdl_mu),
#             lam_diff = diff(ann_lam),
#             volt_diff = diff(voltinism))

meanmm <- mismatch %>% 
  group_by(SiteID) %>% 
  summarise_at(vars(cdl_diff:volt_diff), mean)

mismatch <- left_join(best_yr, opt_yr) %>% 
  mutate(voltmm = opt_volt - voltinism,
         cdlmm = bestcdl - cdl_mu )
mismatch$SiteID <- factor(mismatch$SiteID, levels = levels(res$SiteID))
best_yr$SiteID <- factor(best_yr$SiteID, levels = levels(res$SiteID))

# comparing yearly best choice versus geom mean best choice
plt <- ggplot(mismatch, aes(x = log(ann_lam), y = log(opt_ann_lam), color = cdlmm)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~SiteID, ncol = 2)
plt


plt <- ggplot(best_yr, aes(x = cdl_mu)) +
  geom_histogram(bins = 40, fill = "grey") +
  facet_wrap(~SiteID, ncol = 2) +
  geom_vline(data = opt_site, aes(xintercept = bestcdl), linetype = "dashed") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Mean critical photoperiod")
  
plt




# maximize geometric mean of annual lambda (non including overwintering)
best_cdl <- geom_lam %>% 
  # filter(lambda == 1.5) %>%
  group_by(SiteID, lambda, cdl_sd) %>%
  # mutate(diffbest = mean_annual_lambda - max(mean_annual_lambda)) %>% 
  mutate(diffbest = mean_arith_lambda - max(mean_arith_lambda)) %>% 
  filter(diffbest == 0) %>% 
  arrange(cdl_mu) %>% 
  slice(1)

best_cdl$SiteID <- factor(best_cdl$SiteID, levels = levels(best_cdl$SiteID)[c(1,6,4,5,2,3)])


plt <- ggplot(best_cdl, aes(x = lambda, y = cdl_mu, group = as.factor(cdl_sd), color = as.factor(cdl_sd))) +
  geom_point(  position = position_dodge(width = .2)) +
  scale_color_viridis(name = "Std. dev. CP", discrete = TRUE, begin = .2, end = .8) +
  facet_wrap(~SiteID, ncol = 2) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Optimal mean critical photoperiod") +
  xlab("Between-generation population growth rate")
plt

bestcdltab <- best_cdl %>% 
  filter(cdl_sd == .25, lambda == 2) %>% 
  arrange(SiteID)
write.csv(round_df(bestcdltab), "bestcdl.csv", row.names = FALSE)


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


# extract voltinism/lost/attempted info for the row with best cdl
for (site in unique(res$SiteID)){
  best <- best_cdl %>% filter(SiteID == site)
  tmp <- res %>% 
    filter(SiteID == site, cdl_mu == best$cdl_mu, cdl_sd == best$cdl_sd, lambda == best$lambda)
  
  attempts <- tmp %>% 
    tidyr::gather(gen, num, starts_with("att"))
  
  
}


# site variables from gdd data
site_gdd <- gdd_west %>% 
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
# summary(lm(cdl_mu ~ scale(mean_accumdd), data = dat))
# summary(lm(cdl_mu ~ scale(lat), data = dat))

region <- rgdal::readOGR("./src/ref/ne_50m_admin_1_states_provinces_lakes", 'ne_50m_admin_1_states_provinces_lakes', encoding='UTF-8')
# region <-spTransform(region, CRS(proj4string(template)))
reg.points = fortify(region, region="name_en")
reg.df = left_join(reg.points, region@data, by = c("id" = "name_en"))
theme_set(theme_bw(base_size = 20))


tmpplt <- ggplot() +
  geom_polygon(data = reg.df, aes(x = long, y = lat, group = group), fill = NA, color = "light gray", inherit.aes = FALSE, size = 1, alpha = .3) +
  # geom_point(data = dat, aes(x = lon, y = lat, color = as.factor(round(mean_attvolt))), alpha = 1, size = 3, inherit.aes = FALSE) +
  # scale_color_viridis(discrete = TRUE, option = "C") +
  geom_point(data = dat, aes(x = lon, y = lat, color = cdl_mu), alpha = 1, size = 4, inherit.aes = FALSE) +
  scale_color_viridis(discrete = FALSE) +
  coord_fixed(1.3, xlim = c(-126.5, -119.5), ylim = c(36, 52), expand = FALSE, clip = "on") +
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





