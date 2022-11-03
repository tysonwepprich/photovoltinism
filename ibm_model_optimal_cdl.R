# Speeding up raster-based model with arrays (cohort x pixel)
# Universal code for 3 biocontrol species
# Any number of lifestages allowed
# IBM version


# 1. Setup -------
# packages, options, functions loaded
# TODO: find ways to take out dplyr, purrr, mixsmsn functions?
pkgs <- c("sp", "rgdal", "rgeos", "raster", "lubridate", "mixsmsn", "daymetr",
          "stringr", "prism", "purrr", "foreach", "doParallel", "abind", "ncdf4",
          "ggplot2", "dplyr", "tidyr")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them

# load collection of functions for this model
source('CDL_funcs.R')
source('species_params_ibm.R')

# 2. User input -----

# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model scope
# yr           <- 2018
start_doy    <- 1
end_doy      <- 365
# region_param <- "WEST" # TEST/WEST/EAST/CONUS/SOUTHWEST/NORTHWEST
species      <- "GCA" # GCA/APHA/DCA
biotype      <- "Delta" # TODO: add options for each species, N or S for APHA and GCA

nsim <- 100

# lat <- 43.39
# lon <- -123.35
startyear <- 1992
endyear <- 2018
# site <- "Sutherlin"


# Load sites DCA
sites <- read.csv("data/DCA_sitecoords.csv", header = TRUE)[,c(2, 4, 5)]
names(sites) <- c("ID", "y", "x") #### CHECK THIS
sites <- sites[-which(sites$ID %in% c("Pueblo, CO", "Artesia, NM")), ]
# Load sites GCA
sites <- read.csv("data/GCA_modeling_sites.csv", header = TRUE) [c(4, 21, 22),]
# names(sites) <- c("ID", "y", "x") #### CHECK THIS



# REGION <- assign_extent(region_param = region_param)

# photoperiod decision inclusion
# 2 for logistic, 1 for single value CDL, 0 for none
model_CDL  <- 1


# Derived parameters -----
params <- species_params(mod_type = "ibm", species, biotype, nsim, model_CDL, dd_sd = 5)

ldt <- params$stage_ldt[1]
udt <- params$stage_udt[1]

# Download weather data ----

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


# Calculate degree days and photoperiod ----

# g1 <- readRDS("ngermany.rds")
# g2 <- readRDS("sgermany.rds")
# gdd_all <- bind_rows(g1,g2) %>% 
#   rename(Year = YEAR, yday = YDAY)
# 
# gdd <- readRDS("lythrum_gdd.rds")

gdd_all <- gdd %>%
  dplyr::filter(yday <= end_doy & yday >= start_doy) %>%
  rowwise() %>%
  mutate(degdayOLD = TriDD(tmax..deg.c., tmin..deg.c., ldt, udt),
         degdayVERT = TriDDvert(tmax..deg.c., tmin..deg.c., ldt, udt),
         degdayTRI = degree_days(tmin..deg.c., tmax..deg.c., ldt, udt, method = "single.triangulation"),
         degday = degree_days(tmin..deg.c., tmax..deg.c., ldt, udt, method = "single.sine")) %>%
  ungroup() %>%
  group_by(Site, year) %>%
  arrange(yday) %>%
  mutate(accumdegdayOLD = cumsum(degdayOLD),
         accumdegdayVERT = cumsum(degdayVERT),
         accumdegdayTRI = cumsum(degdayTRI),
         accumdegday = cumsum(degday),
         daylength = photoperiod(Latitude, yday),
         daylength_daymet = dayl..s. / 60 / 60,
         lightfrost = tmin..deg.c. <= 0,
         hardfrost = tmin..deg.c. <= -2)

# try out vertical cutoff
gdd_all <- gdd_all %>% 
  mutate(degday = degdayVERT,
         accumdegday = accumdegdayVERT)


# Summarize sites' seasons----
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
  
write.csv(gdd_seas, "gdd_seas_DCA.csv")  

tapply(gdd_all$daylength, gdd_all$Site, FUN = max)  

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


# naive max voltinism by gdd alone
ow_dd <- params$stage_dd[,1]
gen_dd <- rowSums(params$stage_dd[,-1])

# calculated in different ways, average away individual variation or not?
maxvolt <- gdd_all %>% 
  dplyr::filter_all(all_vars(!is.na(.)))  %>% 
  group_by(Site, year) %>% 
  arrange(yday) %>%
  mutate(nday = n(),
    accumdegday = cumsum(degday)) %>% 
  dplyr::filter(yday == max(yday),
         nday > 350) %>% 
  summarise(maxgen = floor((accumdegday - mean(ow_dd)) / mean(gen_dd)))

# Non-integer output of maxgen accounting for ind. variation
maxvolt <- gdd_all %>%
  dplyr::filter(yday == 365) %>%
  summarise(maxgen = mean(floor((accumdegday - ow_dd) / gen_dd)),
            accumdegday = accumdegday)

# Plot: Voltinism by degree-days (constant) ----
# example years 2011,13,15 for cool to warm
# G&C2015 uses 100dd OW ovip, 497 1st gen, 523 subsequent adult emergence

expand_df <- function(df){
  newdf <- data.frame(gen = 0:df$maxgen) %>% 
    mutate(ddreq = median(ow_dd) + gen * median(gen_dd)) %>% 
    rowwise() %>% 
    mutate(accumdegday = gdd_all %>% 
             dplyr::filter(Site == df$Site[1] & year == df$year[1]) %>% 
             dplyr::filter(accumdegday - ddreq[1] >= 0) %>% 
             arrange(yday) %>% 
             slice(1) %>% pull(accumdegday),
           yday = gdd_all %>% 
             dplyr::filter(Site == df$Site[1] & year == df$year[1]) %>% 
             dplyr::filter(accumdegday - ddreq[1] >= 0) %>% 
             arrange(yday) %>% 
             slice(1) %>% pull(yday))
  return(newdf)
}

milestones <- maxvolt %>% 
  dplyr::filter(year %in% c(2011, 2013, 2015)) %>% 
  group_by(Site, year) %>% 
  do(expand_df(.))
milestones$year <- as.character(milestones$year)
milestones$date <- as.Date(milestones$yday, origin=as.Date("2015-12-31"))




dds <- gdd_all %>% dplyr::filter(year %in% c(2011, 2013, 2015)) %>% 
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


# Plot: Voltinism by photothermograph (constant)----
# add photosensitive time to milestones
cdl <- 15
expand_df <- function(df){
  cycle <- "reproduce"
  newdf <- data.frame(gen = 1, start = median(ow_dd), end = median(ow_dd) + 1 * median(gen_dd),
                      sense = median(ow_dd) + 1 * median(gen_dd) - median(rowSums(params$stage_dd[,c(5:6)]))) %>% 
    mutate(daylength = df %>% dplyr::filter(accumdegday - sense >= 0) %>% 
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
                                  df %>% dplyr::filter(accumdegday - sense >= 0) %>% 
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


# 3. Run lifestage model -----

testsites <- unique(gdd_all$Site)
cdls <- expand.grid(cdl = seq(10, 20, by = .25),
                    cdl_sd = seq(0, 1.5, by = .25),
                    lambda = c(2),
                    site = testsites)
# cdls <- expand.grid(cdl = seq(12, 18, by = 1),
                    # cdl_sd = 1,
                    # lambda = 2,
                    # site = testsites)


# Run choice sites/years
inc_sites <- c("Bellingham, WA", "Sutherlin, OR", "Palermo, CA")
cdls <- cdls %>% 
  filter(site %in% inc_sites, lambda == 1.5, cdl_sd == .5, cdl == 15)

startyear <- min(gdd_all$year)
endyear <- max(gdd_all$year)
startyear <- 2011
endyear <- 2015

ncores <- 20
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
              
              # print(ncdl)
              # print(y)
              
              allyrs <- c(startyear:endyear)
              thisyr <- allyrs[y]
              gdd <- gdd_all %>% 
                dplyr::filter(year == thisyr & Site == cdls[ncdl, "site"])
              # dplyr::filter(year == thisyr & SiteID == sites$ID[site])
              
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
                  # deal with lambda 1.5/2.5 by flipping between integers
                  negg <- round(egg_lam + base::rep(c(0.1, -0.1), length.out = length(egg_lam)))
                  newrows <- newrows[base::rep(row.names(newrows), negg), ]
                  
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
                
                # # daily population counts to save
                # pop_list[[index]] <- ind_df %>% 
                #   group_by(numgen, lifestage, active, diapause) %>% 
                #   tally() %>% 
                #   mutate(doy = index)
              }  # daily loop
              
              # Summarize results -----
              voltinism <- ind_df %>% 
                dplyr::filter(diapause == 1 & active == 0) %>% 
                group_by(numgen) %>% 
                tally() %>% 
                stats::weighted.mean(x = .$numgen, w = .$n)
              
              if(is.nan(voltinism)){ 
                voltinism <- data.frame(voltinism = NA)
              }
              
              attvoltinism <- ind_df %>% 
                dplyr::filter(diapause == 1 | active == 1) %>% 
                # dplyr::filter(numgen > 0) %>% 
                group_by(numgen) %>% 
                tally() %>% 
                stats::weighted.mean(x = .$numgen, w = .$n)

              # redundant with diap below              
              # volt_comp <- ind_df %>%
              #   dplyr::filter(diapause == 1 & active == 0) %>% 
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
                dplyr::filter(numgen > 0) %>% 
                group_by(numgen) %>% 
                tally() %>% 
                tidyr::spread(numgen, n)
              
              names(volt_att) <- paste("att", names(volt_att), sep = "_")
              
              
              lost <- ind_df %>% 
                dplyr::filter(active == 1) %>% 
                group_by(numgen) %>% 
                nrow()
              
              maxgdd <- gdd$accumdegday[which(gdd$yday == 365)]
              diap_gdd <- ind_df %>%
                dplyr::filter(diapause == 1 & active == 0) %>%
                left_join(gdd[,c("yday", "accumdegday")], by = c("TA_DOY" = "yday")) %>% 
                group_by(numgen) %>%
                summarise(meandiapgdd = mean(maxgdd - accumdegday)) %>%
                tidyr::spread(numgen, meandiapgdd)
              
              if(nrow(diap_gdd) == 0){ 
                diap_gdd <- data.frame(diapgdd_1 = NA)
              }else{
                names(diap_gdd) <- paste("diapgdd", names(diap_gdd), sep = "_")
              }
              
              diap_doy <- ind_df %>% 
                dplyr::filter(diapause == 1 & active == 0) %>%
                group_by(numgen) %>% 
                summarise(meandiapdoy = mean(TA_DOY)) %>% 
                tidyr::spread(numgen, meandiapdoy)
              
              if(nrow(diap_doy) == 0){ 
                diap_doy <- data.frame(diapdoy_1 = NA)
              }else{
                names(diap_doy) <- paste("diapdoy", names(diap_doy), sep = "_")
              }

              diap <- ind_df %>%
                dplyr::filter(diapause == 1 & active == 0) %>% 
                group_by(numgen) %>% 
                tally() %>% 
                tidyr::spread(numgen, n)
              
              if(nrow(diap) == 0){
                diap <- data.frame(diap_1 = NA)
              }else{
                names(diap) <- paste("diap", names(diap), sep = "_")
              }
              
              
              repro <- ind_df %>%
                dplyr::filter(diapause == 0) %>% 
                group_by(numgen) %>% 
                tally() %>% 
                tidyr::spread(numgen, n)
              
              if(nrow(repro) == 0){
                repro <- data.frame(repro_1 = NA)
              }else{
                names(repro) <- paste("repro", names(repro), sep = "_")
              }
              
              
              ann_lam <- (ind_df %>% 
                            dplyr::filter(diapause == 1 & active == 0) %>% 
                            nrow()) / nsim
              
              ann_lam_adj_5 <- sum(ind_df %>%
                            dplyr::filter(diapause == 1 & active == 0) %>%
                            left_join(gdd[,c("yday", "accumdegday")], by = c("TA_DOY" = "yday")) %>% 
                            mutate(diap_gdd = maxgdd - accumdegday) %>% 
                            mutate(gdd_penalty = 1 - diap_gdd * .05 / 1000) %>%
                            pull(gdd_penalty))/ nsim
              
              ann_lam_adj_10 <- sum(ind_df %>%
                                     dplyr::filter(diapause == 1 & active == 0) %>%
                                     left_join(gdd[,c("yday", "accumdegday")], by = c("TA_DOY" = "yday")) %>% 
                                     mutate(diap_gdd = maxgdd - accumdegday) %>% 
                                     mutate(gdd_penalty = 1 - diap_gdd * .10 / 1000) %>%
                                     pull(gdd_penalty))/ nsim
              
              
              results <- data.frame(SiteID = cdls[ncdl, "site"], Year = thisyr, cdl_mu, cdl_sd, 
                                    lambda, voltinism, attvoltinism, lost, volt_att,  
                                    diap, diap_gdd, diap_doy, repro, ann_lam, 
                                    ann_lam_adj_5, ann_lam_adj_10)
              
              # thinking about which lineages are doing better or worse
              
              # oa <- ind_df %>%
              #   dplyr::filter(numgen == 0) %>% 
              #   dplyr::select(lineage, OA)
              # cp <- ind_df %>% 
              #   dplyr::filter(numgen > 0) %>% 
              #   group_by(lineage, numgen) %>% 
              #   summarise(meancp = mean(ind_df$cdl[ind_df$indiv %in% parent]))
              # lam <- ind_df %>% 
              #   group_by(lineage) %>% 
              #   summarise(lam = length(which(diapause == 1 & active == 0)))
              # linlost <- ind_df %>% 
              #   group_by(lineage, numgen) %>% 
              #   summarise(linlost = length(which(active == 1))) %>% 
              #   dplyr::filter(linlost > 0)
              # if(length(unique(linlost$numgen)) == 1){
              #   cp <- cp %>% 
              #     dplyr::filter(numgen == unique(linlost$numgen)[1] - 1)
              #   linlost <- linlost %>% 
              #     ungroup() %>% 
              #     dplyr::select(-numgen)
              # }
              # lineages <- right_join(oa, cp) %>% 
              #   left_join(lam) %>% 
              #   left_join(linlost)
              # lineages$linlost[is.na(lineages$linlost)] <- 0
                              
              out <- list()
              out$indiv <- ind_df %>%
                mutate(SiteID = cdls[ncdl, "site"],
                       Year = thisyr,
                       cdl_mu, cdl_sd,  lambda)
              out$results <- results
              # out$pop <- bind_rows(pop_list) %>% 
              #   mutate(SiteID = cdls[ncdl, "site"], 
              #          Year = thisyr, 
              #          cdl_mu, cdl_sd,  lambda)
              return(out)
              

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

saveRDS(outlist, "ibm_outlist_DCA_12sites_daymet.rds")

# quirks: indiv can decide to diapause at beginning of TA, but still be active because is doesn't make it to complete TA stage before end of year
# res <- bind_rows(outlist)
# res <- bind_rows(flatten(outlist)) 
# res <- bind_rows(flatten(outlist)) %>% # if >1 sim per combo
  # group_by(Year, SiteID, cdl_mu, cdl_sd, lambda) %>% 
  # summarise_at(vars(voltinism:diap_2), mean, na.rm = TRUE)
# saveRDS(res, "ibm_results_GCA_6sites_daymet.rds")

# res <- readRDS("ibm_outlist_GCA_6sites_daymet.rds")
# inds <- outlist %>% 
  # flatten() %>%
  # purrr::map(1) %>% dplyr::bind_rows()

firsts <- outlist %>% 
  flatten() %>%
  purrr::map(1) %>% dplyr::bind_rows()
seconds <- outlist %>% 
  flatten() %>%
  purrr::map(2) %>% dplyr::bind_rows()

# firsts$SiteID <- plyr::revalue(firsts$SiteID, c("Yakima Training Center, WA"="Vantage, WA", 
#                                                       "Baskett Slough, OR"="Rickreall, OR"))
# firsts$SiteID <- factor(firsts$SiteID, levels = c("Bellingham, WA", "Vantage, WA", "Rickreall, OR",
#                                                         "Sutherlin, OR", "McArthur, CA", "Palermo, CA"))
# seconds$SiteID <- plyr::revalue(seconds$SiteID, c("Yakima Training Center, WA"="Vantage, WA", 
#                                                 "Baskett Slough, OR"="Rickreall, OR"))
# seconds$SiteID <- factor(seconds$SiteID, levels = c("Bellingham, WA", "Vantage, WA", "Rickreall, OR",
#                                                   "Sutherlin, OR", "McArthur, CA", "Palermo, CA"))





