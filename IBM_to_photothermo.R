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

# Calculate degree days and photoperiod ----

gdd_all <- gdd %>%
  filter(yday <= end_doy & yday >= start_doy) %>%
  rowwise() %>%
  mutate(degday = TriDD(tmax..deg.c., tmin..deg.c., ldt, udt)) %>%
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
maxvolt <- gdd_west %>% 
  filter(yday == 365) %>% 
  summarise(maxgen = mean(floor((accumdegday - ow_dd) / gen_dd)))

maxvolt <- gdd_all %>% 
  filter_all(all_vars(!is.na(.)))  %>% 
  group_by(Site, year) %>% 
  arrange(yday) %>%
  mutate(nday = n(),
         accumdegday = cumsum(degday)) %>% 
  dplyr::filter(yday == max(yday),
                nday > 350) %>% 
  summarise(maxgen = mean(floor((accumdegday - ow_dd) / gen_dd)))


# 3. Run lifestage model -----


sites <- unique(sites$ID)
cdls <- expand.grid(cdl = c(14.2, 14.7, 15.2, 15.7, 16.2),
                    cdl_sd = 0.25,
                    lambda = 2,
                    site = sites[3],
                    sim = 1)


gdd_west <- gdd_all
startyear <- min(gdd_west$year)
endyear <- max(gdd_west$year)

ncores <- 8
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
              ind_df$cdl <- NA
              ind_df$diap_gdd <- 0
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
                newdf$cdl <- Cond(diap != newdf$diapause, 0, CDL)
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
                filter(diapause == 1) %>% 
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
              
              ann_lam_adj <- sum(ind_df %>% 
                                   filter(diapause == 1 & active == 0) %>% 
                                   mutate(gdd_penalty = 1 - diap_gdd / 5000) %>% 
                                   pull(gdd_penalty))/ nsim
              
              
              
              
              results <- data.frame(SiteID = cdls[ncdl, "site"], Year = thisyr, cdl_mu, cdl_sd, 
                                    lambda, voltinism, attvoltinism, lost, volt_att, volt_comp, 
                                    diap, repro, ann_lam, ann_lam_adj, sim = cdls[ncdl, "sim"])
              
              ind_df$Year <- thisyr
              ind_df$sim <- cdls[ncdl, "sim"]
              ind_df$cdl_mu <- cdl_mu
              ind_df$cdl_sd <- cdl_sd
              ind_df$lambda <- lambda
              
              res_list <- list()
              res_list[[1]] <- ind_df
              res_list[[2]] <- results
              return(res_list)
              
            } # close foreach loop
}) # system.time

stopCluster(cl)


res <- flatten(flatten(outlist))


years <- c(1993:2018)
years <- c(2011:2015)


inds <- bind_rows(res[seq(1, length(res), by = 2)]) %>% 
  filter(cdl_mu == 14.2, Year %in% years)
volt <- bind_rows(res[seq(2, length(res), by = 2)])
gdd_yr <- gdd_all %>% filter(year %in% years)


outvolt <- volt %>% 
  filter(Year %in% years)

comp <- outvolt %>% 
  group_by(Year, cdl_mu) %>% 
  dplyr::select(Year, cdl_mu, contains("comp_")) %>% 
  tidyr::gather(gen, num, contains("comp_")) %>% 
  mutate(gen = as.numeric(str_split_fixed(gen, pattern = coll("_"), 2)[,2]),
         param = "completed")

att <- outvolt %>% 
  group_by(Year, cdl_mu) %>% 
  dplyr::select(Year, cdl_mu, contains("att_")) %>% 
  tidyr::gather(gen, num, contains("att_")) %>% 
  mutate(gen = as.numeric(str_split_fixed(gen, pattern = coll("_"), 2)[,2]),
         param = "attempted")

diap <- outvolt %>% 
  group_by(Year, cdl_mu) %>% 
  dplyr::select(Year, cdl_mu, contains("diap_")) %>% 
  tidyr::gather(gen, num, contains("diap_")) %>% 
  mutate(gen = as.numeric(str_split_fixed(gen, pattern = coll("_"), 2)[,2]),
         param = "diapause")

repro <- outvolt %>% 
  group_by(Year, cdl_mu) %>% 
  dplyr::select(Year, cdl_mu, contains("repro_")) %>% 
  tidyr::gather(gen, num, contains("repro_")) %>% 
  mutate(gen = as.numeric(str_split_fixed(gen, pattern = coll("_"), 2)[,2]),
         param = "reproduce")

genvolt <- bind_rows(comp, att, diap, repro)

outvolt2 <- outvolt %>% 
  dplyr::select(SiteID, Year, cdl_mu, voltinism, attvoltinism, lost, ann_lam, ann_lam_adj) %>% 
  left_join(genvolt) %>% 
  tidyr::spread(param, num, fill = 0) %>% 
  mutate(perc_diap = completed / attempted, 
         perc_repro = reproduce / attempted,
         perc_lost = (attempted - completed - reproduce) / (attempted),
         contribution = completed / 500) %>% 
  group_by(SiteID, Year, cdl_mu) %>% 
  arrange(SiteID, cdl_mu, Year, gen) %>%
  ungroup() %>% 
  filter(attempted > 0) %>% 
  dplyr::select(cdl_mu, Year, gen, attempted, perc_diap, perc_repro, perc_lost, contribution, attvoltinism)
  # mutate(obs_lambda = c(attempted[2:(length(attempted))]/attempted[1:(length(attempted)-1)], NA))

write.csv(outvolt2 %>% filter(cdl_mu == 14.2), "table14.2.csv")

outlist <- list()
for (yr in years){
  
  temp <- inds %>% filter(Year == yr)
  gdd <- gdd_yr %>% filter(year == yr) %>% mutate(doy = yday)
  
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
  outlist[[length(outlist) + 1]] <- outcomes
}

outcomes <- bind_rows(outlist)

# Marginal histograms ----
mean_cdl <- mean(inds$cdl, na.rm = TRUE)
maxgdd <- max(gdd_yr$accumdegday, na.rm = TRUE)
maxphoto <- max(c(max(gdd_yr$daylength), max(inds$cdl, na.rm = TRUE)))
minphoto <- min(gdd_yr$daylength)

cdl_hist <- inds %>% 
  ungroup() %>% 
  dplyr::select(cdl) %>%
  filter(complete.cases(.)) %>% 
  mutate(breaks = cut(cdl, breaks=seq(min(cdl),max(cdl),length.out = 20), 
                      labels=seq(min(cdl),max(cdl), length.out = 20)[-1], 
                      include.lowest=TRUE)) %>% 
  mutate(daylength = as.numeric(as.character(breaks))) %>%
  group_by(daylength) %>% 
  summarise(n = n()) %>% 
  mutate(accumdegday = -50 + 50 * n/max(n))

oa <- inds %>% filter(OA_DOY > 0) %>% 
  left_join(gdd_all, by = c("Year" = "year", "OA_DOY" = "yday"))
teneral <- inds %>% filter(TA_DOY > 0) %>% 
  left_join(gdd_all, by = c("Year" = "year", "TA_DOY" = "yday"))
sens_stage <- bind_rows(oa, teneral)
sens_hist <- sens_stage %>%
  ungroup() %>% 
  dplyr::select(accumdegday) %>%
  filter(complete.cases(.)) %>% 
  mutate(breaks = cut(accumdegday, breaks = seq(min(accumdegday), max(accumdegday), length.out = 100),
                      labels = seq(min(accumdegday), max(accumdegday), length.out = 100)[-1],
                      include.lowest = TRUE)) %>% 
  mutate(accumdegday = as.numeric(as.character(breaks))) %>% 
  group_by(accumdegday) %>% 
  summarise(n = n()) %>% 
  mutate(daylength = minphoto - .5 + .5 * n/max(n))

# Photothermograph ----
theme_set(theme_bw(base_size = 12)) 

plt <- ggplot(gdd_yr, aes(x = accumdegday, y = daylength, group = as.factor(year), color = as.factor(year))) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(-50, maxgdd + 100), ylim = c(minphoto - .5, 17.2), expand = FALSE) +
  scale_color_viridis(discrete = TRUE, name = "Year", option = "D") +
  scale_y_continuous(breaks = c(9:17)) +
  guides(color = guide_legend(reverse=FALSE)) +
  geom_segment(data=cdl_hist, size=2.5, show.legend=FALSE,
               aes(x=-50, xend=accumdegday, y=daylength, yend=daylength), inherit.aes = FALSE) +
  geom_segment(data=sens_hist, size=2.5, show.legend=FALSE,
               aes(x=accumdegday, xend=accumdegday, y=minphoto - .5, yend=daylength), inherit.aes = FALSE) +
  geom_hline(yintercept = mean_cdl, linetype = "dashed", size = 1, alpha = .5) +
  # geom_rug(data = results, aes(x = accumdegday), color = "gray", alpha = .1) +
  # geom_rug(data = results, aes(y = cdl), color = "gray", alpha = .1) +
  # ggtitle("Photothermographs with critical photoperiod\nand simulated sensitive stage emergence") +
  xlab("Accumulated degree-days") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(.25, .3)) +
  ylab("Daylength (hours)")
plt



# Voltinism cost ----
conseq <- outcomes %>%
  mutate(Consequence = case_when(outcome == "diapause" ~ "Diapauses",
                                 offspring == "lost" ~ "Reproduces, offspring lost",
                                 offspring != "lost" ~ "Reproduces, offspring fits"),
         date = as.Date(doy, origin=as.Date("2015-12-31")))

conseq$Consequence <- factor(conseq$Consequence, levels = c("Diapauses", "Reproduces, offspring lost", "Reproduces, offspring fits"))


plt2 <- ggplot(conseq, aes(x = accumdegday, y = as.factor(Year), 
                           fill = Consequence, height = ..count..)) +
  # geom_density_ridges(scale = .9, rel_min_height=0.01, stat = "density") +
  geom_density_ridges(scale = .95, rel_min_height=0.001, stat = "binline", bins = 100) +
  scale_fill_viridis(discrete = TRUE, alpha = .25, drop = FALSE) +
  scale_x_continuous(limits = c(-50, maxgdd + 100), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0.1, .1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(.85, .2)) +
  ylab("Year") + xlab("Accumulated degree-days (base 10C)")
# facet_wrap(~year, nrow = 2)
plt2


plt3 <- ggplot(conseq, aes(x = date, y = as.factor(Year), 
                           fill = Consequence, height = ..count..)) +
  # geom_density_ridges(scale = .9, rel_min_height=0.01, stat = "density") +
  geom_density_ridges(scale = .95, rel_min_height=0.001, stat = "binline", bins = 100) +
  scale_fill_viridis(discrete = TRUE, alpha = .25, drop = FALSE) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = c(as.Date("2016-04-01"), as.Date("2016-12-01"))) +
  scale_y_discrete(expand = c(0.1, .1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  ylab("Year") + xlab("Month")
plt3

# library(ggpubr)
# p4 <- ggarrange(plt2, plt3, common.legend = TRUE, legend = "bottom",
                # ncol = 1, nrow = 2, align = "v")
# p4

p5 <- ggarrange(plt + rremove("xlab"), plt2, plt3,
                ncol = 1, nrow = 3, align = "v")
p5


ggsave(plot = p5, filename =  "gca_yak_sim_14.2.png", device = "png", width = 8, height = 12, units = "in")







theme_set(theme_bw(base_size = 12)) 

gdd_all$Site <- factor(gdd_all$Site, levels = levels(gdd_all$Site)[c(1,6,4,5,2,3)])


plt <- ggplot(gdd_all, aes(x = accumdegday, y = daylength, group = as.factor(year), color = as.factor(year))) +
  geom_line(size = 1, color = "black", alpha = .2) +
  # coord_cartesian(xlim = c(-50, maxgdd + 100), ylim = c(minphoto - .5, 17.2), expand = FALSE) +
  # scale_color_viridis(discrete = TRUE, name = "Year", option = "D") +
  scale_y_continuous(breaks = c(9:17)) +
  guides(color = guide_legend(reverse=FALSE)) +
  xlab("Accumulated degree-days") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  ylab("Daylength (hours)") +
  facet_wrap(.~Site, ncol = 2)
plt

