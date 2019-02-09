# Speeding up raster-based model with arrays (cohort x pixel)
# Universal code for 3 biocontrol species
# Any number of lifestages allowed
# IBM version


# 1. Setup -------
# packages, options, functions loaded
# TODO: find ways to take out dplyr, purrr, mixsmsn functions?
pkgs <- c("sp", "rgdal", "rgeos", "raster", "lubridate", "mixsmsn", "dplyr", "daymetr",
          "stringr", "prism", "purrr", "foreach", "doParallel", "abind", "ncdf4")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them

# load collection of functions for this model
source('CDL_funcs.R')
source('species_params.R')

# 2. User input -----

# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model scope
yr           <- 2015
start_doy    <- 1
end_doy      <- 365
# region_param <- "NORTHWEST" # TEST/WEST/EAST/CONUS/SOUTHWEST/NORTHWEST
species      <- "GCA" # GCA/APHA/DCA
biotype      <- "S" # TODO: add options for each species, N or S for APHA and GCA

nsim <- 100
# lambda <- 1.1
lat <- 43.39
lon <- -123.35
startyear <- 2014
endyear <- 2017
site <- "Sutherlin"

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
temp <- download_daymet(site = "default", lat = lat, lon = lon,
                        start = startyear, end = endyear, internal = TRUE,
                        silent = TRUE, force = FALSE)

# Calculate degree days and photoperiod ----
gdd_all <- temp$data %>% 
  filter(yday <= end_doy & yday >= start_doy) %>% 
  rowwise() %>% 
  mutate(degday = TriDD(tmax..deg.c., tmin..deg.c., ldt, udt)) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  arrange(yday) %>% 
  mutate(accumdegday = cumsum(degday),
         daylength = photoperiod(lat, yday),
         lightfrost = tmin..deg.c. <= 0,
         hardfrost = tmin..deg.c. <= -2)




# 3. Run lifestage model -----

cdls <- expand.grid(cdl = seq(min(gdd_all$daylength), 18, length.out = 20),
                    cdl_sd = seq(0, 2, length.out = 4),
                    lambda = seq(0.6, 3, by = 0.4))


ncores <- 20
cl <- makePSOCKcluster(ncores)
registerDoParallel(cl)

test <- system.time({
  
  outlist <- foreach(ncdl = 1:nrow(cdls),
                     .packages = "dplyr",
                     .inorder = FALSE) %:%
    foreach(y = 1:length(startyear:endyear),
            .packages = "dplyr",
            .inorder = FALSE) %dopar%{
              
              allyrs <- c(startyear:endyear)
              thisyr <- allyrs[y]
              gdd <- gdd_all %>% 
                filter(year == thisyr)
              
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
              ind_df$cdl <- 0
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
                
                # Calculate stage-specific degree-days for each cell per day
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
                  negg <- rpois(n = nrow(newrows), lambda = lambda)
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
              
              if(is.nan(voltinism)){ # nobody goes into diapause, unrealistic
                results <- NULL
              }else{
                
                volt_comp <- ind_df %>%
                  filter(diapause == 1 & active == 0) %>% 
                  group_by(numgen) %>% 
                  tally() %>% 
                  tidyr::spread(numgen, n)
                
                names(volt_comp) <- paste("comp", names(volt_comp), sep = "_")
                
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
                
                results <- data.frame(Year = thisyr, cdl_mu, cdl_sd, lambda, voltinism, lost, volt_att, volt_comp, ann_lam)
              }
              
              return(results)
              
            } # close foreach loop
}) # system.time

stopCluster(cl)


# quirks: indiv can decide to diapause at beginning of TA, but still be active because is doesn't make it to complete TA stage before end of year

res <- bind_rows(flatten(outlist))

library(ggplot2)
library(viridis)
plt <- ggplot(res, aes(x = cdl_mu, y = voltinism)) +
  geom_point(aes(color = ann_lam)) +
  scale_x_reverse() +
  scale_color_viridis(begin = 0, end = 1) +
  facet_wrap(cdl_sd~Year, ncol = 4) +
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
  facet_wrap(cdl_sd~Year, ncol = 4) +
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
  facet_wrap(cdl_sd~Year, ncol = 4) +
  theme_bw() +
  ylab("Log(annual lambda)") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated log(Annual Population Growth Rate)", subtitle = "Accounting for lost generations")
plt

# across years


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geom_lam <- res %>% 
  # filter(Year > 2007) %>% 
  group_by(cdl_mu, cdl_sd, lambda) %>% 
  tidyr::complete(Year, cdl_mu, cdl_sd, lambda, fill = list(ann_lam = 0)) %>% 
  summarise(mean_annual_lambda = gm_mean(ann_lam))

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
plt <- ggplot(geom_lam, aes(x = cdl_mu, y = log(mean_annual_lambda), group = lambda, color = lambda)) +
  geom_line() +
  scale_x_reverse() +
  scale_color_viridis(begin = 1, end = 0) +
  facet_wrap(~cdl_sd) +
  theme_bw() +
  ylab("Log(annual lambda)") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated log(Annual Population Growth Rate)", subtitle = "Accounting for lost generations")
plt


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
      dplyr::select(Year, CDL, voltinism, lost) %>% 
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





