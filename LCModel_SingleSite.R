# Single-site lifecycle model
# Adapted from DDRP raster-based model

# Features: 
# 1. Could allow individual variation beyond overwintering emergence, since computing not limited
# 2. Transfers between sites to model common garden / reciprocal transplants

#####
# packages, options, functions loaded
#####
library(sp)
library(rgdal)
library(raster)
library(lubridate)
library(mixsmsn)
library(dplyr)
library(ggplot2)

source('CDL_funcs.R') # load collection of functions for this model

# prism_path <- "prismDL/2014"

sites <- read.csv("data/GCA_modeling_sites.csv", header = TRUE) %>% 
  mutate(ID = as.character(ID)) %>% 
  arrange(ID)
# sites <- sites[c(5,9), ]
sites <- sites[-c(2, 14, 18, 20, 27), ]


# west_sites <- sites[c(1, 3, 4, 6, 7, 20, 21, 22), ] %>% droplevels.data.frame()
# east_sites <- 
# can_sites <- 
#####
#input parameters
#####
# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model extext
start_doy  <- 1
end_doy    <- 365
region_param <- "NW_SMALL"
# life cycle parameters
stgorder   <- c("OA","E","L","P","A","F")
photo_sens <- 3 #c(-1, 3) # integer life stages for now
CDL_mu        <- 16.23 # 15.52
model_CDL  <- 1 # if 1, model photoperiod decision
CDL_log    <- 1 # if 1, model CDL from logistic regression results
owstage    <- "OA"
# Logistic regression photoperiod
# from Fritzi's lab data, two different populations
# coefs <- c(-56.9745, 3.5101) # Northern
coefs <- c(-60.3523, 3.8888) # Southern
# Degree day thresholds
# LDT = lower development threshold, temp at which growth = 0 (using PRISM tmean)
eggLDT     <- 10
eggUDT     <- 37.8  
larvaeLDT  <- 10
larvaeUDT  <- 37.8  
pupaeLDT   <- 10
pupaeUDT   <- 37.8
adultLDT   <- 10 
adultUDT   <- 37.8
# Degree day requirements for life stages
# DD = degree days, number of cumulative heat units to complete that lifestage
OWadultDD_mu  <- 167.6 #based on Oregon counts
#or 108 based on McAvoy 1997 phenology of 3 years
eggDD_mu = 93.3
larvaeDD_mu = 136.4 # 46.9 + 45.8 + 43.7 instars
pupDD_mu = 137.7 
adultDD_mu = 125.9 #time to oviposition
# GDD data and calculation
gdd_data <- "load" # c("load", "calculate")
gdd_file <- "dailygdd_2016_NW_SMALL.grd"
calctype   <-"triangle"
# introducing individual variation, tracked with simulation for each substage
vary_indiv <- 1 # turn on indiv. variation
if (vary_indiv == 1){
  nsim <- 11 # number of substages
}else{
  nsim <- 1
}

#####
# derived parameters
#####
REGION <- assign_extent(region_param = region_param)

nday <- length(start_doy:end_doy)

# Substages with skewed t
arg1 = list(mu = 97.94, sigma2 = 2241.7, shape = 3.92, nu = 9.57)
x = seq(50, 350, length.out = 1000)
y = mixsmsn:::dt.ls(x, loc = arg1$mu, sigma2 = arg1$sigma2, shape = arg1$shape, nu = arg1$nu)
inputdist <- data.frame(x = x, y = y) %>% 
  arrange(x) %>% 
  mutate(CDF = cumsum(y/sum(y)))
substages <- SubstageDistrib(dist = inputdist, numstage = nsim, perc = .99)
# To get observations to fit for overwinter adults and F1 eggs, 
# overwinter pre-oviposition period is only 50 deg days
# substages$means <- substages$means


# template <- raster(files[1])
# template <- crop(raster(files[1]), REGION)
# template[!is.na(template)] <- 0


# # preload PRISM data into RasterBrick
# if (gdd_data == "load"){
#   GDD <- brick(gdd_file)
#   DT_same <- TRUE
# }

if (vary_indiv == 1){
  eggDD = eggDD_mu
  larvaeDD = larvaeDD_mu
  pupDD = pupDD_mu
  adultDD = adultDD_mu
  OWadultDD = substages[, 1] # if using SubstageDistrib function
  cdl_b0 <- coefs[1]
  cdl_b1 <- coefs[2]
  # CDL <- CDL_mu
  params <- data.frame(eggDD, larvaeDD, pupDD, adultDD, OWadultDD, cdl_b0, cdl_b1)
}else{ # doesn't make sense to do this with foreach if no variation, though!
  eggDD = eggDD_mu
  larvaeDD = larvaeDD_mu
  pupDD = pupDD_mu
  adultDD = adultDD_mu
  OWadultDD = OWadultDD_mu
  CDL <- CDL_mu
  params <- data.frame(eggDD, larvaeDD, pupDD, adultDD, OWadultDD, CDL)
}

# new function for photoperiod
SitePhoto <- function(sites, doy, perc_twilight){
  p <- perc_twilight * 6 / 100
  hours <- photoperiod(sites$y, doy, p)
}


#####
#Bring together West, East, Canadian GDD
#####
gdd_year <- 2015
gdd_files <- c("dailygdd_2015_NW_SMALL.grd", "dailygdd_2015_EAST.grd")

site_gdd1 <- extract(brick(gdd_files[1]), sites[, c("x", "y")])
site_gdd1 <- data.frame(site_gdd1)
site_gdd1$ID <- sites$ID
site_gdd2 <- extract(brick(gdd_files[2]), sites[, c("x", "y")])
site_gdd2 <- data.frame(site_gdd2)
site_gdd2$ID <- sites$ID

site_gdd_tidy1 <- site_gdd1 %>% 
  group_by(ID) %>% 
  tidyr::gather(key = "DOY", value = "GDD", layer.1:layer.365) %>% 
  ungroup() %>% 
  dplyr::mutate(DOY = as.numeric(gsub(pattern = "layer.", replacement = "", x = .$DOY))) %>% 
  filter(complete.cases(.))
site_gdd_tidy2 <- site_gdd2 %>% 
  group_by(ID) %>% 
  tidyr::gather(key = "DOY", value = "GDD", layer.1:layer.365) %>% 
  ungroup() %>% 
  dplyr::mutate(DOY = as.numeric(gsub(pattern = "layer.", replacement = "", x = .$DOY))) %>% 
  filter(complete.cases(.))
site_gdd_tidy <- bind_rows(site_gdd_tidy1, site_gdd_tidy2) %>% droplevels.data.frame()
  

if (exists("layer.366", site_gdd_tidy)){
  site_gdd_tidy$layer.366 <- NULL
}

# Add Canadian sites from Daymet
can_gdd <- readRDS("data/gdd_canada_sites.RDS") %>% 
  filter(Year == gdd_year, DOY != 366) %>% 
  dplyr::select(-Year) %>% 
  mutate(ID = as.character(ID))

site_gdd_tidy <- bind_rows(site_gdd_tidy, can_gdd) %>% 
  filter(ID %in% sites$ID) %>% 
  arrange(DOY, ID) %>% 
  group_by(ID) %>% 
  mutate(AccumDD = cumsum(GDD)) %>% 
  ungroup()
  
#####
# Add transfer to common garden
#####
common_garden <- 0
transfer_DD <- 120
transfer_site <- "Corvallis, OR"


#####
# Run model
#####
results <- list()
# old nested foreach loop
for(sim in 1:nsim){
  # Initialize all tracking rasters as zero with the template
  eggDD <- params$eggDD[sim]
  larvaeDD <- params$larvaeDD[sim]
  pupDD <- params$pupDD[sim]
  adultDD <- params$adultDD[sim]
  OWadultDD <- params$OWadultDD[sim]
  # CDL <- params$CDL[sim]
  cdl_b0 <- params$cdl_b0[sim]
  cdl_b1 <- params$cdl_b1[sim]
  
  template <- rep(0, nrow(sites)) 
  DDaccum <- template
  TotalGDD <- template
  transferred <- template # common garden feature, start in own site
  Lifestage <- rep(-1, nrow(sites))  # Need for OW designated as stage -1
  NumGen <- template
  #Lifestage: [0] = egg, [1] = larvae, [2] = pupae, [3] = adult
  
  #### Init Lifestage tracking
  if (owstage == "OE") {
    LSOW0 <- Lifestage == -1 
  } else if (owstage == "OL") {
    LSOW1 <- Lifestage == -1 
  } else if (owstage == "OP") {
    LSOW2 <- Lifestage == -1 
  } else if (owstage == "OA") {
    LSOW3 <- Lifestage == -1 
  }
  LS0 <- Lifestage == 0
  LS1 <- Lifestage == 1
  LS2 <- Lifestage == 2
  LS3 <- Lifestage == 3
  LS4 <- Lifestage == 4
  
  # loop over days
  sublist <- start_doy:end_doy
  outlist <- list()
  diaplist <- list()
  for (d in sublist) {
    print(d)
    index <- which(sublist == d)

    # Common Garden transfer
    # all sites
      tmpGDD <- site_gdd_tidy %>% 
        filter(DOY == d) %>% 
        dplyr::select(GDD) %>% 
        unlist()
      tmpGDD[which(is.na(tmpGDD))] <- 0 # in case of missing GDD (could interpolate)
      # transfer site
      tmpGDDtr <- site_gdd_tidy %>% 
        filter(DOY == d, ID == transfer_site) %>% 
        dplyr::select(GDD) %>% 
        unlist()
      tmpGDDtr[which(is.na(tmpGDDtr))] <- 0
      tmpGDD <- Cond(transferred == 1, tmpGDDtr, tmpGDD)
      
      TotalGDD <- TotalGDD + tmpGDD

      # photoperiod for this day across raster
      if (model_CDL == 1){
        doy <- d
        photo <- SitePhoto(sites, doy, perc_twilight = 25)
        phototr <- SitePhoto(sites[which(sites$ID == transfer_site), ], doy, perc_twilight = 25)
        photo <- Cond(transferred == 1, phototr, photo)
      }
      
      # transferred if threshold reached, otherwise no transfer
      if(common_garden == 1){
        transferred <- Cond(TotalGDD >= transfer_DD, 1, 0)
      }
    
    #### Loop through Stages; order of stages now read from SPP param file ####
    
    for (i in stgorder) {  # Handle stages in the model
      ####  MAIN STEPS FOR EGG STAGE ####
      if (i == "E" | i == "OE") {   # Egg Stage
        dd0tmp <- tmpGDD
        if (i == "OE") { 
          dd0 <- dd0tmp * LSOW0 
        } else if (i == "E") { 
          dd0 <- dd0tmp * LS0 
        }
        DDaccum <- DDaccum + dd0
        
        #Calculate lifestage progression: for dd accum in correct lifestage, is it >= egg DD threshold?
        if (i == "OE") {
          progressOW0 <- (DDaccum * LSOW0) >= OWeggDD
          Lifestage <- Cond(LSOW0 == 1 & progressOW0 == 1, 1, Lifestage)
          #Reset the DDaccum cells to zero for cells that progressed to next lifestage
          DDaccum <- Cond(progressOW0 == 1,(DDaccum - OWeggDD) * LSOW0, DDaccum)
        } else if (i == "E") {
          progress0 <- (DDaccum * LS0) >= eggDD
          Lifestage <- Cond(LS0 == 1 & progress0 == 1, 1, Lifestage)
          #Reset the DDaccum cells to zero for cells that progressed to next lifestage
          DDaccum <- Cond(progress0 == 1,(DDaccum - eggDD) * LS0, DDaccum)
        }

        ####  MAIN STEPS FOR LARVAL STAGE ####
      } else if (i == "L" | i == "OL") {  # Larval Stage
        dd1tmp <- tmpGDD
        if (i == "OL") { 
          dd1 <- dd1tmp * LSOW1
        }else if (i == "L") { 
          dd1 <- dd1tmp * LS1
        }
        DDaccum <- DDaccum + dd1
        
        #Calculate lifestage progression: for dd accum in correct lifestage, is it >= egg DD threshold?
        if (i == "OL") {
          progressOW1 <- (DDaccum * LSOW1) >= OWlarvaeDD
          Lifestage <- Cond(LSOW1 == 1 & progressOW1 == 1, 2, Lifestage)
          #Reset the DDaccum cells to zero for cells that progressed to next lifestage
          DDaccum <- Cond(progressOW1 == 1,(DDaccum - OWlarvaeDD) * LSOW1, DDaccum)
        } else if (i == "L") {
          progress1 <- (DDaccum * LS1) >= larvaeDD
          Lifestage <- Cond(LS1 == 1 & progress1 == 1, 2, Lifestage)
          DDaccum <- Cond(progress1 == 1,(DDaccum - larvaeDD) * LS1, DDaccum)
        }
        
        ####  MAIN STEPS FOR PUPAL STAGE ####
      } else if (i == "P" | i == "OP") {   # Pupal Stage
        dd2tmp <- tmpGDD
        if (i == "OP") { 
          dd2 <- dd2tmp * LSOW2
        } else if (i == "P") { 
          dd2 <- dd2tmp * LS2
        }
        DDaccum <- DDaccum + dd2
        
        #Calculate lifestage progression: for dd accum in correct lifestage, is it >= pupae DD threshold?
        if (i == "OP") {
          progressOW2 <- (DDaccum * LSOW2) >= OWpupDD
          Lifestage <- Cond(LSOW2 == 1 & progressOW2 == 1, 3, Lifestage)
          #Reset the DDaccum cells to zero for cells that progressed to next lifestage
          DDaccum <- Cond(progressOW2 == 1,(DDaccum - OWpupDD) * LSOW2, DDaccum)
        } else if (i == "P") {
          progress2 <- (DDaccum * LS2) >= pupDD
          Lifestage <- Cond(LS2 == 1 & progress2 == 1, 3, Lifestage)
          DDaccum <- Cond(progress2 == 1,(DDaccum - pupDD) * LS2, DDaccum)
        }
        
        ####  MAIN STEPS FOR ADULT STAGE ####
      } else if (i == "A" | i == "OA") {  # Adult stage, or time to 50% oviposition
        dd3tmp <- tmpGDD
        if (i == "OA") { 
          dd3 <- dd3tmp * LSOW3 
        } else if (i == "A") { 
          dd3 <- dd3tmp * LS3 
        }
        DDaccum <- DDaccum + dd3
        
        #Calculate lifestage progression: for dd accum in correct lifestage, is it >= adult DD threshold?
        if (i == "OA") {
          progressOW3 <- (DDaccum * LSOW3) >= OWadultDD
          DDaccum <- Cond(progressOW3 == 1,(DDaccum - OWadultDD) * LSOW3, DDaccum)
          progressOW3[is.na(progressOW3)] <- template[is.na(progressOW3)]
          Lifestage <- Cond(LSOW3 == 1 & progressOW3 == 1, 0, Lifestage)

        } else if (i == "A") {
          progress3 <- (DDaccum * LS3) >= adultDD
          #Reset the DDaccum cells to zero for cells that progressed to next lifestage
          DDaccum <- Cond(progress3 == 1,(DDaccum - adultDD) * LS3, DDaccum)
          #Remove masking effect from progress counter so it doesn't perpetuate through NumGen raster
          Lifestage <- Cond(LS3 == 1 & progress3 == 1,0, Lifestage)

        }

        ####  MAIN STEPS FOR END OF EACH DAY ####
      } else if (i == "F") { # end of the day placeholder 
        
        if (model_CDL == 1){

          sens_mask <- Cond(Lifestage %in% photo_sens, 1, 0)
          prop_diap <- 1 - exp(cdl_b0 + cdl_b1 * photo) /
            (1 + exp(cdl_b0 + cdl_b1 * photo))
          tmpLS4 <- Cond(sens_mask == 1, prop_diap, LS4)
          # need to account for prop_diap declining up until solstice
          # only let proportion diapausing increase over season
          LS4 <- Cond(prop_diap < LS4, LS4, tmpLS4)

        }
        
        # update lifestage masks
        if (owstage == "OE") {
          LSOW <- LSOW0 <- Lifestage == -1
        } else if (owstage == "OL") {
          LSOW <- LSOW1 <- Lifestage == -1
        } else if (owstage == "OP") {
          LSOW <- LSOW2 <- Lifestage == -1
        } else if (owstage == "OA") {
          LSOW <- LSOW3 <- Lifestage == -1
        }
        LS0 <- Lifestage == 0
        LS1 <- Lifestage == 1
        LS2 <- Lifestage == 2
        LS3 <- Lifestage == 3
        
      } # close F
    } # close stgorder
    outlist[[index]] <- data.frame(Lifestage)
    diaplist[[index]] <- data.frame(LS4)
  } # close day
  outdf <- bind_cols(outlist) 
  names(outdf) <- paste("DOY", 1:365, sep = "_")
  outdf$ID <- sites$ID
  outdf <- outdf %>% 
    tidyr::gather(DOY, Lifestage, DOY_1:DOY_365) %>% 
    mutate(DOY = as.numeric(stringr::str_split_fixed(string = DOY, pattern = "_", n = 2)[, 2]),
           sim = sim) 
  
  diap <- bind_cols(diaplist)
  names(diap) <- paste("DOY", 1:365, sep = "_")
  diap$ID <- sites$ID
  diap <- diap %>% 
    tidyr::gather(DOY, Diap_perc, DOY_1:DOY_365) %>% 
    mutate(DOY = as.numeric(stringr::str_split_fixed(string = DOY, pattern = "_", n = 2)[, 2]),
           sim = sim)
  
  outdf <- left_join(outdf, diap)
  results[[sim]] <- outdf
} #close for loop

allresults <- bind_rows(results)

weighted_lifestage <- substages %>% 
  mutate(sim = 1:nrow(substages)) %>% 
  right_join(allresults) %>% 
  group_by(ID, DOY, Lifestage) %>% 
  summarise(Relprop = sum(weights) / sum(substages$weights))


weighted_diapause <- substages %>% 
  mutate(sim = 1:nrow(substages)) %>% 
  right_join(allresults) %>% 
  dplyr::select(sim, weights, ID, DOY, Diap_perc) %>% 
  group_by(ID, DOY) %>% 
  summarise(Prop_diapause = sum(weights * Diap_perc) / sum(substages$weights))

allresults2 <- weighted_lifestage %>% 
  ungroup() %>% 
  mutate(Lifestage = as.factor(Lifestage)) %>% 
  tidyr::complete(DOY, ID, Lifestage, fill = list(Relprop = 0)) %>% 
  left_join(site_gdd_tidy)

allresults <- allresults2 %>% 
  left_join(weighted_diapause) %>% 
  mutate(Relprop = Relprop * (1 - Prop_diapause),
         Date = as.Date(DOY, origin=as.Date(paste0((gdd_year - 1), "-12-31"))))
diap <- weighted_diapause %>%
  ungroup() %>%
  left_join(allresults[, c("ID", "DOY", "AccumDD", "Date")]) %>%
  distinct()

# transfer_date <- site_gdd_tidy %>%
#   group_by(ID) %>%
#   filter(AccumDD >= transfer_DD) %>%
#   arrange(DOY) %>%
#   filter(1:n() == 1) %>%
#   mutate(Date = as.Date(DOY, origin=as.Date(paste0((gdd_year - 1), "-12-31"))))

levels(allresults$Lifestage) <- c("Overwinter", "Egg", "Larva", "Pupa", "Adult")

# plot to check results
plt <- ggplot(allresults, aes(x = Date, y = Relprop, group = Lifestage, color = Lifestage)) +
  geom_line() +
  scale_x_date(date_breaks = "2 month", date_labels = "%b") +
  facet_wrap(~ID, ncol = 5) +
  ggtitle("No transfer, 2015, Southern CDL") +
  theme_bw() +  
  geom_line(data = diap, aes(x = Date, y = Prop_diapause, group = ID, color = "diapause")) 
  # geom_vline(data = transfer_date, aes(xintercept = Date))
plt


ggsave("GCAsites_notransfer_2015_SCDL.png",
       plot = plt, device = "png", width = 12, height = 8, units = "in")


plt <- ggplot(site_gdd_tidy, aes(x = DOY, y = GDD)) +
  geom_point() +
  facet_wrap(~ID, ncol = 5) +
  theme_bw()
plt

plt <- ggplot(test2, aes(x = Date, y = Relprop, group = Lifestage, color = Lifestage)) +
  geom_line() +
  # geom_vline(xintercept = 120) +
  # scale_x_date(date_breaks = "2 month", date_labels = "%b") +
  # facet_wrap(~ID, ncol = 5) +
  theme_bw()
plt


test <- left_join(allresults, site_gdd_tidy)
plt <- ggplot(test, aes(x = AccumDD, y = Lifestage)) +
  geom_line() +
  facet_wrap(~ID, ncol = 5) +
  theme_bw()
plt

justgdd <- test %>% 
  dplyr::select(-Lifestage, sim, Diap_perc) %>% 
  distinct()
plt <- ggplot(justgdd, aes(x = DOY, y = AccumDD)) +
  geom_point() +
  facet_wrap(~ID, ncol = 5) +
  geom_hline(yintercept = 100) +
  theme_bw()
plt
