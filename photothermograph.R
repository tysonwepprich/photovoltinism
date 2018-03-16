
# This code will:
# 1. Download PRISM data if needed (unfortunately, does entire CONUS for each year).
# 2. Calculate accumulated degree-days for each day/year at multiple sites.
# 3. Calculate photoperiod for each day/year at multiple sites.
# 4. Plot photothermograph with species critical photoperiod
# 5. Expand to individual variation in critical photoperiod and/or required degree-days.

library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(lubridate)
library(sp)
library(rgdal)
library(raster)
library(purrr)
library(stringr)
library(prism)

# Input development thresholds
LDT <- 10
UDT <- 37

# Input where your PRISM downloaded data are
# base_path <- "/data/PRISM/"
base_path <- "prismDL"

# Input years to analyze
years <- c(2013:2018)

# Input sites with ID, x, y columns
sites <- read.csv("data/GCA_modeling_sites.csv", header = TRUE)


# Utility functions
source('CDL_funcs.R')

# # PRISM download
# # using ropensci 'prism' package to access webservice
# # downloads entire CONUS, so files are large
# 
# startdate <- "2018-02-28"
# enddate <- "2018-03-14"
# yr <- year(ymd(startdate))
# 
# # need to set download directory
# options(prism.path = paste("prismDL", yr, sep = "/"))
# 
# get_prism_dailys("tmin", minDate = startdate, 
#                  maxDate = enddate, keepZip = FALSE)
# get_prism_dailys("tmax", minDate = startdate, 
#                  maxDate = enddate, keepZip = FALSE)
# 
# 
# 
# # Site based growing degree-days from PRISM
# # Also can do this on rasters to make maps
# outlist <- list()
# for (yr in years){
#   prism_path <- paste(base_path, yr, sep = "/")
#   
#   # Search pattern for PRISM daily temperature grids. Load them for processing.
#   # TMIN
#   pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
#   tminfiles <- list.files(path = prism_path, pattern=pattern, 
#                           all.files=FALSE, full.names=TRUE, recursive = TRUE)
#   
#   # PRISM has different names for files with different quality
#   # this function selects according to hierarchy
#   tminfiles <- ExtractBestPRISM(tminfiles, yr)
#   
#   r <- raster(tminfiles[1])
#   tminstack <- stack(r)
#   tminstack@layers <- sapply(tminfiles, function(x) { r@file@name=x; r } ) 
#   
#   # TMAX
#   pattern = paste("(PRISM_tmax_)(.*)(_bil.bil)$", sep="")
#   tmaxfiles <- list.files(path = prism_path, pattern=pattern, 
#                           all.files=FALSE, full.names=TRUE, recursive = TRUE)
#   tmaxfiles <- ExtractBestPRISM(tmaxfiles, yr)
#   
#   r <- raster(tmaxfiles[1])
#   tmaxstack <- stack(r)
#   tmaxstack@layers <- sapply(tmaxfiles, function(x) { r@file@name=x; r } ) 
#   
#   # by site
#   tmin <- extract(tminstack, sites[, c("x", "y")])
#   tmax <- extract(tmaxstack, sites[, c("x", "y")])
#   
#   GDD <- data.frame(TriDD(tmax, tmin, LDT, UDT))
#   names(GDD) <- paste("yday", 1:ncol(GDD), sep = "_")
#   gdd <- sites %>% 
#     bind_cols(GDD) %>% 
#     tidyr::gather(key = yday, value = DD, -ID, -x, -y) %>% 
#     mutate(yday = unlist(map(stringr::str_split(string = yday, pattern = coll("_"), 2), 2)),
#            year = yr)
#   outlist[[length(outlist)+1]] <- gdd 
# }
# 
# gdddf <- bind_rows(outlist)
# 
# saveRDS(gdddf, "gdddf.rds")

gdddf <- readRDS("gdddf.rds")


# Map range of GDD at different sites, 100-150DD dates
cutoffday <- 73 # day of year when PRISM observations switch to forecasts

# Add new predictions for current year
preddf <- gdddf %>% 
  filter(complete.cases(.),
         year != 2018) %>% 
  mutate(yday = as.numeric(yday)) %>% 
  arrange(yday) %>% 
  group_by(ID, yday, x, y) %>%   
  summarise(DD = mean(DD)) %>% 
  filter(yday > cutoffday) %>% 
  mutate(year = "2018+5yrmean")

# replace predictions after cutoff day
obs <- gdddf %>% 
  mutate(yday = as.numeric(yday)) %>% 
  filter(complete.cases(.),
         year == 2018, 
         yday <= cutoffday) %>% 
  mutate(year = "2018+5yrmean")

preds <- bind_rows(preddf, obs)

gdddf$year <- as.factor(as.character(gdddf$year))
levels(gdddf$year)[6] <- "2018+10yrmean"

resdf <- gdddf %>% 
  filter(complete.cases(.)) %>% 
  mutate(yday = as.numeric(yday)) %>% 
  bind_rows(preds) %>% 
  arrange(yday) %>% 
  group_by(year, ID) %>% 
  mutate(AccumDD = cumsum(DD)) %>% 
  group_by(ID) %>% 
  mutate(photo = photoperiod(lat = y, doy = yday),
         photo2 = geosphere::daylength(lat = y, doy = yday)) %>% 
  # filter(ID %in% levels(gdddf$ID)[c(1, 3, 4, 5, 10, 12, 15, 16, 17, 21, 22, 25, 26, 28, 30)]) # East and West
  filter(ID %in% levels(gdddf$ID)[c(1, 3, 5, 10, 15, 17, 21, 22, 25, 26, 30)])  # Just NW sites

maxgdd <- max(resdf$AccumDD)
owgdd <- 100
gengdd <- 493
sens_gdd <- seq(from = owgdd, to = maxgdd, by = gengdd)

theme_set(theme_bw(base_size = 16)) 


plt <- ggplot(resdf, aes(x = AccumDD, y = photo, group = as.factor(year), color = as.factor(year))) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, name = "Year", option = "D") +
  guides(color = guide_legend(reverse=TRUE)) +
  geom_hline(yintercept = 16.2, linetype = "dashed") +
  geom_vline(xintercept = sens_gdd, linetype = "dotted") +
  facet_wrap(~ID, ncol = 3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.09)) +
  ggtitle("Galerucella photothermographs with Northern-adapted critical photoperiod") +
  xlab("Accumulated degree-days") +
  ylab("Photoperiod (1.5 degrees twilight included)")
plt

ggsave(filename = "Galerucella_phototherm_north.png", plot = plt,
       device = "png", width = 16, height = 10, units = "in")


# output data to csv for import into Google Calendar
outdf <- resdf %>% 
  filter(year == "2018+5yrmean") %>% 
  dplyr::select(ID, start, end) %>% 
  mutate(start = start + duration(17, units = "years"),
         end = end + duration(17, units = "years"))
names(outdf) <- c("Subject", "Start Date", "End Date")
write.csv(outdf, file = "GoogleEvents.csv", row.names = FALSE)



