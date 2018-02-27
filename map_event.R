# Make pest event maps
# given a degree-day value, plot date it occurs by raster
# could also extract distribution of dates by shapefile polygons


library(sp)
library(rgdal)
library(raster)

library(ggplot2)
library(viridis)
library(mapdata)
library(dplyr)
library(tidyr)
library(geofacet)
theme_set(theme_bw(base_size = 16)) 

# library(gridExtra)
# library(grid)
gdd_cutoff <- 100

source('CDL_funcs.R')
region_param <- "EAST"
years <- c(2014:2017)
gdd_files <- paste0("dailygdd_", years, "_", region_param, ".grd")

REGION <- assign_extent(region_param = region_param)
states <- map_data("state", xlim = c(REGION@xmin, REGION@xmax),
                   ylim = c(REGION@ymin, REGION@ymax), lforce = "e")
names(states)[1:2] <- c("x", "y")
GDD <- brick(gdd_files[1])

template <- GDD[[1]]
template[!is.na(template)] <- 0
template <- crop(template, REGION)

dflist <- list()
for (yr in years){
  filename <- gdd_files[grep(pattern = yr, x = gdd_files, fixed = TRUE)]
  res <- brick(filename)
  # fix NA problem, just assigned as zero
  res <- res + template # template is all zeros and NA
  accum <- calc(res, fun = cumsum)
  accum <- calc(accum, fun = function(x) {min(which(x >= gdd_cutoff))})
  
  df <- as.data.frame(accum, xy=TRUE)
  names(df)[3] <- "doy"
  df$doy[-which(is.finite(df$doy))] <- NA
  df$year <- yr
  dflist[[length(dflist)+1]] <- df
}
resdf <- dplyr::bind_rows(dflist)

resdf <- resdf %>% 
  group_by(year) %>% 
  mutate(date = format(as.Date(doy, origin=as.Date(paste0((year[1] - 1), "-12-31"))),
                       format = "%m-%d"),
         week = lubridate::week(as.Date(doy, origin=as.Date(paste0((year[1] - 1), "-12-31")))))

plt <- ggplot(resdf, aes(x, y)) +
  geom_raster(aes(fill = week)) +
  scale_fill_viridis(na.value = "white") + 
  geom_polygon(data = states, aes(group = group), fill = NA, color = "black", size = .1) +
  facet_wrap(~year, nrow = 2) +
  coord_fixed(1.3)
plt

# extract by sites

sites <- data.frame(ID = c("Corvallis, OR", "Richland, WA", "JB Lewis-McChord, WA", "Palermo, CA",
                           "Ephrata, WA", "Yakima Training Center, WA", "Camp Rilea, OR",
                           "Ft Drum, NY", "West Point, NY", "Kellogg LTER, MI",
                           "The Wilds, OH", "Duluth, MN", "Coeburn, VA", "Mountain Home AFB, ID",
                           "Quantico MCB, VA", "Hanscom AFB, MA", "Ft Bragg, NC",
                           "Ogden, UT", "Buckley AFB, CO", "S Portland, OR",
                           "Sutherlin, OR", "Bellingham, WA", "Wentzville, MO"),
                    x = c(-123.263, -119.283, -122.53, -121.625360, -119.555424, -120.461073,
                          -123.934759, -75.763566, -73.962210, -85.402260, -81.733314,
                          -92.158597, -82.466417, -115.865101, -77.311254, -71.276231,
                          -79.083248, -112.052908, -104.752266, -122.658887,
                          -123.315854, -122.479482,  -90.852),
                    y = c(44.564, 46.275, 47.112, 39.426829, 47.318546, 46.680138,
                          46.122867, 44.055684, 41.388456, 42.404749, 39.829447,
                          46.728247, 36.943103, 43.044083, 38.513995, 42.457068,
                          35.173401, 41.252509, 39.704018, 45.470532,
                          43.387721, 48.756105, 38.816))



dflist <- list()
for (yr in years){
  filename <- gdd_files[grep(pattern = yr, x = gdd_files, fixed = TRUE)]
  res <- brick(filename)
  # fix NA problem, just assigned as zero
  res <- res + template # template is all zeros and NA
  accum <- calc(res, fun = cumsum)
  accum <- calc(accum, fun = function(x) {min(which(x >= gdd_cutoff))})
  
  e <- raster::extract(accum, sites[, c(2:3)], buffer=1000, fun = mean, na.rm = TRUE)
  
  df <- sites
  df$doy <- e
  df$year <- yr
  dflist[[length(dflist)+1]] <- df
}
resdf <- dplyr::bind_rows(dflist)

resdf <- resdf %>% 
  filter(complete.cases(.)) %>% 
  group_by(year) %>% 
  mutate(date = as.Date(doy, origin=as.Date("2000-12-31")))



plt <- ggplot(resdf, aes(x = reorder(ID, date), y = date, color = as.factor(year))) +
  geom_jitter(size = 3, width = .15) +
  scale_y_date(date_breaks = "1 week", date_labels = "%b-%d") +
  coord_flip() +
  ggtitle("Date when 100 degree-days recorded for Galerucella")
plt

ggsave(filename = "Galerucella_100DD_EAST.png", plot = plt, device = "png", width = 15, height = 8, units = "in")


# Map range of GDD at different sites, 100-150DD dates
cutoffday <- 55 # when PRISM observations switch to forecasts

# Add new predictions for current year
preddf <- gdddf %>% 
  filter(complete.cases(.),
         year != 2018) %>% 
  mutate(yday = as.numeric(yday)) %>% 
  arrange(yday) %>% 
  group_by(ID, yday) %>%   
  summarise(DD = mean(DD)) %>% 
  filter(yday > cutoffday) %>% 
  mutate(year = "2018+5yrmean")


# replace predictions after cutoff day
obs <- gdddf %>% 
  select(-x, -y) %>% 
  mutate(yday = as.numeric(yday)) %>% 
  filter(complete.cases(.),
         year == 2018, 
         yday <= cutoffday) %>% 
  mutate(year = "2018+5yrmean")

preds <- bind_rows(preddf, obs)

gdddf$year <- as.factor(as.character(gdddf$year))
levels(gdddf$year)[6] <- "2018+10yrmean"

resdf <- gdddf %>% 
  select(-x, -y) %>% 
  filter(complete.cases(.)) %>% 
  mutate(yday = as.numeric(yday)) %>% 
  bind_rows(preds) %>% 
  arrange(yday) %>% 
  group_by(year, ID) %>% 
  mutate(AccumDD = cumsum(DD)) %>% 
  filter(AccumDD >= 100) %>% 
  filter(AccumDD <= 150) %>% 
  filter(row_number() == 1L | row_number() == n()) %>% 
  mutate(date = as.Date(yday, origin=as.Date("2000-12-31")),
         bound = c("start", "end")) %>% 
  dplyr::select(ID, year, date, bound) %>% 
  tidyr::spread(bound, date) %>% 
  # filter(ID %in% levels(gdddf$ID)[c(1, 3, 4, 5, 10, 12, 15, 16, 17, 21, 22, 25, 26, 28, 30)]) %>%  # East and West
  filter(ID %in% levels(gdddf$ID)[c(1, 3, 5, 10, 15, 17, 21, 22, 25, 26, 30)]) %>%  # Just NW sites
  group_by(ID) %>% 
  mutate(meanstart = mean(start))






theme_set(theme_bw(base_size = 16)) 


plt <- ggplot(resdf, aes(x = reorder(ID, meanstart), ymin = start, ymax = end, group = as.factor(year), color = as.factor(year))) +
  # geom_jitter(size = 3, width = .15) +
  geom_linerange(position = position_dodge(width = .6), size = 2) +
  scale_color_viridis(discrete = TRUE, name = "Year", option = "D") +
  scale_y_date(date_breaks = "1 week", date_labels = "%b-%d") +
  coord_flip() +
  ggtitle("Date range for Galerucella overwintering adult emergence (100-150 degree-days)") +
  xlab("Sites ordered by\nmean start of emergence") +
  guides(color = guide_legend(reverse=TRUE)) 

plt

ggsave(filename = "Galerucella_daterange.png", plot = plt, device = "png", width = 15, height = 8, units = "in")




