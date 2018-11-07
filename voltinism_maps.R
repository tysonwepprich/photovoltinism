# Map raster results from DDRP models
# 3 maps similar to Grevstad & Coop 2015
# - Potential voltinism based on degree-days only
# - Attempted voltinism based on photoperiod cues
# - Mismatch between the two

# 1. Setup -------
# packages, options, functions loaded
library(sp)
library(rgdal)
library(raster)
library(ggplot2)
library(viridis)
library(scales)
# library(mapdata)
library(dplyr)
library(tidyr)
# plotting options
theme_set(theme_bw(base_size = 20)) 
source('CDL_funcs.R')
source('species_params.R')

# 2. User input -----

# input directory with DDRP results
newname <- "APHA_2017_ALL"

# make map template for CRS, extent, etc.
template <- brick(paste0(newname, "/diap_all.grd"))[[1]]
template[!is.na(template)] <- 0

# download state/province boundaries for mapping
if (!file.exists("./src/ref/ne_50m_admin_1_states_provinces_lakes/ne_50m_admin_1_states_provinces_lakes.dbf")){
  download.file(file.path('http://www.naturalearthdata.com/http/',
                          'www.naturalearthdata.com/download/50m/cultural',
                          'ne_50m_admin_1_states_provinces_lakes.zip'), 
                f <- tempfile())
  unzip(f, exdir = "./src/ref/ne_50m_admin_1_states_provinces_lakes")
  rm(f)
}

region <- readOGR("./src/ref/ne_50m_admin_1_states_provinces_lakes", 'ne_50m_admin_1_states_provinces_lakes', encoding='UTF-8')
region <-spTransform(region, CRS(proj4string(template)))
reg.points = fortify(region, region="name_en")
reg.df = left_join(reg.points, region@data, by = c("id" = "name_en"))


# 3. Summary RasterBrick of voltinism ----
# In model_lifecycle.R, when results are processed, each generation
# is split into its own rasterbrick. This section brings all of them
# together for a map of the mean generation number for each day at each cell
# Might be a better way to do this...

# Also, depending on whether SplitMap was used in model_lifecycle.R, 
# these files may be named *_weighted.grd or *_all.grd

returnwd <- getwd()
setwd(newname)
f <-list.files()
rasfiles <- f[grep(pattern = "NumGen", x = f, fixed = TRUE)]
# rasfiles <- rasfiles[grep(pattern = "_weighted.grd", x = rasfiles, fixed = TRUE)]
rasfiles <- rasfiles[grep(pattern = "_all.grd", x = rasfiles, fixed = TRUE)]

ls_index <- stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,1]
ls_index <- as.numeric(gsub(pattern = "NumGen", replacement = "", fixed = TRUE, x = ls_index))

for (day in 1:365){
  
  numgenstack <- mosaic(template, stack(rasfiles[1])[[day]], fun = max, na.rm = TRUE) * ls_index[1]
  for (m in 2:length(rasfiles)){
    numgenstack <- addLayer(numgenstack, ls_index[m] * mosaic(template, stack(rasfiles[m])[[day]], fun = max, na.rm = TRUE))
    
  }
  maxnumgen <- sum(numgenstack, na.rm = TRUE) / 1000 + template
  
  if (!exists("numgen_stack")){
    numgen_stack <- stack(maxnumgen)
  } else {
    numgen_stack <- addLayer(numgen_stack, maxnumgen)
  }
}

# see what is looks like
plot(numgen_stack[[350]])

# save for later
writeRaster(numgen_stack, filename = "NumGen", overwrite = TRUE)

# 4. Processing rasters for maps ----
# load RasterBrick of generations
numgen_stack <- brick("NumGen.grd")[[1:365]]

# Raster prep for plots
# lastday is a raster that constrains the degree-day map
# when used with "stackSelect" function (each cell could have different end of season)
# it could be based on frost, etc. See Future Directions.
# for here, it's just a raster of 365 so that it does not
# constrain the end of season
lastday <- template + 365

# raster of potential number of generations with degree-days only
gddvolt <- stackSelect(calc(numgen_stack, function(x){cummax(x)}), lastday)

# calculate raster of generations constrained by diapause choices
diap_brick <- stack("diap_all.grd")
diap <- stackSelect(diap_brick, lastday)
threshold <- .95 # assumption that if 95% of population in diapause, stop counting new generations
volt <- Cond(diap_brick <= threshold * 1000, numgen_stack, 1) # just assigns diapause cells to one so they don't count as maximum in next line
volt2 <- stackSelect(calc(volt, function(x){cummax(x)}), lastday)

# this accounts for partial generations where some don't go into diapause
# it counts partial generation as decimals, but this isn't ideal
volt3 <- volt2 + (1 - diap/1000)

# deadmask just assigns places where 0 generations complete to NA
# helps with maps so these are blank spaces
dead <- Cond(volt[[365]] == 0, 0, volt[[365]])
deadmask <- Cond(dead > 0, 1, 0)
deadmask[deadmask == 0] <- NA
volt4 <- volt3 * deadmask
gddvolt <- gddvolt * deadmask
# mismatch
mmvolt <- (gddvolt - volt4) * deadmask
plot(mmvolt)


# 5. Pretty maps ----
# TODO: make naming files automatic
# continuous vs discrete color scales not automatic yet either, choose below
map_type <- "mismatch" # "potential", "attempted"
plttitle <- "Aphalara voltinism mismatch with photoperiod cue"
pltname <- "APHA_mismatch_2017_discrete.png"


to_map <- switch(map_type,
                 "mismatch" = mmvolt,
                 "potential" = gddvolt,
                 "attempted" = volt4)
df <- as.data.frame(to_map, xy=TRUE)

names(df)[3] <- "Voltinism"

# DISCRETE OR CONTINUOUS COLOR SCALE?
# discrete voltinism classes for visual
df$Voltinism <- round(df$Voltinism)
maxvolt <- max(df$Voltinism, na.rm = TRUE)
minvolt <- min(df$Voltinism, na.rm = TRUE)
df$Voltinism <- factor(df$Voltinism, levels = c(minvolt:maxvolt))
# change levels manually as desired
levels(df$Voltinism) <- c(levels(df$Voltinism)[1:5], rep("4+", 7))
df <- df %>% 
  filter(!is.na(Voltinism))


# continuous voltinism for visual
df <- df %>% 
  filter(!is.na(Voltinism)) %>% 
  mutate(Voltinism = ifelse(Voltinism >= 4, 4, Voltinism))


# Plot, several options for scaling, choose one
# viridis needs to be told discrete or continuous
# gradient2 only does continuous
# brewer only does discrete classes
theme_set(theme_void(base_size = 18)) 
tmpplt <- ggplot(data = df, aes(x, y, fill = Voltinism)) +
  geom_raster() +
  geom_polygon(data = reg.df, aes(x = long, y = lat, group = group), fill = NA, color = "dark gray", inherit.aes = FALSE, size = .25) +
  # scale_fill_viridis(breaks = seq(0, 10, by = 1), name = "Generations", na.value = "white", begin = 0, end = 1, discrete = TRUE) +
  # scale_fill_gradient2(low=muted("red"), high=muted("blue"), midpoint = 0) +
  scale_fill_brewer(guide = guide_legend(reverse = TRUE), type = "div", palette = "RdBu") +
  ggtitle(plttitle) +
  coord_fixed(1.3, xlim = c(min(df$x), max(df$x)), ylim = c(min(df$y), max(df$y)), expand = FALSE, clip = "on") +
  theme(legend.position = c(0.1, 0.15), plot.title = element_text(hjust = 0.5))
tmpplt
ggsave(pltname,
       plot = tmpplt, device = "png", width = 18, height = 11, units = "in")




# 6. Future directions ----

# what about frost limiting?
# using tmin object from model_lifecycle.R for daily min temperatures
tmin <- crop(aggregate(tminfile[[180:365]], fact = 2, fun = mean, na.rm = TRUE, expand = TRUE), template)

# find hard frosts less than -2C
test <- Cond(tmin <= -2, 1, 0)

# find and plot first frost
first <- calc(test, fun=cumsum)
ffrost <- 179 + which.min(Cond(first == 0, 2, first)) # what I want is cumsum == 1 for first day
ffrost <- Cond(ffrost == 180, 365, ffrost) # set no frost places to 365
plot(ffrost)

# takes voltinism raster and finds the generation present
# at each cell when first frost occurs at the cell
numgen_frost <- stackSelect(numgen_stack, ffrost)




# this code is a kludge to add a pre-diapause requirement for aphalara
# need tmin and tmax rasters to calculate gdd raster
# start from end of season
# add pre-diapause gdd requirement for aphalara post-model run
for (index in 365:165){
  tmin <- crop(aggregate(tminfile[[index]], fact = 2, fun = mean, na.rm = TRUE, expand = TRUE), template)
  tmax <- crop(aggregate(tmaxfile[[index]], fact = 2, fun = mean, na.rm = TRUE, expand = TRUE), template)
  dd_tmp <- overlay(tmax, tmin, template + 6.9, template + 32, fun = TriDD)
  
  if (!exists("dd_stack")){
    dd_stack <- stack(dd_tmp)
  } else {
    dd_stack <- addLayer(dd_stack, dd_tmp)
  }
}

last <- calc(dd_stack, fun=cumsum)
lastday <- 365 - which.max(Cond(last > 75, 1, last))
