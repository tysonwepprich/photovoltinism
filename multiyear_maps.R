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
dat <- readRDS("APHA_output_2016to2018_NWSMALL_N.rds")
yr <- 2016

# weather_path <- "/home/macav2metdata/IPSL_rcp85/"
weather_path <- "data/PRISM" # PRISM data on grub server (needs to have stable files downloaded)
# weather_data_source <- "macav2" # could also have daymet or macav2
weather_data_source <- "prism"


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


# Derived parameters -----
REGION <- assign_extent(region_param = region_param)

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
  
  tmin <- shift(tmin, x = -360)
  tmax <- shift(tmax, x = -360)
  
  geo_template <- crop(tmin[[1]], REGION)
  # template <- crop(aggregate(raster(files[1]), fact = 2), REGION)
  geo_template[!is.na(geo_template)] <- 0
}

template <- geo_template

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


# summarise years of MACA simulations
br <- brick(replicate(3, template))
potential <- attempt <- mismatch <- diapause <- gdd <- br

for (y in 1:3){
  attem <- setValues(br, dat[, , 6, y])
  poten <- setValues(br, dat[, , 7, y])
  weivolt <- setValues(br, dat[, , 8, y])
  diap <- setValues(br, dat[, , 9, y])
  VoltDiap <- function(x, y){
    ss <- (subset(x, 2:nlayers(x)) - subset(x, 1:(nlayers(x)-1))) * (1 - y)[[2:nlayers(x)]]
    tmp <- calc(ss, fun = cumsum)
    return(tmp)
  }
  voltdiap <- VoltDiap(attem, diap)
  potential <- setValues(potential, values = getValues(poten[[53]]), layer = y)
  attempt <- setValues(attempt, values = getValues(voltdiap[[52]]), layer = y)
  mismatch <- setValues(mismatch, values = getValues(voltdiap[[52]] - poten[[53]]), layer = y)
  diapause <- setValues(diapause, values = getValues(diap[[53]]), layer = y)
  # gdd <- setValues(gdd, values = getValues(dd[[53]]), layer = y)
  
}

avg_pot <- calc(potential, mean)
avg_att <- calc(attempt, mean)
avg_mm <- calc(mismatch, mean)
avg_diap <- calc(diapause, mean)
# avg_gdd <- calc(gdd, mean)

worst_diap <- calc(diapause, min)
worst_mm <- calc(diapause, min)

# # Raster prep for plots
# # lastday is a raster that constrains the degree-day map
# # when used with "stackSelect" function (each cell could have different end of season)
# # it could be based on frost, etc. See Future Directions.
# # for here, it's just a raster of 365 so that it does not
# # constrain the end of season
# lastday <- template + 365
# 
# # raster of potential number of generations with degree-days only
# gddvolt <- stackSelect(calc(numgen_stack, function(x){cummax(x)}), lastday)
# 
# # calculate raster of generations constrained by diapause choices
# diap_brick <- stack("diap_all.grd")
# diap <- stackSelect(diap_brick, lastday)
# threshold <- .95 # assumption that if 95% of population in diapause, stop counting new generations
# volt <- Cond(diap_brick <= threshold * 1000, numgen_stack, 1) # just assigns diapause cells to one so they don't count as maximum in next line
# volt2 <- stackSelect(calc(volt, function(x){cummax(x)}), lastday)
# 
# # this accounts for partial generations where some don't go into diapause
# # it counts partial generation as decimals, but this isn't ideal
# volt3 <- volt2 + (1 - diap/1000)
# 
# # deadmask just assigns places where 0 generations complete to NA
# # helps with maps so these are blank spaces
# dead <- Cond(volt[[365]] == 0, 0, volt[[365]])
# deadmask <- Cond(dead > 0, 1, 0)
# deadmask[deadmask == 0] <- NA
# volt4 <- volt3 * deadmask
# gddvolt <- gddvolt * deadmask
# # mismatch
# mmvolt <- (gddvolt - volt4) * deadmask
# plot(mmvolt)


# 5. Pretty maps ----
# TODO: make naming files automatic
# continuous vs discrete color scales not automatic yet either, choose below
species <- "Aphalara itadori"
map_type <- "attempted" #"diapause" # "attempted" # "potential", "mismatch"
if (map_type == "diapause"){
  plttitle <- paste(species, map_type, "proportion with photoperiodic cue", sep = " ")
  pltsub <- paste(min(yr), "to", max(yr), "worst year", sep = " ")
}else{
  if (map_type == "potential"){
    plttitle <- paste(species, map_type, "voltinism without photoperiodic cue", sep = " ")
    pltsub <- paste(min(yr), "to", max(yr), "mean", sep = " ")
    
  }else{
    plttitle <- paste(species, map_type, "voltinism with photoperiodic cue", sep = " ")
    pltsub <- paste(min(yr), "to", max(yr), "mean", sep = " ")
  }
}
pltname <- paste(species, map_type, biotype, ".png", sep = "_") 


to_map <- switch(map_type,
                 "mismatch" = avg_mm,
                 "potential" = avg_pot,
                 "attempted" = avg_att,
                 "diapause" = avg_diap)
df <- as.data.frame(to_map, xy=TRUE)

names(df)[3] <- "Voltinism"

# DISCRETE OR CONTINUOUS COLOR SCALE?
# discrete voltinism classes for visual
df$Voltinism <- round(df$Voltinism)
maxvolt <- max(df$Voltinism, na.rm = TRUE)
minvolt <- min(df$Voltinism, na.rm = TRUE)
df$Voltinism <- factor(df$Voltinism, levels = c(minvolt:maxvolt))
# change levels manually as desired
# levels(df$Voltinism) <- c(levels(df$Voltinism)[1:6], rep("6+", 2))
df <- df %>% 
  filter(!is.na(Voltinism))


# continuous voltinism for visual
df <- df %>%
  filter(!is.na(Voltinism)) %>%
  mutate(Voltinism = ifelse(Voltinism >= 5, 5, Voltinism))


# Plot, several options for scaling, choose one
# viridis needs to be told discrete or continuous
# gradient2 only does continuous
# brewer only does discrete classes
pts <- data.frame(x = -123.2326, y = 44.8478)

theme_set(theme_void(base_size = 20))
tmpplt <- ggplot(data = df, aes(x, y, fill = Voltinism)) +
  geom_raster() +
  geom_polygon(data = reg.df, aes(x = long, y = lat, group = group), fill = NA, color = "black", inherit.aes = FALSE, size = 1) +
  # scale_fill_viridis(breaks = seq(0, 7, by = 1), name = ifelse(map_type == "diapause", "% Diapause",  "Generations"), na.value = "white", begin = 0, end = ifelse(map_type == "diapause", 1, maxvolt/7), discrete = ifelse(map_type == "diapause", FALSE, TRUE)) +
  # scale_fill_gradient2(low=muted("red"), high=muted("blue"), midpoint = 0) +
  scale_fill_viridis(na.value = "white", begin = 0, end = maxvolt / 3.58) +
  # scale_fill_brewer(guide = guide_legend(reverse = TRUE), type = "div", palette = "RdBu") +
  geom_point(data = pts, aes(x, y), inherit.aes = FALSE, size = 2) +
  # ggtitle(plttitle, subtitle = pltsub) +
  coord_fixed(1.3, xlim = c(min(df$x), max(df$x)), ylim = c(min(df$y), max(df$y)), expand = FALSE, clip = "on") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
tmpplt
ggsave(pltname,
       plot = tmpplt, device = "png", width = 6, height = 10, units = "in")




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
