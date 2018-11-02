# Summarize and map raster results from DDRP models

# 1. Setup -------
# packages, options, functions loaded
library(sp)
library(rgdal)
library(raster)
library(ggplot2)
library(viridis)
library(mapdata)
library(dplyr)
library(tidyr)
# library(geofacet) # didn't install correctly
library(gridExtra)
library(grid)
# plotting options
theme_set(theme_bw(base_size = 20)) 
source('CDL_funcs.R')
source('species_params.R')

# 2. User input -----

# input results directory
newname <- "APHA_2017_ALL"

# Pest Specific, Multiple Life Stage Phenology Model Parameters:
# model scope
yr           <- 2017
start_doy    <- 1
end_doy      <- 365
region_param <- "ALL"
species      <- "APHA" # GCA/APHA/DCA
biotype      <- "S" # TODO: add options for each species

# introducing individual variation, tracked with simulation for each substage
# assign to 1 to match previous model versions
nsim <- 7 # number of substages/cohorts to approximate emergence distribution

# photoperiod decision inclusion
# 2 for logistic, 1 for single value CDL, 0 for none
model_CDL  <- 2 

# Derived parameters -----
REGION <- assign_extent(region_param = region_param)
params <- species_params(species, biotype, nsim, model_CDL)

# states <- map_data("state", xlim = c(REGION@xmin, REGION@xmax),
#                    ylim = c(REGION@ymin, REGION@ymax), lforce = "e")
# names(states)[1:2] <- c("x", "y")
# 
# us <- getData("GADM",country="USA",level=1)
# canada <- getData("GADM",country="CAN",level=1)
# 
# states <- crop(us, REGION)
# provs <- crop(canada, REGION)

if (!file.exists("./src/ref/ne_50m_admin_1_states_provinces_lakes/ne_50m_admin_1_states_provinces_lakes.dbf")){
  download.file(file.path('http://www.naturalearthdata.com/http/',
                          'www.naturalearthdata.com/download/50m/cultural',
                          'ne_50m_admin_1_states_provinces_lakes.zip'), 
                f <- tempfile())
  unzip(f, exdir = "./src/ref/ne_50m_admin_1_states_provinces_lakes")
  rm(f)
}

region <- readOGR("./src/ref/ne_50m_admin_1_states_provinces_lakes", 'ne_50m_admin_1_states_provinces_lakes', encoding='UTF-8')
reg.points = fortify(region, region="name_en")
reg.df = left_join(reg.points, region@data, by = c("id" = "name_en"))


# # is GDD already calculated if lifestages have same parameters?
# GDD <- brick(gdd_file)
# template <- GDD[[1]]
# template[!is.na(template)] <- 0
# template <- crop(template, REGION)

# Sites -----
# coordinates as examples
# # Galerucella
# sites <- data.frame(ID = c("Corvallis, OR", "Richland, WA", "JB Lewis-McChord, WA", "Palermo, CA",
#                            "Ephrata, WA", "Yakima Training Center, WA", "Camp Rilea, OR",
#                            "Ft Drum, NY", "West Point, NY", "Kellogg LTER, MI",
#                            "The Wilds, OH", "Duluth, MN", "Coeburn, VA", "Mountain Home AFB, ID",
#                            "Quantico MCB, VA", "Hanscom AFB, MA", "Ft Bragg, NC",
#                            "Ogden, UT", "Buckley AFB, CO", "S Portland, OR",
#                            "Sutherlin, OR", "Bellingham, WA"),
#                     x = c(-123.263, -119.283, -122.53, -121.625360, -119.555424, -120.461073,
#                           -123.934759, -75.763566, -73.962210, -85.402260, -81.733314,
#                           -92.158597, -82.466417, -115.865101, -77.311254, -71.276231,
#                           -79.083248, -112.052908, -104.752266, -122.658887,
#                           -123.315854, -122.479482),
#                     y = c(44.564, 46.275, 47.112, 39.426829, 47.318546, 46.680138,
#                           46.122867, 44.055684, 41.388456, 42.404749, 39.829447,
#                           46.728247, 36.943103, 43.044083, 38.513995, 42.457068,
#                           35.173401, 41.252509, 39.704018, 45.470532,
#                           43.387721, 48.756105))

# Diorhabda
sites <- data.frame(ID = c("TopockMarsh", "Lovelock", "GoldButte", "Delta", "BigBend",
                           "Bighorn", "Ft Carson", "Pinon canyon", "Yuma PG", "Nyssa"),
                    x = c(-114.5387, -118.5950, -114.2188, -112.9576, -114.6479,
                          -108.249, -104.756, -103.822, -114.265, -116.984),
                    y = c(34.7649, 40.04388, 36.73357, 39.14386, 35.10547,
                          44.998, 38.755, 37.417, 32.6837, 43.876))
sites$ID <- factor(sites$ID, c("Nyssa", "Bighorn", "Lovelock", "Delta", "Ft Carson", "Pinon canyon",
                               "GoldButte", "BigBend", "TopockMarsh", "Yuma PG"))


# Define grid order of sites -----
# mygrid <- data.frame(
#   code = c("JBLM", "DRUM", "EPHR", "BOIS", "DULU", "MASS", "WSPT", "YAKI", "RILE",
#            "UTAH", "KBST", "COEB", "WILD", "CORV", "DENV", "QUAN", "PLMO", "BRAG"),
#   name = c("JB Lewis-McChord, WA", "Ft Drum, NY", "Ephrata, WA", "Mountain Home AFB, ID",
#            "Duluth, MN", "Hanscom AFB, MA", "West Point, NY", "Yakima Training Center, WA",
#            "Camp Rilea, OR", "Ogden, UT", "Kellogg LTER, MI", "Coeburn, VA", "The Wilds, OH",
#            "Corvallis, OR", "Buckley AFB, CO", "Quantico MCB, VA", "Palermo, CA", "Ft Bragg, NC"),
#   row = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3),
#   col = c(1, 5, 2, 3, 4, 6, 6, 2, 1, 3, 4, 5, 4, 1, 3, 6, 2, 5),
#   stringsAsFactors = FALSE
# )
# 
# mygrid <- data.frame(
#   code = c("JBLM", "EPHR", "YAKI", "RILE",
#            "CORV", "PLMO", "PORT", "BELL", "SUTH"),
#   name = c("JB Lewis-McChord, WA", "Ephrata, WA",
#            "Yakima Training Center, WA",
#            "Camp Rilea, OR", 
#            "Corvallis, OR", "Palermo, CA", "S Portland, OR", 
#            "Bellingham, WA", "Sutherlin, OR"),
#   row = c(2, 1, 1, 2, 3, 3, 2, 1, 3),
#   col = c(3, 2, 3, 1, 1, 3, 2, 1, 2),
#   stringsAsFactors = FALSE
# )
# # geofacet::grid_preview(mygrid)



# 3. Quick default plots ------

res <- brick(paste(newname, "/", "NumGen2_001_weighted.grd", sep = ""))
# NAvalue(res) <- 200


plot(res[[365]])
plot(res[[seq(95, 130, 5)]])
plot(res[[seq(100, 360, 50)]])
plot(res[[seq(200, 300, 20)]])
plot(res[[seq(180, 280, 20)]])


# 4. Diapause/voltinism time series -----
# shows rapid change, only if within sensitive stage
# plot CDL as it moves through time and space

res <- brick(gdd_file)
res <- brick(paste(newname, "/", "LS4_001_weighted.grd", sep = ""))

test <- raster::extract(res, y = sites[, 2:3])
gdd <- raster::extract(GDD, y = sites[, 2:3])
plot(cumsum(gdd), test)


# plot weighted lifestages on same scale
res <- brick(paste(newname, "/", "diap_all.grd", sep = ""))
template <- res[[1]]
template[!is.na(template)] <- 0


names(res) <- paste("d", formatC(1:nlayers(res), width = 3, format = "d", flag = "0"), sep = "")
df = as.data.frame(res, xy=TRUE)
df1 <- df %>% 
  tidyr::gather(key = "DOY", value = "Perc_diapause", d001:d365) %>% 
  dplyr::mutate(DOY = as.numeric(gsub(pattern = "d", replacement = "", x = .$DOY))) %>% 
  dplyr::filter(DOY %in% seq(100, 365, 5))
df2 <- df1 %>% 
  filter(DOY %in% seq(190, 270, 10))


CDL <- (logit(.5) - params$CDL[2])/params$CDL[3]
FindApproxLat4CDL <- function(ras, doy, cdl, degtwil){
  ys <- seq(extent(ras)@ymin, extent(ras)@ymax, .01)
  hours <- photoperiod(lat = ys, doy = doy, p = degtwil)
  lat <- ys[which(abs(hours - cdl) == min(abs(hours - cdl)))]
  return(lat)
}
# Make this function select doy where photo time series crossing zero multiple times
# FindApproxDOY4CDL <- function(ras, lat, cdl, degtwil){
#   doy <- seq(1, 365, 1)
#   hours <- photoperiod(lat = lat, doy = doy, p = degtwil)
#   doys <- doy[which(abs(hours - cdl) == min(abs(hours - cdl)))]
#   return(doys)
# }


photo_df <- data.frame(DOY = seq(100, 365, 5))
photo_df <- photo_df %>% 
  rowwise() %>% 
  mutate(lat = FindApproxLat4CDL(template, DOY, CDL, 6)) %>% 
  filter(DOY %in% seq(190, 270, 10))


plt <- ggplot(df2, aes(x, y)) +
  geom_raster(aes(fill = Perc_diapause)) +
  scale_fill_viridis(na.value = "white") + 
  geom_polygon(data = states, aes(group = group), fill = NA, color = "black", size = .1) +
  geom_hline(data = photo_df, aes(yintercept = lat), color = "white") +
  facet_wrap(~DOY, nrow = 3) +
  coord_fixed(1.3) +
  # annotate("text", x = -73, y = 30, label = "Earliest substage") +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
plt

ggsave(paste("NEW","CDL", ".png", sep = ""),
       plot = plt, device = "png", width = 10, height = 6, units = "in")



# 5. Lifestage time series by biotype CDL --------
f <-list.files(path = newname, recursive = TRUE, full.names = TRUE)
rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]
rasfiles <- rasfiles[grep(pattern = "weighted", x = rasfiles, fixed = TRUE)]
# rasfiles <- rasfiles[-grep(pattern = "LS4", x = rasfiles, fixed = TRUE)]

tmp <- list()
for (i in rasfiles){
  fil <- stringr::str_split_fixed(i, pattern = "/", 3)
  ls <- fil[grep(pattern = "weighted", x = fil, fixed = TRUE)]
  ls <- unique(stringr::str_split_fixed(ls, pattern = "_", 2)[,1])
  photo_resp <- fil[-grep(pattern = "weighted", x = fil, fixed = TRUE)]
  photo_resp <- photo_resp[-grep(pattern = newname, x = photo_resp, fixed = TRUE)]
  
  
  res <- brick(i)
  names(res) <- paste("d", formatC(1:nlayers(res), width = 3, format = "d", flag = "0"), sep = "")
  
  stage_prop <- raster::extract(res, y = sites[, 2:3])
  stage_prop <- cbind(sites, stage_prop)
  stage_prop <- stage_prop %>% 
    tidyr::gather(key = "DOY", value = "Proportion", d001:d365) %>% 
    dplyr::mutate(DOY = as.numeric(gsub(pattern = "d", replacement = "", x = .$DOY)))
  gdd <- raster::extract(GDD, y = sites[, 2:3])
  gdd <- cbind(sites, gdd)
  gdd <- gdd %>% 
    tidyr::gather(key = "DOY", value = "GDD", layer.1:layer.365) %>% 
    dplyr::mutate(DOY = as.numeric(gsub(pattern = "layer.", replacement = "", x = .$DOY))) %>% 
    group_by(ID) %>% 
    arrange(DOY) %>% 
    dplyr::mutate(Accum_GDD = cumsum(GDD))
  outdf <- left_join(stage_prop, gdd)
  outdf$Lifestage <- ls
  outdf$photo_resp <- photo_resp
  tmp[[length(tmp) + 1]] <- outdf
}
tsdat <- bind_rows(tmp)


tsdat$Lifestage <- factor(tsdat$Lifestage, c("LSOW", "LS0", "LS1", "LS2", "LS3", "LS4", "NumGen"))
levels(tsdat$Lifestage) <- c("Overwinter", "Egg", "Larva", "Pupa", "Adult", "Diapause", "Voltinism")
# tsdat$ID <- factor(tsdat$ID, c("JB Lewis-McCord", "Richland", "Corvallis", "Yuba City"))

pltdat <- tsdat # %>% 
#   filter(Lifestage %in% c("Egg", "Diapause"))

# plt <- ggplot(pltdat, aes(x = DOY, y = Proportion, group = Lifestage, color = Lifestage)) +
#   geom_line(size = 2) +
#   facet_wrap(~ID, ncol = 1)
# plt
# ggsave(paste("NEW","LifestageTS", ".png", sep = ""),
#        plot = plt, device = "png", width = 8, height = 6, units = "in")

# multiply other lifestages with diapause proportion
diap <- pltdat %>% 
  dplyr::filter(Lifestage == "Diapause") %>% 
  dplyr::select(ID, DOY, Proportion, photo_resp) %>% 
  rename(Diap_Prop = Proportion)
pltdat1 <- pltdat %>% 
  dplyr::select(-photo_resp) %>%
  dplyr::filter(Lifestage %!in% c("Voltinism", "Diapause")) %>% 
  left_join(diap) %>% 
  # group_by(Lifestage) %>% 
  mutate(Proportion = Proportion * (1 - Diap_Prop))
# not sure how to do this: get voltinism from adult
pltdat2 <- tsdat %>% filter(Lifestage == "Diapause")
# # figuring out Voltinism from results
# pltdat3 <- tsdat %>%
#   filter(Lifestage == "Voltinism") %>% 
#   mutate(diapperc = diap) %>% 
#   group_by(ID) %>% 
#   mutate(Proportion = ifelse(diapperc < .8, Proportion, 0),
#          Proportion = cummax(Proportion))

pltdat <- bind_rows(pltdat1, pltdat2) %>% 
  filter(Lifestage %in% c("Egg", "Diapause")) %>%
  filter(ID %in% unique(pltdat1$ID)[c(1, 3, 4, 6, 7, 20:22)]) %>% 
  mutate(Date = as.Date(DOY, origin=as.Date("2015-12-31")),
         photo_resp = paste("Photo", photo_resp, sep = "_")) #%>% 
# filter(photo_resp %in% c("Photo_North_sol", "Photo_South_sol"))
# pltdat$photo_resp <- factor(pltdat$photo_resp, 
#                             c("Photo_Lovelock", "Photo_Delta", "Photo_GoldButte",
#                               "Photo_BigBend", "Photo_TopockMarsh"))

pltdat$ID <- droplevels(pltdat$ID)
pltdat$ID <- factor(pltdat$ID, levels(pltdat$ID)[c(1, 4, 8, 2, 6, 3, 7, 5)])

# pltdat <- tsdat %>% 
#   filter(Lifestage %in% c("Pupa", "Diapause"))
theme_set(theme_bw(base_size = 10)) 

plt <- ggplot(pltdat, aes(x = Accum_GDD, y = Proportion, group = Lifestage, color = Lifestage)) +
  geom_line(size = 2) +
  # scale_x_date(date_breaks = "2 month", date_labels = "%b") +
  # geom_rect(data = subset(pltdat, site == 'Big Bend State Park'), 
  #           aes(fill = site), xmin = -Inf, xmax = Inf,
  #           ymin = -Inf, ymax = Inf, alpha = 0.2)
  # facet_geo(~ID, grid = mygrid) +
  # geom_vline(xintercept = 100) +
  # coord_cartesian(xlim = c(75, 300)) 
  ggtitle(newname) +
  facet_grid(photo_resp~ID)
plt

ggsave(paste(newname,"LifestageTSgdd", ".png", sep = ""),
       plot = plt, device = "png", width = 12, height = 4, units = "in")


# plot of daily GDD by site
plt <- ggplot(gdd, aes(x = DOY, y = GDD, group = ID)) +
  geom_line(size = 2) +
  # facet_geo(~ID, grid = mygrid, scales = "free_y")
  # coord_cartesian(xlim = c(75, 300)) 
  facet_wrap(~ID, ncol = 1)
plt

# example plot of lifestages

pltdat <- bind_rows(pltdat1, pltdat2) %>% 
  filter(Lifestage %in% c("Egg", "Diapause")) %>%
  mutate(Date = as.Date(DOY, origin=as.Date("2015-12-31")),
         photo_resp = paste("Photo", photo_resp, sep = "_"))
pltdat$photo_resp <- factor(pltdat$photo_resp, 
                            c("Photo_Lovelock", "Photo_Delta", "Photo_GoldButte",
                              "Photo_BigBend", "Photo_TopockMarsh"))
exdat <- pltdat %>% 
  filter(ID == "S Portland, OR") %>% 
  filter(Lifestage != "Voltinism") %>% 
  mutate(Date = as.Date(DOY, origin=as.Date("2014-12-31")))
plt <- ggplot(exdat, aes(x = Date, y = Proportion, group = Lifestage, color = Lifestage)) +
  geom_line(size = 2) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b") +
  facet_grid(photo_resp~ID)
plt

ggsave(paste(newname,"Portland", ".png", sep = ""),
       plot = plt, device = "png", width = 8, height = 6, units = "in")



# 6. All lifestages by day ------
newname <- "new_west_SCDL_100twil"
# loop to grab all lifestage results on a certain day to plot together
# save in one giant data.frame to use for plots
returnwd <- getwd()
setwd(newname)

f <-list.files()
days <- seq(100, 350, 50)
rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]
rasfiles <- rasfiles[grep(pattern = "weighted", x = rasfiles, fixed = TRUE)]

dflist <- list()
for (i in rasfiles){
  res <- brick(i)
  # fix NA problem, just assigned as zero
  res <- res + template # template is all zeros and NA
  ls <- unique(stringr::str_split_fixed(i, pattern = "_", 2)[,1])
  for (j in days){
    df <- as.data.frame(res[[j]], xy=TRUE)
    names(df)[3] <- "Percent_of_simulations"
    df$doy <- j
    df$lifestage <- ls
    dflist[[length(dflist)+1]] <- df
  }
}
resdf <- dplyr::bind_rows(dflist)
names(resdf)[3] <- "Present"
setwd(returnwd)


pltlist <- list()
lifestages <- sort(unique(resdf$lifestage))[-6] # remove OW for even 6
ls_labels <- c("egg", "larva", "pupa", "adult", "diapause", "voltinism")
for (d in days){
  for (p in 1:length(lifestages)){
    pltdf <- resdf %>% 
      filter(lifestage == lifestages[p], 
             doy == d)
    tmpplt <- ggplot(pltdf, aes(x, y, fill = Present)) +
      geom_raster() +
      geom_polygon(data = states, aes(group = group), fill = NA, color = "black", size = .1) +
      theme_bw() +
      ggtitle(ls_labels[p])
    if(max(pltdf$Present, na.rm = TRUE) == 0){
      tmpplt <- tmpplt +
        scale_fill_viridis(na.value = "white", begin = 0, end = 0) 
    }else{
      tmpplt <- tmpplt +
        scale_fill_viridis(na.value = "white", begin = 0, end = 1)  
    }
    pltlist[[length(pltlist)+1]] <- tmpplt
  }
}

for (i in 1:length(days)){
  index <- (i*length(lifestages) - (length(lifestages) - 1)):(i*length(lifestages))
  plt <- grid.arrange(grobs = pltlist[index],
                      ncol = 3,
                      top = paste("DOY", days[i], sep = " "))
  ggsave(paste("NEWDOY", days [i], ".png", sep = ""),
         plot = plt, device = "png", width = 16, height = 9, units = "in")
}


# 7. Annual voltinism -------

setwd(newname)

# Diorhabda
sites <- sites[1:5, ]

days <- 365
f <-list.files(recursive = TRUE, full.names = TRUE)
rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]
rasfiles <- rasfiles[grep(pattern = "weighted", x = rasfiles, fixed = TRUE)]
photosite <- rasfiles[grep(pattern = "LS4", x = rasfiles, fixed = TRUE)]

for (ps in photosite){
  sitename <- stringr::str_split_fixed(string = ps, pattern = "/", n = 3)[2]
  volt <- brick(rasfiles[grep(pattern = "NumGen", x = rasfiles, fixed = TRUE)])
  template <- volt[[1]]
  template[!is.na(template)] <- 0
  
  volt <- volt + template
  diap <- brick(ps)
  diap <- diap + template
  threshold <- .75
  volt <- Cond(diap <= threshold, volt, 0)
  volt2 <- calc(volt, function(x){cummax(x)})
  
  df <- as.data.frame(volt2[[days]], xy=TRUE)
  names(df)[3] <- "Voltinism"
  
  
  # discrete voltinism classes for visual
  # df$Voltinism <- round(df$Voltinism/.5)*.5
  df$Voltinism <- round(df$Voltinism)
  maxvolt <- max(df$Voltinism, na.rm = TRUE)
  df$Voltinism <- as.factor(as.character(df$Voltinism))
  df <- df %>% 
    filter(!is.na(Voltinism))
  sites_plt <- sites
  sites_plt$show <- "A"
  sites_plt$show[sites_plt$ID == sitename] <- "B"
  
  tmpplt <- ggplot(data = df, aes(x, y, fill = Voltinism)) +
    geom_raster() +
    geom_polygon(data = states, aes(group = group), fill = NA, color = "black", inherit.aes = TRUE, size = .1) +
    theme_bw() +
    coord_fixed(1.3, xlim = c(-120.1925, -108.2658), ylim = c(31.51583, 42.3175)) +
    scale_fill_viridis(na.value = "white", begin = 0, end = maxvolt * .125, discrete = TRUE) + 
    geom_point(data = sites_plt, aes(x = x, y = y, color = show), size = 2, inherit.aes = FALSE) +
    scale_color_manual(values=c("white", "red")) +
    ggtitle(paste(sitename, "photoperiod response", sep = " ")) +
    guides(color = FALSE)
  tmpplt
  ggsave(paste(sitename,"_NumGen", ".png", sep = ""),
         plot = tmpplt, device = "png", width = 8, height = 4.8, units = "in")
}


# with partial generations considered
pltdf <- resdf %>% 
  filter(lifestage == "NumGen") %>% 
  mutate(Voltinism = Present)
filtdf <- resdf %>% 
  filter(lifestage == "LS4")
pltdf$Voltinism[which(filtdf$Present >= .75)] <- 0
pltdf <- pltdf %>% 
  group_by(x, y) %>% 
  summarise(Voltinism = max(Voltinism))


tmpplt <- ggplot(data = pltdf, aes(x, y, fill = Voltinism)) +
  geom_raster() +
  geom_polygon(data = states, aes(group = group), fill = NA, color = "black", inherit.aes = TRUE, size = .1) +
  theme_bw() +
  scale_fill_viridis(na.value = "white", begin = 0, end = 1)  +
  geom_point(data = sites, aes(x = x, y = y), size = 3, color = "black", inherit.aes = FALSE)
# geom_text(data = sites, aes(x = x, y = y, label = ID), inherit.aes = FALSE) + 
# ggtitle(ls_labels[p])
tmpplt
ggsave(paste("NEWPartial","NumGen", ".png", sep = ""),
       plot = tmpplt, device = "png", width = 6, height = 6, units = "in")


# APHALARA
# Compare N and S CDL and timing of overwinter emergence

years <- 2016 # c(2014:2016)
sims <- "sim1" # paste0("sim", c(1, 4, 7))
cdl <- "SCDL" # c("NCDL", "SCDL")
emrg <- c("GDD", "PHO")
params <- expand.grid(years, sims, cdl, emrg)
names(params) <- c("year", "sim", "cdl", "emrg")

outlist <- list()
for (i in 1:nrow(params)){
  yr <- params$year[i]
  sim <- params$sim[i]
  cdl <- params$cdl[i]
  emrg <- params$emrg[i]
  
  newname <- paste0("APHA_CONUS_", cdl, "_", yr, "_", emrg)
  
  days <- 365
  f <-list.files(path = newname, recursive = TRUE, full.names = TRUE)
  rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]
  rasfiles <- rasfiles[grep(pattern = sim, x = rasfiles, fixed = TRUE)]
  photosite <- rasfiles[grep(pattern = "LS4", x = rasfiles, fixed = TRUE)]
  
  volt <- brick(rasfiles[grep(pattern = "NumGen", x = rasfiles, fixed = TRUE)])
  volt <- volt + template
  diap <- brick(photosite)
  diap <- diap + template
  threshold <- .75
  volt <- Cond(diap <= threshold, volt, 0)
  volt2 <- calc(volt, function(x){cummax(x)})
  
  df <- as.data.frame(volt2[[days]], xy=TRUE)
  names(df)[3] <- "Voltinism"
  
  # discrete voltinism classes for visual
  df$Voltinism <- round(df$Voltinism)
  df$Voltinism <- as.factor(as.character(df$Voltinism))
  df$year <- yr
  df$sim <- sim
  df$cdl <- cdl
  df$emrg <- emrg
  
  outlist[[i]] <- df
}
outdf <- bind_rows(outlist)
levels(outdf$emrg) <- c("Degree-days", "Photoperiod")
# levels(outdf$sim) <- c("Early", "Mid", "Late")
levels(outdf$cdl) <- c("Northern", "Southern")

voltlevs <- unique(outdf$Voltinism)
voltlevs <- sort(voltlevs[is.na(voltlevs) == FALSE])
voltcols <- viridis(n = length(voltlevs), begin = 0, end = 1)
names(voltcols) <- voltlevs

for (yr in years){
  df <- outdf %>% filter(year == yr)
  tmpplt <- ggplot(data = df, aes(x, y, fill = Voltinism)) +
    geom_raster() +
    geom_polygon(data = states, aes(group = group), fill = NA, color = "black", inherit.aes = TRUE, size = .1) +
    theme_bw() +
    coord_fixed(1.3) +
    scale_fill_manual(values = voltcols) +
    # scale_fill_viridis(na.value = "white", begin = 0, end = numcol * .1, discrete = TRUE) + 
    ggtitle(paste0("Aphalara photoperiod response in ", yr)) +
    guides(color = FALSE) +
    facet_wrap( ~ emrg, ncol = 1)
  tmpplt
  ggsave(paste0("APHA_", yr,"_NumGen_PHO", ".png", sep = ""),
         plot = tmpplt, device = "png", width = 10, height = 8, units = "in")
  
}


# plot diapause over time
newnames <- c("APHA_CONUS_SCDL_2016_PHO", "APHA_CONUS_SCDL_2016_GDD", "APHA_CONUS_NCDL_2016_GDD")
outlist <- list()
for (nn in newnames){
  days <- 365
  f <-list.files(path = nn, recursive = TRUE, full.names = TRUE)
  rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]
  rasfiles <- rasfiles[grep(pattern = sim, x = rasfiles, fixed = TRUE)]
  
  eggs <- brick(rasfiles[grep(pattern = "LS0", x = rasfiles, fixed = TRUE)])
  eggs <- eggs + template
  eggs2 <- calc(eggs, function(x){min(which(x == 1))})
  
  diap <- brick(rasfiles[grep(pattern = "LS4", x = rasfiles, fixed = TRUE)])[[365]]
  diap <- diap + template
  
  df <- as.data.frame(eggs2, xy=TRUE)
  names(df)[3] <- "FirstEggs"
  
  df2 <- as.data.frame(diap, xy = TRUE)
  names(df2)[3] <- "Diap"
  
  outdf <- merge(df, df2)
  outdf$name <- nn
  outlist[[length(outlist)+1]] <- outdf
}
outdf <- bind_rows(outlist)
outdf$FirstEggs[which(is.infinite(outdf$FirstEggs))] <- NA

tmpplt <- ggplot(data = outdf, aes(x, y, fill = Diap)) +
  geom_raster() +
  geom_polygon(data = states, aes(group = group), fill = NA, color = "black", inherit.aes = TRUE, size = .1) +
  theme_bw() +
  coord_fixed(1.3) +
  # scale_fill_manual(values = voltcols) +
  scale_fill_viridis(na.value = "white", begin = 0, end = 1, discrete = FALSE) +
  ggtitle("Diapause at end of year") +
  guides(color = FALSE) +
  facet_wrap( ~ name, ncol = 1)
tmpplt
ggsave(paste0("APHA_Eggs_GDD", ".png", sep = ""),
       plot = tmpplt, device = "png", width = 10, height = 8, units = "in")




# 8. Simple plot of voltinism and diapause -----
doy <- 365
newname <- "APHA_2017_ALL"
f <- "NumGen1_all.grd"
res <- brick(paste(newname, "/", f, sep = ""))
df <- as.data.frame(res[[doy]], xy=TRUE)
df$var <- stringr::str_split_fixed(string = f, pattern = "_", n = 3)[1]
res <- brick(paste(newname, "/", "diap_all.grd", sep = ""))
df1 <- as.data.frame(res[[doy]], xy=TRUE)
df1$var <- "diapause"

outdf <- bind_rows(df, df1) %>% 
  mutate(proportion = layer.365 / 1000)

reg.df <- reg.df %>% 
  filter(long >= min(outdf$x), long <= max(outdf$x),
         lat >= min(outdf$y), lat <= max(outdf$y))

tmpplt <- ggplot(data = outdf, aes(x, y, fill = proportion)) +
  geom_raster() +
  geom_polygon(data = reg.df, aes(x = long, y = lat, group = group), fill = NA, color = "black", inherit.aes = FALSE, size = .3) +
  theme_bw() +
  coord_fixed(1.3) +
  # scale_fill_manual(values = voltcols) +
  scale_fill_viridis(na.value = "white", begin = 0, end = 1, discrete = FALSE) +
  ggtitle(paste0("Diapause and generation at DOY ", doy), subtitle = newname) +
  guides(color = FALSE) +
  facet_wrap( ~ var, ncol = 1) +
  scale_y_continuous(expand = c(0, 0) ) +
  scale_x_continuous(expand = c(0, 0) )
tmpplt

ggsave(paste0(newname, df$var[1], ".png"),
       plot = tmpplt, device = "png", width = 10, height = 12, units = "in")


# 9. More voltinism ----


# Need a voltinism map that averages the different generations

newname <- "APHA_2017_ALL"

returnwd <- getwd()
setwd(newname)
f <-list.files()
rasfiles <- f[grep(pattern = "NumGen", x = f, fixed = TRUE)]
rasfiles <- rasfiles[grep(pattern = "_all.grd", x = rasfiles, fixed = TRUE)]

ls_index <- stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,1]
ls_index <- as.numeric(gsub(pattern = "NumGen", replacement = "", fixed = TRUE, x = ls_index))

# most voltinism will be zero at a particular place, few with more than one value need to be averaged

# from template in model_lifecycle.R
# aggregates the numgen raster for each generation into one voltinism map at end of year
numgenstack <- mosaic(template, stack(rasfiles[1])[[365]], fun = max, na.rm = TRUE) * ls_index[1]
for (m in 2:length(rasfiles)){
  numgenstack <- addLayer(numgenstack, ls_index[m] * mosaic(template, stack(rasfiles[m])[[365]], fun = max, na.rm = TRUE))
  # numgenstack[numgenstack == 0] <- NA
  
}

numgen <- sum(numgenstack, na.rm = TRUE) / 1000 + template
plot(numgen)
writeRaster(numgen, filename = "numgen_365_gdd", datatype = "FLT4S")

# make numgen rasterbrick that aggregates the single generation bricks
# why did I remove this from model_lifecycle processing?



newname <- "APHA_2017_ALL"

returnwd <- getwd()
setwd(newname)
f <-list.files()
rasfiles <- f[grep(pattern = "NumGen", x = f, fixed = TRUE)]
rasfiles <- rasfiles[grep(pattern = "_all.grd", x = rasfiles, fixed = TRUE)]

ls_index <- stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,1]
ls_index <- as.numeric(gsub(pattern = "NumGen", replacement = "", fixed = TRUE, x = ls_index))

for (day in 1:365){
  
  numgenstack <- mosaic(template, stack(rasfiles[1])[[day]], fun = max, na.rm = TRUE) * ls_index[1]
  for (m in 2:length(rasfiles)){
    # numgenstack <- addLayer(numgenstack, mosaic(template, stack(rasfiles[m])[[day]], fun = max, na.rm = TRUE))
    numgenstack <- addLayer(numgenstack, ls_index[m] * mosaic(template, stack(rasfiles[m])[[day]], fun = max, na.rm = TRUE))
    
    }
  maxnumgen <- sum(numgenstack, na.rm = TRUE) / 1000 + template
  
  # maxnumgen <- which.max(numgenstack)
  # maxnumgen  <- setValues(maxnumgen, ls_index[getValues(maxnumgen)]) # always number files with leading zeros!
  
  if (!exists("numgen_stack")){
    numgen_stack <- stack(maxnumgen)
  } else {
    numgen_stack <- addLayer(numgen_stack, maxnumgen)
  }
}


plot(numgen_stack[[350]])

writeRaster(numgen_stack, filename = "NumGen", overwrite = TRUE)


numgen_stack <- brick("NumGen.grd")
numgen_stack <- numgen_stack[[1:365]]







# what about frost limiting?
# tmin from model_lifecycle.R
tmin <- crop(aggregate(tminfile[[180:365]], fact = 2, fun = mean, na.rm = TRUE, expand = TRUE), template)

test <- Cond(tmin <= -2, 1, 0)

first <- calc(test, fun=cumsum)
ffrost <- 179 + which.min(Cond(first == 0, 2, first)) # what I want is cumsum == 1 for first day
ffrost <- Cond(ffrost == 180, 365, ffrost) # set no frost places to 365
plot(ffrost)

numgen_frost <- stackSelect(numgen_stack, ffrost)

# numgen with diapause
volt <- numgen_stack
diap <- stack("diap_all.grd")
threshold <- .75
volt <- Cond(diap <= threshold * 1000, volt, 0)
volt2 <- calc(volt, function(x){cummax(x)})

# mismatch
ng_gdd <- brick("numgen_365_gdd.grd")
mmvolt <- volt2[[365]] - ng_gdd
plot(mmvolt)

# diverging color option
library(scales) # for muted
d + scale_colour_gradient2(low="red", high="blue")
d + scale_colour_gradient2(low=muted("red"), high=muted("blue"))



plttitle <- "Aphalara potential voltinism (no photoperiod) constrained by first hard frost (-2C)"
pltname <- "APHA_numgen_frost.png"
df <- as.data.frame(numgen_frost, xy=TRUE)
names(df)[3] <- "Voltinism"


# discrete voltinism classes for visual
# df$Voltinism <- round(df$Voltinism/.5)*.5
df$Voltinism <- round(df$Voltinism)
maxvolt <- max(df$Voltinism, na.rm = TRUE)
df$Voltinism <- factor(df$Voltinism, levels = c(1:maxvolt))
levels(df$Voltinism) <- c(levels(df$Voltinism)[1:6], rep("7+", 5))

df <- df %>% 
  filter(!is.na(Voltinism))


reg.df <- reg.df %>% 
  filter(long >= min(df$x), long <= max(df$x),
         lat >= min(df$y), lat <= max(df$y))
theme_set(theme_void(base_size = 20)) 

tmpplt <- ggplot(data = df, aes(x, y, fill = Voltinism)) +
  geom_raster() +
  geom_polygon(data = reg.df, aes(x = long, y = lat, group = group), fill = NA, color = "dark gray", inherit.aes = FALSE, size = .25) +
  coord_fixed(1.3) +
  scale_fill_viridis(name = "Generations", na.value = "white", begin = 0, end = maxvolt * .09, discrete = TRUE) +
  ggtitle(plttitle) +
  theme(legend.position = c(0.85, 0.4), plot.title = element_text(hjust = 0.5))
tmpplt
ggsave(pltname,
       plot = tmpplt, device = "png", width = 18, height = 11, units = "in")

