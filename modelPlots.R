# summarize raster results from DDRP models

library(sp)
library(rgdal)
library(raster)

rasterOptions(overwrite = FALSE, 
              chunksize = 1e+07,
              maxmemory = 1e+08,
              tmpdir = "~/REPO/photovoltinism/rastertmp/")

# for parallel simulations with control over seed for reproducibility
library(doRNG)
library(foreach) # for parallelized loops
library(doMC) # for LINUX
library(doSNOW) # for WINDOWS


library(ggplot2)
library(viridis)
library(mapdata)
library(dplyr)
library(tidyr)
library(geofacet)
theme_set(theme_bw(base_size = 12)) 

library(gridExtra)
library(grid)
# results directory
# newname <- "GCA_WEST_SCDL_2017"
# newname <- "DCA_SW_2016"
newname <- "APHA_CONUS_NCDL_GDD"

source('CDL_funcs.R')
region_param <- "CONUS"
gdd_file <- "dailygdd.grd"

REGION <- assign_extent(region_param = region_param)

states <- map_data("state", xlim = c(REGION@xmin, REGION@xmax),
                   ylim = c(REGION@ymin, REGION@ymax), lforce = "e")
names(states)[1:2] <- c("x", "y")

GDD <- brick(gdd_file)

template <- GDD[[1]]
template[!is.na(template)] <- 0
template <- crop(template, REGION)

# coordinates as examples
# # Galerucella
sites <- data.frame(ID = c("Corvallis, OR", "Richland, WA", "JB Lewis-McChord, WA", "Palermo, CA",
                           "Ephrata, WA", "Yakima Training Center, WA", "Camp Rilea, OR",
                           "Ft Drum, NY", "West Point, NY", "Kellogg LTER, MI",
                           "The Wilds, OH", "Duluth, MN", "Coeburn, VA", "Mountain Home AFB, ID",
                           "Quantico MCB, VA", "Hanscom AFB, MA", "Ft Bragg, NC",
                           "Ogden, UT", "Buckley AFB, CO", "S Portland, OR",
                           "Sutherlin, OR", "Bellingham, WA"),
                    x = c(-123.263, -119.283, -122.53, -121.625360, -119.555424, -120.461073,
                          -123.934759, -75.763566, -73.962210, -85.402260, -81.733314,
                          -92.158597, -82.466417, -115.865101, -77.311254, -71.276231,
                          -79.083248, -112.052908, -104.752266, -122.658887,
                          -123.315854, -122.479482),
                    y = c(44.564, 46.275, 47.112, 39.426829, 47.318546, 46.680138,
                          46.122867, 44.055684, 41.388456, 42.404749, 39.829447,
                          46.728247, 36.943103, 43.044083, 38.513995, 42.457068,
                          35.173401, 41.252509, 39.704018, 45.470532,
                          43.387721, 48.756105))

# # Diorhabda
# sites <- data.frame(ID = c("TopockMarsh", "Lovelock", "GoldButte", "Delta", "BigBend"),
#                     x = c(-114.5387, -118.5950, -114.2188, -112.9576, -114.6479),
#                     y = c(34.7649, 40.04388, 36.73357, 39.14386, 35.10547))
# sites$ID <- factor(sites$ID, c("Lovelock", "Delta", "GoldButte", "BigBend", "TopockMarsh"))


nsim <- 7

# Try 2 with skewed t
arg1 = list(mu = 97.94, sigma2 = 2241.7, shape = 3.92, nu = 9.57)
x = seq(50, 350, length.out = 1000)
y = mixsmsn:::dt.ls(x, loc = arg1$mu, sigma2 = arg1$sigma2, shape = arg1$shape, nu = arg1$nu)
inputdist <- data.frame(x = x, y = y) %>% 
  arrange(x) %>% 
  mutate(CDF = cumsum(y/sum(y)))
substages <- SubstageDistrib(dist = inputdist, numstage = 7, perc = .99)
# To get observations to fit for overwinter adults and F1 eggs, 
# overwinter pre-oviposition period is only 50 deg days
substages$means <- substages$means + 15


# APHALARA SUBSTAGES

# skewed t ballpark estimation from Len's OV percentiles
arg1 = list(mu = 250, sigma2 = 25000, shape = 3, nu = 10)
x = seq(150, 950, length.out = 1000)
y = mixsmsn:::dt.ls(x, loc = arg1$mu, sigma2 = arg1$sigma2, shape = arg1$shape, nu = arg1$nu)
plot(x, y)
inputdist <- data.frame(x = x, y = y) %>% 
  arrange(x) %>% 
  mutate(CDF = cumsum(y/sum(y)))
substages <- SubstageDistrib(dist = inputdist, numstage = nsim, perc = .99)


newdirs <- c(
  # "GCA_NWSMALL_2014", "GCA_NWSMALL_2015", "GCA_NWSMALL_2016",
             # "GCA_NWSMALL_2014/North", "GCA_NWSMALL_2015/North", "GCA_NWSMALL_2016/North",
             # "GCA_NWSMALL_2014/South", "GCA_NWSMALL_2015/South", "GCA_NWSMALL_2016/South",
             # "GCA_NWSMALL_2014/North_sol", "GCA_NWSMALL_2015/North_sol", "GCA_NWSMALL_2016/North_sol",
             "GCA_NWSMALL_2014/South_sol", "GCA_NWSMALL_2015/South_sol", "GCA_NWSMALL_2016/South_sol")
for (newname in newdirs){
  ####################################
  # Weighted results by substage sizes
  # Not needed for older model with only one parameter per stage
  returnwd <- getwd()
  setwd(newname)
  
  f <-list.files()
  rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]
  rasfiles <- rasfiles[grep(pattern = "sim", x = rasfiles, fixed = TRUE)]
  
  # for each sim with unique information to save
  ls <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,1])
  maps <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,2])
  sims <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,3])
  sims <- gsub(pattern = ".grd", replacement = "", x = sims)
  
  # parallel backend for foreach loop
  if(.Platform$OS.type == "unix"){
    ncores <- length(ls) * length(maps) / 2
    registerDoMC(cores = ncores)
  }else if(.Platform$OS.type == "windows"){
    ncores <- parallel::detectCores() - 1
    cl <- makeCluster(ncores)
    registerDoSNOW(cl)
  }
  
  # run inside foreach now
  # tmppath <- paste0("~/REPO/photovoltinism/rastertmp/", "run", newname)
  # dir.create(path = tmppath, showWarnings = FALSE)
  # #sets temp directory
  # rasterOptions(tmpdir=file.path(tmppath)) 
  
  # this loop takes a lot of time
  
  system.time({
    foreach(i = ls,
            .packages= "raster",
            .export = c("rasfiles", "substages")) %:% 
      foreach(m = maps,
              .packages = "raster",
              .export = c("rasfiles", "substages")) %dopar%{
                # .inorder = FALSE) %do%{
                # foreach(i = ls, .packages = "raster") %dopar% {
                tmppath <- paste0("~/REPO/photovoltinism/rastertmp/", i, "run", m)
                dir.create(path = tmppath, showWarnings = FALSE)
                #sets temp directory
                rasterOptions(tmpdir=file.path(tmppath)) 
                
                fs <- sort(rasfiles[grep(pattern = i, x = rasfiles, fixed = TRUE)])
                fs <- sort(fs[grep(pattern = m, x = fs, fixed = TRUE)])
                template <- brick(fs[1])[[1]]
                template[!is.na(template)] <- 0
                
                ll <- replicate(nlayers(brick(fs[1])), template)
                blank <- brick(ll)
                for (j in 1:length(fs)){
                  ras_weighted <- brick(fs[j]) * substages[j, 2] # weights for each substage size
                  blank <- overlay(blank, ras_weighted, fun=function(x,y) x + y)
                }
                outras <- writeRaster(blank, filename = paste(i, m, "weighted", sep = "_"),
                                      overwrite = TRUE)
                removeTmpFiles(h = 0)
                unlink(tmppath, recursive = TRUE)
                
              }
  })
  
  stopCluster(cl) #WINDOWS
  setwd(returnwd)
  
}
# # # remove non-weighted results to save disk space
# # however, need results of each sim to model different photoperiod without total rerun
# cleanup <- list.files()
# cleanup <- cleanup[-grep("weighted", x = cleanup)]
# lapply(cleanup, FUN = file.remove) # CAREFUL HERE!

# mosaic maps together if CONUS
# maybe hold off on this since its files so large
# # LINUX
# ncores <- length(ls)
# registerDoMC(cores = ncores)
# WINDOWS
ncores <- 4
cl <- makeCluster(ncores)
registerDoSNOW(cl)
f <-list.files()
rasfiles <- f[grep(pattern = "weighted", x = f, fixed = TRUE)]
lstage <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,1])

system.time({
  foreach(i = lstage,
         .packages = "raster") %dopar%{
           fs <- rasfiles[grep(pattern = i, x = rasfiles, fixed = TRUE)]
           
           bricklist <- list()
           for (m in 1:length(fs)){
             bricklist[[m]] <- brick(fs[m])
           }
           
           bricklist$filename <- paste(i, "all", sep = "_")
           bricklist$overwrite <- TRUE
           test <- do.call(merge, bricklist)
         }
  
  cleanup <- list.files()
  cleanup <- cleanup[-grep("all", x = cleanup)]
  lapply(cleanup, FUN = file.remove) # CAREFUL HERE!
})
stopCluster(cl) #WINDOWS
setwd(returnwd)
#####################
# quick plots of results

res <- brick(paste(newname, "/", "LS_001_si.grd", sep = ""))
res <- brick(paste(newname, "/", "NumGen_001_sim7.grd", sep = ""))
NAvalue(res) <- 200

res <- brick(paste(newname, "/", "LS4_001_sim4.grd", sep = ""))
res <- brick(paste(newname, "/", "NumGen_002_weighted.grd", sep = ""))

# res <- crop(res, extent(-125.1,-103.8,40.6,49.2)) #NORTHWEST
# res <- crop(res, extent(-124.7294, -116.2949, 41.7150, 46.4612)) #OREGON
plot(res[[365]])
plot(res[[seq(100, 360, 50)]])
plot(res[[seq(200, 300, 20)]])
plot(res[[seq(180, 280, 20)]])
################

# diapause/numgen time series
# shows rapid change, only if within sensitive stage
# plot CDL as it moves through time and space

res <- brick(gdd_file)
res <- brick(paste(newname, "/", "LS4_001_weighted.grd", sep = ""))
# 
# test <- vector()
# gdd <- vector()
# cell <- 100000
# tmpGDD <- crop(GDD, template)
# 
# # first try
# for (i in 1:nlayers(res)){
#   test[i] <- res[[i]][cell]
#   gdd[i] <- tmpGDD[[i]][cell]
# }
# plot(cumsum(gdd), test)

# with extract, much easier
test <- raster::extract(res, y = sites[, 2:3])
gdd <- raster::extract(GDD, y = sites[, 2:3])
plot(cumsum(gdd), test)


# plot weighted lifestages on same scale
states <- map_data("state", xlim = c(REGION@xmin, REGION@xmax),
                   ylim = c(REGION@ymin, REGION@ymax), lforce = "e")
names(states)[1:2] <- c("x", "y")

names(res) <- paste("d", formatC(1:nlayers(res), width = 3, format = "d", flag = "0"), sep = "")
df = as.data.frame(res, xy=TRUE)
df1 <- df %>% 
  tidyr::gather(key = "DOY", value = "Perc_diapause", d001:d365) %>% 
  dplyr::mutate(DOY = as.numeric(gsub(pattern = "d", replacement = "", x = .$DOY))) %>% 
  dplyr::filter(DOY %in% seq(100, 365, 5))
df2 <- df1 %>% 
  filter(DOY %in% seq(120, 230, 15))


logit <- function(p){log(p/(1-p))}
CDL <- (logit(.5) - 90.092)/-6.115 #Northern
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
  filter(DOY %in% seq(120, 230, 15))


plt <- ggplot(df2, aes(x, y)) +
  geom_raster(aes(fill = Perc_diapause)) +
  scale_fill_viridis(na.value = "white") + 
  geom_polygon(data = states, aes(group = group), fill = NA, color = "black", size = .1) +
  geom_hline(data = photo_df, aes(yintercept = lat), color = "white") +
  facet_wrap(~DOY, nrow = 2) +
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



######
# Number of generations when diapause is a %???




######
mygrid <- data.frame(
  code = c("JBLM", "DRUM", "EPHR", "BOIS", "DULU", "MASS", "WSPT", "YAKI", "RILE",
           "UTAH", "KBST", "COEB", "WILD", "CORV", "DENV", "QUAN", "PLMO", "BRAG"),
  name = c("JB Lewis-McChord, WA", "Ft Drum, NY", "Ephrata, WA", "Mountain Home AFB, ID",
           "Duluth, MN", "Hanscom AFB, MA", "West Point, NY", "Yakima Training Center, WA",
           "Camp Rilea, OR", "Ogden, UT", "Kellogg LTER, MI", "Coeburn, VA", "The Wilds, OH",
           "Corvallis, OR", "Buckley AFB, CO", "Quantico MCB, VA", "Palermo, CA", "Ft Bragg, NC"),
  row = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3),
  col = c(1, 5, 2, 3, 4, 6, 6, 2, 1, 3, 4, 5, 4, 1, 3, 6, 2, 5),
  stringsAsFactors = FALSE
)

mygrid <- data.frame(
  code = c("JBLM", "EPHR", "YAKI", "RILE",
           "CORV", "PLMO", "PORT", "BELL", "SUTH"),
  name = c("JB Lewis-McChord, WA", "Ephrata, WA",
           "Yakima Training Center, WA",
           "Camp Rilea, OR", 
           "Corvallis, OR", "Palermo, CA", "S Portland, OR", 
           "Bellingham, WA", "Sutherlin, OR"),
  row = c(2, 1, 1, 2, 3, 3, 2, 1, 3),
  col = c(3, 2, 3, 1, 1, 3, 2, 1, 2),
  stringsAsFactors = FALSE
)
# geofacet::grid_preview(mygrid)

newname <- "GCA_NWSMALL_2015"
gdd_file <- paste(newname, "dailygdd.grd", sep = "/")
GDD <- brick(gdd_file)
######
# Plots of lifestage time series with different photoperiod responses
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


######
# Plot all lifestages for a particular day
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


# plot just NumGen

setwd(newname)

days <- 365
f <-list.files(path = newname, recursive = TRUE, full.names = TRUE)
rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]
rasfiles <- rasfiles[grep(pattern = "weighted", x = rasfiles, fixed = TRUE)]
photosite <- rasfiles[grep(pattern = "LS4", x = rasfiles, fixed = TRUE)]

for (ps in photosite){
  sitename <- stringr::str_split_fixed(string = ps, pattern = "/", n = 3)[2]
  volt <- brick(rasfiles[grep(pattern = "NumGen", x = rasfiles, fixed = TRUE)])
  volt <- volt + template
  diap <- brick(ps)
  diap <- diap + template
  threshold <- .75
  volt <- Cond(diap <= threshold, volt, 0)
  volt2 <- calc(volt, fun = function(x){cummax(x)})
  
  df <- as.data.frame(volt2[[days]], xy=TRUE)
  names(df)[3] <- "Voltinism"
  
  # discrete voltinism classes for visual
  # df$Voltinism <- round(df$Voltinism/.5)*.5
  df$Voltinism <- round(df$Voltinism)
  df$Voltinism <- as.factor(as.character(df$Voltinism))
  
  sites_plt <- sites
  sites_plt$show <- "A"
  sites_plt$show[sites_plt$ID == sitename] <- "B"
  
  tmpplt <- ggplot(data = df, aes(x, y, fill = Voltinism)) +
    geom_raster() +
    geom_polygon(data = states, aes(group = group), fill = NA, color = "black", inherit.aes = TRUE, size = .1) +
    theme_bw() +
    coord_fixed(1.3) +
    scale_fill_viridis(na.value = "white", begin = 0, end = 1, discrete = TRUE) + 
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



