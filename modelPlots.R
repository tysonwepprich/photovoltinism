# summarize raster results from DDRP models

library(sp)
library(rgdal)
library(raster)

rasterOptions(overwrite = FALSE, 
              chunksize = 1e+07,
              maxmemory = 1e+08)

# for parallel simulations with control over seed for reproducibility
library(doRNG)
library(foreach) # for parallelized loops
library(doMC) 

library(ggplot2)
library(viridis)
library(mapdata)
library(dplyr)
library(tidyr)
theme_set(theme_bw(base_size = 14)) 

library(gridExtra)
library(grid)
# results directory
newname <- "new_west_SCDL_25twil"

source('CDL_funcs.R')
region_param <- "WEST"
gdd_file <- "meanGDD_07_13.grd"

REGION <- switch(region_param,
                 "CONUS"        = extent(-125.0,-66.5,24.0,50.0),
                 "NORTHWEST"    = extent(-125.1,-103.8,40.6,49.2),
                 "OR"           = extent(-124.7294, -116.2949, 41.7150, 46.4612),
                 "TEST"         = extent(-124, -122.5, 44, 45),
                 "WEST"         = extent(-125.14, -109, 37, 49.1))

states <- map_data("state", xlim = c(REGION@xmin, REGION@xmax),
                   ylim = c(REGION@ymin, REGION@ymax), lforce = "e")
names(states)[1:2] <- c("x", "y")

GDD <- brick(gdd_file)

template <- GDD[[1]]
template[!is.na(template)] <- 0
template <- crop(template, REGION)

# coordinates as examples
sites <- data.frame(ID = c("Corvallis", "Richland", "JB Lewis-McCord", "Yuba City"),
                    x = c(-123.263, -119.283, -122.53, -121.615),
                    y = c(44.564, 46.275, 47.112, 39.14))

# Take empirical distribution and calculate substages
# OW oviposition distribution
eggdist <- dbeta(x = seq(0, 1, length.out = 1000), 
                 shape1 = 3.888677, shape2 = 2.174208)
inputdist <- data.frame(x = seq(59.6, 223.3677, length.out = 1000),
                        y = eggdist)
inputdist$CDF <- cumsum(inputdist$y) / sum(inputdist$y, na.rm = TRUE)

nsim <- 7
substages <- SubstageDistrib(dist = inputdist, numstage = nsim, perc = .99)


####################################
# Weighted results by substage sizes
# Not needed for older model with only one parameter per stage
returnwd <- getwd()
setwd(newname)

f <-list.files()
rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]

# for each sim with unique information to save
ls <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,1])
# maps <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,2])
# sims <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,3])
# sims <- gsub(pattern = ".grd", replacement = "", x = sims)
ncores <- length(ls)
registerDoMC(cores = ncores)
# this loop takes a lot of time

foreach(i = ls, .packages = "raster") %dopar% {
  fs <- sort(rasfiles[grep(pattern = i, x = rasfiles, fixed = TRUE)])
  ll <- replicate(nlayers(brick(fs[1])), template)
  blank <- brick(ll)
  for (j in 1:length(fs)){
    ras_weighted <- brick(fs[j]) * substages[j, 2] # weights for each substage size
    blank <- overlay(blank, ras_weighted, fun=function(x,y) x + y)
  }
  outras <- writeRaster(blank, filename = paste(i, "weighted", sep = "_"),
                        overwrite = TRUE)
}

setwd(returnwd)
#####################
# quick plots of results

res <- brick(paste(newname, "/", "LS2_001_sim1.grd", sep = ""))
res <- brick(paste(newname, "/", "NumGen_001_sim7.grd", sep = ""))
NAvalue(res) <- 200

res <- brick(paste(newname, "/", "LS4_weighted.grd", sep = ""))
res <- brick(paste(newname, "/", "NumGen_weighted.grd", sep = ""))

# res <- crop(res, extent(-125.1,-103.8,40.6,49.2)) #NORTHWEST
# res <- crop(res, extent(-124.7294, -116.2949, 41.7150, 46.4612)) #OREGON
plot(res[[seq(100, 360, 50)]])
plot(res[[seq(150, 250, 20)]])
plot(res[[seq(180, 220, 5)]])
################

# diapause/numgen time series
# shows rapid change, only if within sensitive stage
# plot CDL as it moves through time and space


res <- brick(paste(newname, "/", "LS4_weighted.grd", sep = ""))
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
# CDL <- (logit(.5) - -56.9745)/3.5101 #Northern
CDL <- (logit(.5) - -60.3523)/3.8888 #Southern
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
# plot time series of each stage
returnwd <- getwd()
setwd(newname)

f <-list.files()
rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]
rasfiles <- rasfiles[grep(pattern = "weighted", x = rasfiles, fixed = TRUE)]

tmp <- list()
for (i in rasfiles){
  ls <- unique(stringr::str_split_fixed(i, pattern = "_", 2)[,1])
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
  tmp[[length(tmp) + 1]] <- outdf
}
tsdat <- bind_rows(tmp) %>% 
  filter(Lifestage != "NumGen")

setwd(returnwd)


tsdat$Lifestage <- factor(tsdat$Lifestage, c("LSOW", "LS0", "LS1", "LS2", "LS3", "LS4"))
levels(tsdat$Lifestage) <- c("Overwinter", "Egg", "Larva", "Pupa", "Adult", "Diapause")
tsdat$ID <- factor(tsdat$ID, c("JB Lewis-McCord", "Richland", "Corvallis", "Yuba City"))

pltdat <- tsdat %>% 
  filter(Lifestage %in% c("Egg", "Adult", "Diapause"))
                          
# plt <- ggplot(pltdat, aes(x = DOY, y = Proportion, group = Lifestage, color = Lifestage)) +
#   geom_line(size = 2) +
#   facet_wrap(~ID, ncol = 1)
# plt
# ggsave(paste("NEW","LifestageTS", ".png", sep = ""),
#        plot = plt, device = "png", width = 8, height = 6, units = "in")

# multiply egg and diapause
diap <- pltdat$Proportion[pltdat$Lifestage == "Diapause"]
pltdat1 <- pltdat %>% 
  filter(Lifestage == "Adult") %>%
  # filter(Lifestage %in% c("Egg", "Adult")) %>% 
  mutate(Proportion = Proportion * (1 - diap))
pltdat2 <- pltdat %>% filter(Lifestage == "Diapause") 
pltdat <- rbind(pltdat1, pltdat2)

plt <- ggplot(pltdat, aes(x = Accum_GDD, y = Proportion, group = Lifestage, color = Lifestage)) +
  geom_line(size = 2) +
  facet_wrap(~ID, ncol = 1)
plt
ggsave(paste(newname,"LifestageTS", ".png", sep = ""),
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

pltdf <- resdf %>% 
  filter(lifestage == "NumGen", 
         doy == 350) %>% 
  mutate(Voltinism = Present)
tmpplt <- ggplot(data = pltdf, aes(x, y, fill = Voltinism)) +
  geom_raster() +
  geom_polygon(data = states, aes(group = group), fill = NA, color = "black", inherit.aes = TRUE, size = .1) +
  theme_bw() +
  scale_fill_viridis(na.value = "white", begin = 0, end = 1)  +
  geom_point(data = sites, aes(x = x, y = y), size = 3, color = "black", inherit.aes = FALSE)
  # geom_text(data = sites, aes(x = x, y = y, label = ID), inherit.aes = FALSE) + 
  # ggtitle(ls_labels[p])
tmpplt
ggsave(paste("NEW","NumGen", ".png", sep = ""),
       plot = tmpplt, device = "png", width = 6, height = 6, units = "in")

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



