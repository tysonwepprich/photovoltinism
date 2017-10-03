# summarize raster results from DDRP models


####################################
# Weighted results by substage sizes

returnwd <- getwd()
setwd(newname)

f <-list.files()
rasfiles <- f[grep(pattern = ".grd", x = f, fixed = TRUE)]

# for each sim with unique information to save
ls <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 2)[,1])
# maps <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,2])
# sims <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,3])
# sims <- gsub(pattern = ".grd", replacement = "", x = sims)

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


res <- brick(paste(newname, "/", "LS2_001_sim4.grd", sep = ""))
res <- brick(paste(newname, "/", "NumGen_001_sim7.grd", sep = ""))
NAvalue(res) <- 200

res <- brick(paste(newname, "/", "LS4_weighted.grd", sep = ""))
res <- brick(paste(newname, "/", "NumGen_weighted.grd", sep = ""))

# res <- crop(res, extent(-125.1,-103.8,40.6,49.2)) #NORTHWEST
# res <- crop(res, extent(-124.7294, -116.2949, 41.7150, 46.4612)) #OREGON
plot(res[[seq(100, 360, 50)]])
plot(res[[seq(150, 250, 20)]])
plot(res[[seq(180, 220, 5)]])


# diapause/numgen time series
# shows rapid change, only if within sensitive stage
test <- vector()
gdd <- vector()
cell <- 100000
tmpGDD <- crop(GDD, template)

for (i in 1:nlayers(res)){
  test[i] <- res[[i]][cell]
  gdd[i] <- tmpGDD[[i]][cell]
}
plot(cumsum(gdd), test)

# plot weighted lifestages on same scale
library(ggplot2)
library(viridis)
library(mapdata)
library(dplyr)
library(tidyr)
states <- map_data("state", xlim = c(REGION@xmin, REGION@xmax),
                   ylim = c(REGION@ymin, REGION@ymax), lforce = "e")
names(states)[1:2] <- c("x", "y")


df = as.data.frame(res, xy=TRUE)
df1 <- df %>% 
  gather(key = "DOY", value = "Perc_diapause", layer.1:layer.365) %>% 
  mutate(DOY = as.numeric(stringr::str_split_fixed(DOY, pattern = "layer.", n = 2)[, 2])) %>% 
  filter(DOY %in% seq(100, 365, 5))
df2 <- df1 %>% 
  filter(DOY %in% seq(150, 220, 5))


logit <- function(p){log(p/(1-p))}
# CDL <- (logit(.5) - -56.9745)/3.5101 #Northern
CDL <- (logit(.5) - -60.3523)/3.8888 #Southern
FindApproxLat4CDL <- function(ras, doy, cdl){
  ys <- seq(extent(ras)@ymin, extent(ras)@ymax, .01)
  hours <- photoperiod(lat = ys, doy = doy, p = 6)
  lat <- ys[which(abs(hours - cdl) == min(abs(hours - cdl)))]
  return(lat)
}


photo_df <- data.frame(DOY = seq(100, 365, 5))
photo_df <- photo_df %>% 
  rowwise() %>% 
  mutate(lat = FindApproxLat4CDL(template, DOY, CDL)) %>% 
  filter(DOY %in% seq(150, 220, 5))


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




######
# Number of generations when diapause is a %???








