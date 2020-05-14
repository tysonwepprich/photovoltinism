library(OpenStreetMap)
library(rgdal)
map <- openmap(c(70,-179), c(-70,179))
plot(map)

map <- openmap(c(52,-105),
               c(26,-125),
               type="stamen-terrain")
plot(map)


library(ggmap)
library(ggplot2)
library(ggsn)
library(ggrepel)
#  Loading required package: ggplot2
#  Registered S3 methods overwritten by 'ggplot2':
#    method         from 
#    [.quosures     rlang
#    c.quosures     rlang
#    print.quosures rlang
#  Google's Terms of Service: https://cloud.google.com/maps-platform/terms/.
#  Please cite ggmap if you use it! See citation("ggmap") for details.

us <- c(left = -125, bottom = 30, right = -105, top = 50)
west <- get_stamenmap(us, zoom = 5, maptype = "terrain-background") 

sites # from DCA_chamber_diapause2.R
sites2 <- data.frame(site = c(
                           "Ephrata, WA", "Yakima TC, WA",
                           "Sutherlin, OR", "Bellingham, WA",
                           "McArthur, CA", "Palermo, CA", "Rickreall, OR"),
                    lon = c(-119.655253, -119.986530,
                          -123.315854, -122.479482,
                          -121.41, -121.58, -123.27),
                    lat = c(47.161647, 46.756318,
                          43.387721, 48.756105,
                          41.10, 39.41, 44.98),
                    Species = "Galerucella calmariensis")
allsites <- bind_rows(sites, sites2) %>% data.frame()
allsites$long <- allsites$lon
# 
anc <- c(x = -122, y = 32)
# n_anc <- c(x = -116.8, y = 48)

plt <- ggmap(west) +
  geom_point(data = allsites, aes(x = long, y = lat, shape = Species), size = 3) +
  # scale_shape(name = NULL) +
  geom_text_repel(data = allsites, aes(x = long, y = lat, label = site), hjust = "left", nudge_x = 0.2, size = 5) +
  xlab("Longitude") + ylab("Latitude") +
  theme(legend.position = c(.75, .9),
        legend.background = element_rect(fill=alpha('white', 0.5)),
        legend.key = element_rect(fill = alpha("white", 0.0)),
        legend.text = element_text(face = "italic")) +
  ggsn::scalebar(data = allsites, transform = TRUE, dist = 200, st.size=5, 
                 height=0.025, dist_unit = "km", model = 'WGS84', location = "bottomleft")
plt



# what about point of invasive plants?
library(tidyverse)
dfs <- list.files("C:/Users/wepprict/Box Sync/Box Sync/SERDP Phenology Photoperiod/EDDMAPS", full.names = TRUE)
dat <- read.csv(dfs[5], header = TRUE, skip = 3, stringsAsFactors = FALSE)

ggplot(data = dat, aes(x = Longitude, y = Latitude, color = Status)) +
  geom_point()
  








rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}


# Diorhabda

us <- c(left = -118.5, bottom = 30, right = -111, top = 40)
west <- get_stamenmap(us, zoom = 5, maptype = "terrain-background") 

sites <- data.frame(ID = c(
  "Ephrata, WA", "Yakima TC, WA",
  "Sutherlin, OR", "Bellingham, WA",
  "McArthur, CA", "Palermo, CA", "Rickreall, OR"),
  long = c(-119.655253, -119.986530,
           -123.315854, -122.479482,
           -121.41, -121.58, -123.27),
  lat = c(47.161647, 46.756318,
          43.387721, 48.756105,
          41.10, 39.41, 44.98))

anc <- c(x = -116.8, y = 39.1)
n_anc <- c(x = -116.8, y = 48)

sites <- sites %>% filter(lat < 40) %>% dplyr::select(-year, -rowid) %>% distinct()


plt <- ggmap(west) +
  geom_point(data = sites, aes(x = lon, y = lat), size = 3) +
  geom_text(data = sites, aes(x = lon, y = lat, label = site), hjust = "left", nudge_x = 0.2, size = 5) +
  xlab("Longitude") + ylab("Latitude") +
  # coord_proj(
  #   paste0("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96",
  #          " +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")) +
  scale_y_continuous(limits = c(32, 40), expand = c(0, 0))
  # ggsn::scalebar(data = sites, anchor = anc,
  #                dist = 75, dist_unit = "km", st.size = 4,
  #                transform = TRUE, model = "WGS84")
plt


