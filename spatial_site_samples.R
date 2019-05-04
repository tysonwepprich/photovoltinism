# temperature data for NA loosestrife sites
library(dplyr)
library(ggplot2)
library(raster)
library(lubridate)
library(stringr)

occur <- read.csv("data/edd_loose.csv", header = TRUE, skip = 3)
occur2 <- read.delim("data/gbif_loose.csv", header = TRUE, sep = "\t")

df <- occur %>% 
  filter(OccStatus != "Negative") %>% 
  dplyr::select(Latitude, Longitude, ObsDate) %>% 
  filter(complete.cases(.)) %>% 
  mutate(Date = dmy(ObsDate),
         Source = "EDDMAPS") %>% 
  dplyr::select(-ObsDate)

df2 <- occur2 %>% 
  mutate(Latitude = decimalLatitude,
         Longitude = decimalLongitude,
         Date = ymd(stringr::str_split_fixed(eventDate, fixed("T"), n = 2)[,1])) %>% 
  dplyr::select(Latitude, Longitude, Date) %>% 
  mutate(Source = "GBIF")

dat <- bind_rows(df, df2) %>% 
  filter(Longitude < -10, Latitude > 10) %>% 
  mutate(Year = year(Date))

dat <- dat %>% 
  dplyr::select(Latitude, Longitude) %>% 
  distinct()




region <- rgdal::readOGR("./src/ref/ne_50m_admin_1_states_provinces_lakes", 'ne_50m_admin_1_states_provinces_lakes', encoding='UTF-8')
# region <-spTransform(region, CRS(proj4string(template)))
reg.points = fortify(region, region="name_en")
reg.df = left_join(reg.points, region@data, by = c("id" = "name_en"))
theme_set(theme_bw(base_size = 20))


tmpplt <- ggplot(data = datcut, aes(x = Longitude, y = Latitude)) +
  geom_point(alpha = .3, color = "red") +
  geom_polygon(data = reg.df, aes(x = long, y = lat, group = group), fill = NA, color = "black", inherit.aes = FALSE, size = 1, alpha = .3) +
  coord_fixed(1.3, xlim = c(-127, -50), ylim = c(25, 55), expand = FALSE, clip = "on")
tmpplt


ji <- function(xy, origin=c(0,0), cellsize=c(.75,.75)) {
  t(apply(xy, 1, function(z) cellsize/2+origin+cellsize*(floor((z - origin)/cellsize))))
}
JI <- ji(cbind(dat$Longitude, dat$Latitude))
dat$X <- JI[, 1]
dat$Y <- JI[, 2]
dat$Cell <- paste(dat$X, dat$Y)

counts <- by(dat, dat$Cell, function(d) c(d$X[1], d$Y[1], nrow(d)))
counts.m <- matrix(unlist(counts), nrow=3)
rownames(counts.m) <- c("X", "Y", "Count")

count.max <- max(counts.m["Count",])
colors = sapply(counts.m["Count",], function(n) hsv(sqrt(n/count.max), .7, .7, .5))
plot(counts.m["X",] + 1/2, counts.m["Y",] + 1/2, #cex=sqrt(counts.m["Count",]/100),
     pch = 19, col=colors,
     xlab="Longitude of cell center", ylab="Latitude of cell center",
     main="Event counts within one-degree grid cells")


datcut <- dat %>% 
  group_by(Cell) %>% 
  slice(base::which.min(pointDistance(p1 = cbind(Longitude, Latitude), p2 = cbind(X, Y), lonlat = TRUE, allpairs = FALSE)))

coordinates(datcut) <- c("Longitude", "Latitude")

dists <- pointDistance(p1 = datcut, lonlat = TRUE, allpairs = TRUE)

points_matrix <- rgeos::gWithinDistance(datcut, dist = .25, byid = TRUE)
points_matrix[lower.tri(points_matrix, diag=TRUE)] <- NA
points_matrix

colSums(points_matrix, na.rm=TRUE) == 0
v <- colSums(points_matrix, na.rm=TRUE) == 0
datcut <- as.data.frame(datcut[v, ])

outdat <- datcut %>% filter(Latitude < 55, Longitude < -59.5)
saveRDS(outdat, "lythrum.rds")
