# plot rasters using ggplots, and consistent scaling
library(raster)
library(ggplot2)
library(viridis)
library(ggthemes)

prism_path <- "/data/PRISM/2014"
source('CDL_funcs.R')


pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
files <- list.files(path = prism_path, pattern=pattern, all.files=FALSE, full.names=TRUE)
test <- raster(files[1]) # from FCM code
prismdate <- lubridate::ymd(stringr::str_split_fixed(files[1], pattern = "_", 6)[5])

test_spdf <- as(test, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")

testmap <- ggplot() +  
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8) + 
  # geom_polygon(data=OR, aes(x=long, y=lat, group=group), 
  #              fill=NA, color="grey50", size=0.25) +
  scale_fill_viridis(name = "PRISM TMIN(C)") +
  # scale_fill_viridis() +
  coord_equal() +
  theme_map() +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm")) +
  ggtitle(prismdate)
  # theme(legend.title = element_text(paste(prismdate, "PRISM TMIN(C)", sep = " ")))
ggsave("testmap.png", plot = testmap, device = "png", dpi = 300, 
       width = 5, height = 4, units = "in")
