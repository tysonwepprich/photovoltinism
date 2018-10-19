# process Daymet data to extend N and S boundary of climate data past CONUS
library(raster)
library(prism)
library(daymetr)
library(maps)
library(maptools)

weather_path <- "daymet"

can <- map("world", "canada", fill=TRUE, col="transparent",
           plot=FALSE)
range(can$x, na.rm=TRUE)
range(can$y, na.rm=TRUE)
can$names
IDs <- sapply(strsplit(can$names, ":"), function(x) x[1])
can_sp <- map2SpatialPolygons(can, IDs=IDs,
                              proj4string=CRS("+proj=longlat +datum=WGS84"))

class(tile_outlines)

test <- intersect(tile_outlines, can_sp)
test2 <- test[test$Ymax < 55, ]



mex <- map("world", "mexico", fill=TRUE, col="transparent",
           plot=FALSE)
range(mex$x, na.rm=TRUE)
range(mex$y, na.rm=TRUE)
mex$names
IDs <- sapply(strsplit(mex$names, ":"), function(x) x[1])
mex_sp <- map2SpatialPolygons(mex, IDs=IDs,
                              proj4string=CRS("+proj=longlat +datum=WGS84"))

class(tile_outlines)

testmex <- intersect(tile_outlines, mex_sp)
testmex2 <- testmex[testmex$YMin > 23.5, ]

# list all the tiles included in canada and mexico to use
dl_tiles <- c(unique(test2$TileID), unique(testmex2$TileID))


# download daymet netcdf tiles, reproject to PRISM, save as R Raster (.grd)

yrs <- c(2013:2017)
prism_example <- raster("daymet/PRISM_tmax_stable_4kmD1_20170101_bil/PRISM_tmax_stable_4kmD1_20170101_bil.bil")
library(foreach)
library(doParallel)
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

outfiles <- foreach(year = 1:length(yrs),
                    .packages= c("raster", "daymetr"),
                    .inorder = FALSE) %:% 
  foreach(tile = 1:length(dl_tiles),
          .packages= c("raster", "daymetr"),
          .inorder = FALSE)%dopar%{
            
            yr <- yrs[year]
            id <- dl_tiles[tile]
            tmin <- paste0(weather_path, "/tmin_", yr, "_", id, ".nc")
            tmax <- paste0(weather_path, "/tmax_", yr, "_", id, ".nc")
            
            if(!file.exists(tmin)){
              
              download_daymet_tiles(tiles = id,
                                    start = yr, end = yr, path = weather_path, param = "tmin",
                                    silent = FALSE, force = FALSE)
            }
            if(!file.exists(tmax)){
              download_daymet_tiles(tiles = id,
                                    start = yr, end = yr, path = weather_path, param = "tmax",
                                    silent = FALSE, force = FALSE)
            }
            
            
            # reproject to prism CRS
            prism_res <- 25.875 / 621
            newproj <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
            
            x <- brick(tmin, varname = "tmin", level = 1)
            projectRaster(x, res = prism_res, crs = newproj, method="bilinear", 
                          filename = paste0(weather_path, "/tmin_", yr, "_", id, ".grd"))
            
            x <- brick(tmax, varname = "tmax", level = 1)
            projectRaster(x, res = prism_res, crs = newproj, method="bilinear", 
                          filename = paste0(weather_path, "/tmax_", yr, "_", id, ".grd"))
            
            # remove daymet netcdf
            if(file.exists(paste0(weather_path, "/tmin_", yr, "_", id, ".grd"))){
              file.remove(tmin)
            }
            if(file.exists(paste0(weather_path, "/tmax_", yr, "_", id, ".grd"))){
              file.remove(tmax)
            }
          }

if(exists("cl")){
  stopCluster(cl)
}



# stitch together canada and mexico tiles by year and tmin/tmax
library(stringr)
library(purrr)
tiles <- list.files(weather_path, recursive = FALSE)
tiles <- tiles[grep(pattern = ".grd", x = tiles, fixed = TRUE)]
tile_params <- str_split_fixed(tiles, pattern = coll("."), n = 2)[, 1]
temp <- str_split_fixed(tile_params, pattern = coll("_"), n = 3)
temps <- temp[, 1]
tileyears <- temp[, 2]
tilenames <- temp[, 3]
country <- dplyr::case_when(tilenames > 11500 ~ "canada",
                            tilenames < 11500 ~ "mexico")

cl <- makePSOCKcluster(20)
registerDoParallel(cl)

outfiles <- foreach(year = 1:length(yrs),
                    .packages= c("raster", "daymetr"),
                    .inorder = FALSE) %:% 
  foreach(nation = 1:2,
          .packages= c("raster", "daymetr"),
          .inorder = FALSE)%dopar%{
            
  
            yr <- yrs[year]
            cntry <- unique(country)[nation]
            index_tmin <- tiles[which(temps == "tmin" &
                                        tileyears == yr &
                                        country == cntry)]
            index_tmax <- tiles[which(temps == "tmax" &
                                        tileyears == yr &
                                        country == cntry)]
            
            setwd(weather_path)
            tmin <- paste0("tmin_", yr, "_", cntry, ".grd")
            tmax <- paste0("tmax_", yr, "_", cntry, ".grd")
            
            rasterOptions(tolerance = 0.5)
            
            AutoShift <- function(ras){
              r <- shift(ras, x = -1*origin(ras)[1], y = -1*origin(ras)[2])
            }
            
            if(!file.exists(tmin)){
              rlist <- lapply(index_tmin, brick)
              rlist <- lapply(rlist, AutoShift)
              rlist$filename <- tmin
              rlist$overwrite <- TRUE
              m <- do.call(raster::merge, rlist)
            }
            if(!file.exists(tmax)){
              rlist <- lapply(index_tmax, brick)
              rlist <- lapply(rlist, AutoShift)
              rlist$filename <- tmax
              rlist$overwrite <- TRUE
              rlist$tolerance <- 0.1
              m <- do.call(raster::merge, rlist)
            }
            
            

            # # remove daymet tiles
            # # DIDNT WORK on .gri files associated with .grd
            # if(file.exists(tmin)){
            #     lapply(index_tmin, FUN = file.remove) # CAREFUL HERE!
            #   
            # }
            # if(file.exists(tmax)){
            #     lapply(index_tmax, FUN = file.remove) # CAREFUL HERE!
            # }
            # setwd("../")
          }

if(exists("cl")){
  stopCluster(cl)
}