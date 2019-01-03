# download daymet, aggregate cells, save as raster

# process Daymet data to extend N and S boundary of climate data past CONUS
library(raster)
library(prism)
library(daymetr)
library(maps)
library(maptools)
library(dplyr)

weather_path <- "daymet/2009"

class(tile_outlines)

# select which tiles you want by coordinates
test <- tile_outlines[which(tile_outlines$Ymax < 55 & 
                              tile_outlines$Ymax > 24 & 
                              tile_outlines$XMax > -140), ]
plot(test)

# cut a few loose tiles
test <- test[-which(test$TileID %in% c(11038, 11548, 11940, 11941, 10651, 12664)), ]

# list all the tiles included in canada and mexico to use
dl_tiles <- unique(test$TileID)


# download daymet netcdf tiles, reproject to PRISM, save as R Raster (.grd)

yrs <- 2009
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

# for some reason, the above loop hangs up and says its still working
# even when all tiles have been downloaded succesfully
# have to kill R sessions in command line



if(exists("cl")){
  stopCluster(cl)
}


# test if all tiles downloaded and converted successfully
grdfiles <- list()
for(year in 1:length(yrs)){
  for(tile in 1:length(dl_tiles)){
    
    yr <- yrs[year]
    id <- dl_tiles[tile]
    grdfiles[[length(grdfiles)+1]] <- paste0(weather_path, "/tmax_", yr, "_", id, ".grd")
    grdfiles[[length(grdfiles)+1]] <- paste0(weather_path, "/tmin_", yr, "_", id, ".grd")
    
  }
}

grddf <- data.frame(row = 1:length(grdfiles), name = unlist(grdfiles))

exist <- list.files(weather_path, full.names = TRUE)
existdf <- data.frame(name = exist[grep(pattern = ".grd", x = exist, fixed = TRUE)])

rerundf <- grddf %>% 
  anti_join(existdf, by = "name")

cl <- makePSOCKcluster(6)
registerDoParallel(cl)

# rerun ones that didn't download
outfiles <- foreach(rr = 1:nrow(rerundf),
                    .packages= c("raster", "daymetr"),
                    .inorder = FALSE)%dopar%{
                      
                      tm <- as.character(rerundf$name[rr])
                      tempvar <- substr(tm, 8, 11)
                      id <- substr(tm, 18, 22)
                      yr <- substr(tm, 13, 16)
                      
                      if(!file.exists(tm)){
                        
                        download_daymet_tiles(tiles = id,
                                              start = yr, end = yr, path = weather_path, param = tempvar,
                                              silent = FALSE, force = FALSE)
                      }
                      
                      tm_nc <- sub(".grd", replacement = ".nc", x = tm, fixed = TRUE)
                      
                      # reproject to prism CRS
                      prism_res <- 25.875 / 621
                      newproj <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
                      
                      
                      x <- brick(tm_nc, varname = tempvar, level = 1)
                      projectRaster(x, res = prism_res, crs = newproj, method="bilinear", 
                                    filename = tm)
                      
                      # remove daymet netcdf
                      if(file.exists(tm)){
                        file.remove(tm_nc)
                      }
                      
                    }



# stitch together canada and mexico tiles by year and tmin/tmax
# has to first stitch together each row of tiles, otherwise uses too much memory
library(stringr)
library(purrr)
tiles <- list.files(weather_path, recursive = FALSE, full.names = TRUE)
tiles <- tiles[grep(pattern = ".grd", x = tiles, fixed = TRUE)]
tile_params <- str_split_fixed(tiles, pattern = coll("."), n = 2)[, 1]
# below, needs adjustment for number of slashes in weather_path
# TODO: change this so automatic
tile_params <- str_split_fixed(tile_params, pattern = coll("/"), n = 3)[, 3]

temp <- str_split_fixed(tile_params, pattern = coll("_"), n = 3)
temps <- temp[, 1]
tileyears <- temp[, 2]
tilenames <- temp[, 3]
tilerow <- substr(tilenames, 1, 3)

# rmfiles <- tiles[grep(pattern = "11548", x = tiles, fixed = TRUE)]

# cl <- makePSOCKcluster(2)
# registerDoParallel(cl)

outfiles <- foreach(year = 1:length(yrs),
                    .packages= c("raster", "daymetr"),
                    .inorder = FALSE)%:%
  foreach(varname = 1:2,
          .packages= c("raster", "daymetr"),
          .inorder = FALSE)%dopar%{
                      
                      varname = 1
                      yr <- yrs[year]
                      tvar <- unique(temps)[varname]
                      index <- tiles[which(temps == tvar &
                                                  tileyears == yr)]
                      trows <- tilerow[which(temps == tvar &
                                               tileyears == yr)]
                      tfile <- paste0("daymet/", tvar, "_", yr, ".grd")
                      
                      rasterOptions(tolerance = 0.1)
                      
                      AutoShift <- function(ras){
                        r <- raster::shift(ras, x = -1*origin(ras)[1], y = -1*origin(ras)[2])
                      }
                      
                      if(!file.exists(tfile)){
                        row_raster <- list()
                        # merge by daymet row as 1st step
                        for (d in 1:length(unique(trows))){
                          drow <- unique(trows)[d]
                          index_tmp <- index[which(trows == drow)]
                          rlist <- lapply(index_tmp, brick)
                          rlist <- lapply(rlist, AutoShift)
                          
                          updating_raster <- rlist[[1]]
                          
                          for (i in 2:length(rlist)) {
                            next_raster <- rlist[[i]]
                            updating_raster <- merge(x = updating_raster, 
                                                     y = next_raster)
                          }
                          row_raster[[d]] <- updating_raster                          
                        } # close row loop
                        
                        row_raster <- lapply(row_raster, AutoShift)
                        updating_raster <- row_raster[[1]]
                        
                        for (i in 2:length(row_raster)) {
                          next_raster <- row_raster[[i]]
                          updating_raster <- merge(x = updating_raster, 
                                                   y = next_raster, 
                                                   filename = tfile,
                                                   overwrite = TRUE)
                        }
                      } # close if file exists
                        
                      # remove tiled daymet files 
                      outfiles <- gsub(pattern = ".grd", replacement = ".gri", x = index, fixed = TRUE)
                      lapply(outfiles, FUN = file.remove)
                      lapply(index, FUN = file.remove)
                      
                      }

if(exists("cl")){
  stopCluster(cl)
}

