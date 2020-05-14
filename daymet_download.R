# download daymet, aggregate cells, save as raster

# process Daymet data to extend N and S boundary of climate data past CONUS
library(raster)
library(prism)
library(daymetr)
library(maps)
library(maptools)
library(dplyr)
library(ncdf4)
# library(gdalUtils)
# library(velox)
# library(stars)

weather_path <- "daymet/process"

class(tile_outlines)

# select which tiles you want by coordinates
test <- tile_outlines[which(tile_outlines$Ymax < 63 & 
                              tile_outlines$Ymax > 24 & 
                              tile_outlines$XMax > -140), ]
plot(test)

# cut a few loose tiles
test <- test[-which(test$TileID %in% c(11038, 11548, 11940, 11941, 10651, 12664, 13559)), ]

# list all the tiles included in canada and mexico to use
dl_tiles <- unique(test$TileID)


# download daymet netcdf tiles, reproject to PRISM, save as R Raster (.grd)

yrs <- c(2016:2018)
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
                
                if(length(rlist) > 1){
                  for (i in 2:length(rlist)) {
                    next_raster <- rlist[[i]]
                    updating_raster <- merge(x = updating_raster, 
                                             y = next_raster)
                  }
                }
                row_raster[[d]] <- updating_raster                          
              } # close row loop
              
              row_raster <- lapply(row_raster, AutoShift)
              updating_raster <- row_raster[[1]]
              
              for (i in 2:length(row_raster)) {
                print(i)
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




# GRAVEYARD

# fixing gaps in rasters when merged
# aggregating pre-merge doesn't work well
for (tt in c(13522:13526)){
  download_daymet_tiles(tiles = tt,
                      start = 2016, end = 2016, path = "daymet", param = "tmin",
                      silent = FALSE, force = FALSE)
}

fs <- list.files("daymet", pattern = ".nc", full.names = TRUE)

rlist <- list()
for (i in 1:length(fs)){
  rlist[[i]] <- brick(fs[i], varname = "tmin", level = 1)
  }

# reproject to prism CRS
prism_res <- 25.875 / 621
newproj <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"

tfile <- "daymet/test.grd"

rlist$filename <- tfile
rlist$overwrite <- "TRUE"

rmerg <- do.call(raster::merge, rlist)

rmerg_res <- projectRaster(rmerg, res = prism_res, crs = newproj, method="bilinear", 
                           filename = "daymet/testres.grd")


rlist_res <- list()
for (i in 1:5){
  rlist_res[[i]] <- projectRaster(rlist[[i]], res = prism_res, crs = newproj, method = "bilinear")
}

AutoShift <- function(ras){
  r <- raster::shift(ras, x = -1*origin(ras)[1], y = -1*origin(ras)[2])
}
rlist_res <- lapply(rlist_res, AutoShift)

rmerg_res2 <- do.call(raster::merge, rlist_res)
plot(rmerg_res2[[1]])

# Just download whole goddam continent (3gb)
dl <- "https://thredds.daac.ornl.gov/thredds/fileServer/ornldaac/1328/2012/daymet_v3_tmax_2012_na.nc4"
fn <- stringr::str_split_fixed(dl, pattern = stringr::coll("/"), n = 9)[,9]

download.file(dl, 
              destfile = paste0("../photovoltinism/daymet/", fn))

r <- brick("daymet/daymet_v3_tmin_2017_na.nc4", varname = "tmin", level = 1)[[1]]
# reproject to prism CRS
prism_res <- 25.875 / 621
newproj <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
rc <- crop(r, extent(r)*.9)

tt <- system.time({
  # ragg <- gdalwarp(r, f2 <- tempfile(fileext='.tif'), output_Raster = TRUE, multi = TRUE, r = 'med', tr=res(r)*4)
  ragg <- raster::aggregate(r, fact = 500, fun = max, expand = FALSE, na.rm = TRUE)
rproj <- projectRaster(r, res = prism_res, crs = newproj, method="bilinear")
                           # filename = "daymet/tmax_2017.grd", overwrite = TRUE)
})


# Steps
# 1. Continent brick
# 2. Process in parallel by day
# 3. projectRaster by day slice and write file
# 4. Rebrick files


library(foreach)
library(doParallel)
library(stringr)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)


fs <- list.files("daymet", pattern = '.nc4', full.names = TRUE)

outfiles <- foreach(f = 1:length(fs),
          .packages= c("raster", "stringr"),
          .inorder = FALSE)%:%
  foreach(day = 1:365,
          .packages= c("raster", "stringr"),
          .inorder = FALSE)%dopar%{
            
            fl <- fs[f]
            newname <- stringr::str_split_fixed(string = fl, pattern = "_", n = 3)[,3]
            varname <- stringr::str_split_fixed(string = newname, pattern = "_", n = 3)[,1]
            
            newname <- stringr::str_split_fixed(string = newname, pattern = "_", n = 3)[,1:2]
            e <- extent(-142, -52, 24, 60)
            
            fslice <- brick(fl, varname = varname, level = 1)[[day]]
            proj4string(fslice) <- CRS("+proj=lcc +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +lat_1=25 +lat_2=60 +datum=WGS84 +units=m")
            
            # reproject to prism CRS
            prism_res <- 25.875 / 621
            newproj <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
            
            outr <- crop(projectRaster(fslice, res = prism_res, crs = newproj, method="bilinear"), e)
            
            dat <- gsub(pattern = "X", replacement = "", x = names(fslice), fixed = TRUE) 
            dat <- gsub(pattern = ".", replacement = "", x = dat, fixed = TRUE) 
            
            fname <-  paste0("/home/tyson/REPO/photovoltinism/daymet/", newname[2], "/DAYMET_", varname, "_stable_4kmD1_", dat, "_bil.grd")
            writeRaster(outr, filename = fname, overwrite = TRUE)
          }
            
if(exists("cl")){
  stopCluster(cl)
}
# 
# 
# newfiles <- list.files("daymet", pattern = ".grd", full.names = TRUE)
# for (f in 1:length(fs)){
#   fl <- fs[f]
#   newname <- stringr::str_split_fixed(string = fl, pattern = "_", n = 3)[,3]
#   newname <- stringr::str_split_fixed(string = newname, pattern = "_", n = 3)[,1:2]
#   dayfiles <- newfiles[grep(pattern = paste(newname[1], newname[2], sep = "_"), x = newfiles, fixed = TRUE)]
# 
#   yearbrick <- brick(as.list(dayfiles))
#   writeRaster(yearbrick, filename = paste0("daymet/", newname[1], "_", newname[2], ".grd"),
#               overwrite = TRUE)
#   unlink(dayfiles)
#   }
# 
# r <- brick("daymet/tmin_2017.grd")
# 
# 
# # Make Daymet files like PRISM directory??
# fs <- list.files("/home/tyson/REPO/photovoltinism/daymet", pattern = ".grd", full.names = TRUE, recursive = TRUE)
# 
# todo <- fs[4]
# tmp <- brick(todo)
# 
# for(day in 1:365){
#   tmp2 <- tmp[[day]]
#   dat <- gsub(pattern = "X", replacement = "", x = names(tmp2), fixed = TRUE) 
#   dat <- gsub(pattern = ".", replacement = "", x = dat, fixed = TRUE) 
#   
#   writeRaster(tmp2, filename = paste0("/home/tyson/REPO/photovoltinism/daymet/2018/DAYMET_tmin_stable_4kmD1_", dat, "_bil.bil"))
#   }
# 
# 
# # add prj files???
# fs <- list.files("/home/tyson/REPO/photovoltinism/daymet", pattern = "bil.bil", full.names = TRUE, recursive = TRUE)
# for (i in 1:length(fs)){
#   tmp <- raster(fs[i])
#   crs(tmp) <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
#   writeRaster(tmp, filename = fs[i], overwrite = TRUE, prj = TRUE)
# }
# 
# 
# 
# 
# region <-spTransform(region, CRS(proj4string(template)))
# reg.points = fortify(region, region="name_en")
# reg.df = left_join(reg.points, region@data, by = c("id" = "name_en")) %>% 
#   filter(geonunit %in% c("Canada", "Mexico", "United States of America"), long < -50)
#   
# 
# 
