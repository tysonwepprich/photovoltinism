# photoperiod functions

photoperiod <- function(lat, doy, p = 6){
  theta <- 0.2163108 + 
    2 * atan(0.9671396 * tan(0.00860 * (doy - 186)))
  phi <- asin(0.39795 * cos(theta))
  D <- 24 - (24 / pi) * acos(
    (sin(p * pi / 180) + sin(lat * pi / 180) * sin(phi))/
      (cos(lat * pi / 180) * cos(phi))
  )
}



# calculate photoperiod values for raster
# day of year from 'd', index of daily PRISM data used for DD accumulation
RasterPhoto <- function(rast, doy){
  xy <- coordinates(rast)
  hours <- photoperiod(xy[, 2], doy)
  outrast <- setValues(rast, hours)
  rast[!is.na(rast)] <- 0
  outrast <- rast + outrast
}

