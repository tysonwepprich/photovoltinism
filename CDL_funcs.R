########                       BEGINNING function definitions                      #########
#### if then else raster function [sim. to GRASS r.mapcalc if(x,a,b)]:
Cond=function(condition, trueValue, falseValue){
  return(condition * trueValue + (!condition)*falseValue)
}

#####
# Simple way to assign weights for substages
# Choose number of substages and the decimal (0, 1) of population to include
# Outputs standard deviations from mean and weights for substage
SubstageVals <- function(numstage, perc){
  low <- qnorm((1 - perc)/2)
  high <- qnorm(1 - (1 - perc) / 2)
  bounds <- seq(low, high, length.out = numstage + 1)
  means <- (bounds[1:numstage] + bounds[2:(numstage + 1)]) / 2
  weights <- diff(pnorm(bounds), lag = 1)
  return(data.frame(means, weights))
}

#### DD Calculation Methods: Tmean (no upper threshold, Tmax+Tmin AVG w/horiz upper threshold, 
#   and Single Triangle w/horiz. upper threshold

####  Simple Mean Temp DD Calc method: ((tmean > LDT) * (tmean - LDT))
# same result as  max((tmax + tmin)/2 - LDT,0)  
# so no need for tmean PRISM data. 
SimpDD=function(tmax,tmin,LDT){
  return(max((tmax + tmin)/2 - LDT,0))
}

#### Averaging DD Calc method (max + min/2 - tlow) but with horizontal (substitution) upper threshold:
AvgDD=function(tmax, tmin, LDT, UDT){
  return(Cond(tmax < LDT, 0, Cond(tmin > UDT, 0, Cond(tmax > UDT, (UDT + tmin)/2 - LDT, Cond((tmax + tmin)/2-LDT < 0,0,(tmax + tmin)/2 - LDT)))))
}

#### Single triangle with upper threshold (Sevachurian et al. 1977) - also a good substitution for 
#  single sine method
TriDD=function(tmax, tmin, LDT, UDT){
  Tmp1=6*((tmax-LDT)*(tmax-LDT))/(tmax-tmin)
  Tmp2=6*((tmax-UDT)*(tmax-UDT))/(tmax-tmin)
  Cond(tmax < LDT,0,
       Cond(tmin >= UDT,UDT-LDT,
            Cond((tmax < UDT) & (tmin <= LDT), Tmp1/12,
                 Cond((tmin <= LDT) & (tmax >= UDT), (Tmp1-Tmp2)/12,
                      Cond((tmin > LDT) & (tmax >= UDT), 6*(tmax+tmin-2*LDT)/12 - (Tmp2/12),
                           Cond((tmin > LDT) & (tmax < UDT), 6*(tmax+tmin-2*LDT)/12,0))))))
} 




# photoperiod functions

photoperiod <- function(lat, doy, p = 1.5){
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



# from https://stackoverflow.com/questions/29784829/r-raster-package-split-image-into-multiples

# The function spatially aggregates the original raster
# it turns each aggregated cell into a polygon
# then the extent of each polygon is used to crop
# the original raster.
# The function returns a list with all the pieces
# in case you want to keep them in the memory. 
# it saves and plots each piece
# The arguments are:
# raster = raster to be chopped            (raster object)
# ppside = pieces per side                 (integer)
# save   = write raster                    (TRUE or FALSE)
# plot   = do you want to plot the output? (TRUE or FALSE)
SplitRas <- function(raster,ppside,save,plot){
  h        <- ceiling(ncol(raster)/ppside)
  v        <- ceiling(nrow(raster)/ppside)
  agg      <- aggregate(raster,fact=c(h,v))
  agg[]    <- 1:ncell(agg)
  agg_poly <- rasterToPolygons(agg)
  names(agg_poly) <- "polis"
  r_list <- list()
  for(i in 1:ncell(agg)){
    e1          <- extent(agg_poly[agg_poly$polis==i,])
    r_list[[i]] <- crop(raster,e1)
  }
  if(save==T){
    for(i in 1:length(r_list)){
      ii <- formatC(i, width = 2, format = "d", flag = "0")
      writeRaster(r_list[[i]],filename=paste("SplitRas", ii, sep=""),
                  format="GTiff", datatype="FLT4S", overwrite=TRUE)  
    }
  }
  if(plot==T){
    par(mfrow=c(ppside,ppside))
    for(i in 1:length(r_list)){
      plot(r_list[[i]],axes=F,legend=F,bty="n",box=FALSE)  
    }
  }
  return(r_list)
}


# from https://www.r-bloggers.com/notes-on-multivariate-gaussian-quadrature-with-r-code/
## compute Gauss-Hermite quadrature points and weights
## for a one-dimensional integral.
## points -- number of points
## interlim -- maximum number of Newton-Raphson iterations

hermite <- function (points, z) {
  p1 <- 1/pi^0.4
  p2 <- 0
  for (j in 1:points) {
    p3 <- p2
    p2 <- p1
    p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
  }
  pp <- sqrt(2 * points) * p2
  c(p1, pp)
}

gauss.hermite <- function (points, iterlim = 50) {
  x <- w <- rep(0, points)
  m <- (points + 1)/2
  for (i in 1:m) {
    z <- if (i == 1) 
      sqrt(2 * points + 1) - 2 * (2 * points + 1)^(-1/6)
    else if (i == 2) 
      z - sqrt(points)/z
    else if (i == 3 || i == 4) 
      1.9 * z - 0.9 * x[i - 2]
    else 2 * z - x[i - 2]
    for (j in 1:iterlim) {
      z1 <- z
      p <- hermite(points, z)
      z <- z1 - p[1]/p[2]
      if (abs(z - z1) <= 1e-15) 
        break
    }
    if (j == iterlim) 
      warning("iteration limit exceeded")
    x[points + 1 - i] <- -(x[i] <- z)
    w[i] <- w[points + 1 - i] <- 2/p[2]^2
  }
  r <- cbind(x * sqrt(2), w/sum(w))
  colnames(r) <- c("Points", "Weights")
  r
}



CombineMaps <- function(rasfiles, tmpdir, newdir){
  
  rasfiles <- rasfiles[grep(pattern = ".grd", x = rasfiles, fixed = TRUE)]
  
  # for each sim with unique information to save
  ls <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,1])
  maps <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,2])
  sims <- unique(stringr::str_split_fixed(rasfiles, pattern = "_", 3)[,3])
  sims <- sort(gsub(pattern = ".grd", replacement = "", x = sims))
  
  reslist <- list()
  # for each sim with unique information to save
  for (i in ls){
    fs <- rasfiles[grep(pattern = i, x = rasfiles, fixed = TRUE)]
    for (j in sims){
      fs2 <- fs[grep(pattern = j, x = fs, fixed = TRUE)]
      
      bricklist <- list()
      for (m in 1:length(fs2)){
        bricklist[[m]] <- brick(fs2[m])
      }
      
      bricklist$filename <- paste(newdir, 
                                  paste(i, j, "all", sep = "_"), sep = "/")
      bricklist$overwrite <- TRUE
      test <- do.call(merge, bricklist)
      reslist[[length(reslist) + 1]] <- bricklist$filename
    }
  }
  return(reslist)
  cleanup <- list.files(path = tmpdir)
  cleanup <- cleanup[-grep("all", x = cleanup)]
  lapply(cleanup, FUN = file.remove) # CAREFUL HERE!
  
}


