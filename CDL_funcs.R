# Function definitions----

'%!in%' <- function(x,y)!('%in%'(x,y))

isFALSE <- function (x) {
  identical(x, FALSE)
}

logit <- function(p){log(p/(1-p))}

#### if then else raster function [sim. to GRASS r.mapcalc if(x,a,b)]:
Cond=function(condition, trueValue, falseValue){
  return(condition * trueValue + (!condition)*falseValue)
}


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
  # if tmax == tmin, gives NaN
  tmax <- Cond(tmax == tmin, tmax + 0.01, tmax)
  
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
RasterPhoto <- function(rast, doy, perc_twilight){
  xy <- coordinates(rast)
  p <- perc_twilight * 6 / 100
  hours <- photoperiod(xy[, 2], doy, p)
  outrast <- setValues(rast, hours)
}



SubstageDistrib <- function(dist, numstage, perc){
  # this is an approximation from GAM predictions, could make more precise
  # with more predicted values or interpolation, but doesn't seem necessary
  ReturnClosestValue <- function(dist, xval){
    out <- dist$CDF[which(dist$x > xval)[1]]
  }
  
  low <- (1 - perc)/2
  high <- 1 - (1 - perc) / 2
  low <- dist$x[which(dist$CDF > low)[1]]
  high <- dist$x[which(dist$CDF > high)[1]]
  
  bounds <- seq(low, high, length.out = numstage + 1)
  means <- (bounds[1:numstage] + bounds[2:(numstage + 1)]) / 2
  weights <- diff(sapply(X = bounds, FUN = ReturnClosestValue, dist = dist), lag = 1)
  return(data.frame(means, weights))
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



# Add new extent definitions here for use in models and plots
assign_extent <- function(region_param = c("CONUS", "NORTHWEST", "NW_SMALL", "OR","EAST", "TEST", "WEST", "SOUTHWEST")){
  REGION <- switch(region_param,
                   "ALL"          = extent(-133.5, -51.5, 24, 54.7),
                   "SW_MEX"       = extent(-122, -95, 24, 48),
                   "CONUS"        = extent(-125.0,-66.5,24.0,50.0),
                   "NORTHWEST"    = extent(-125.1,-103.8,40.6,49.2),
                   "NW_SMALL"      = extent(-125.1, -116, 38.5, 49.2),
                   "OR"           = extent(-124.7294, -116.2949, 41.7150, 46.4612),
                   "TEST"         = extent(-124, -122.5, 44, 45),
                   "WEST"         = extent(-121.4, -102.5, 31.23, 47.33),
                   "SOUTHWEST"    = extent(-120.17, -108.25, 31.5, 42.3),
                   "EAST"         = extent(-93, -71, 35, 47))
  return(REGION)
}

# Take .bil files from PRISM yearly directories
# Return best data for each day
# Remove leap year if not needed
ExtractBestPRISM <- function(prismfiles, yr, leap = "keep"){
  numsplits <- str_count(string = prismfiles[1], pattern = "/")
  pfile <- str_split(string = prismfiles, pattern = coll("/"), numsplits) %>% map(numsplits)
  qa <- str_split(string = pfile, pattern = coll("_"), 6) %>% map(3) %>% unlist()
  
  dates <- regexpr(pattern = "[0-9]{8}", text = prismfiles)
  
  df <- data.frame(dates = regmatches(prismfiles, dates),
                   quality = substr(qa, start=1, stop=4),
                   rownum = 1:length(qa))
  
  # sorting backwards matches data hierarchy
  # stable > provisional > early > 10yr
  df2 <- df %>% 
    mutate(quality = as.character(quality)) %>% 
    group_by(dates) %>% 
    dplyr::arrange(desc(quality)) %>% 
    mutate(datarank = 1:n()) %>% 
    filter(datarank == 1)
  
  if (leap == "remove"){
    # still has issue with leap year
    if (yr %% 4 != 0) {
      df2 <- df2[!grepl(pattern = paste0(yr, "0229"), x = df2$dates), ]
    }
  }
  
  best <- prismfiles[df2$rownum]
  dates <- regexpr(pattern = "[0-9]{8}", text = best)
  
  fileorder <- order(regmatches(best, dates))
  prismfiles <- best[fileorder]
  return(prismfiles)
}



# From TrenchR on github: https://github.com/trenchproject/TrenchR/blob/master/R/DDFunctions.R
#' Calculate degree days
#' 
#' 
#' @details Calculate degree days
#' @description This function allows you to calculate degree days using single or double sine wave and single or double triangulation approximation. Source: http://ipm.ucanr.edu/WEATHER/ddfigindex.html. Double methods assume symmetry, that is that next day's thermal minima is equal to previous day. Double Sine wave approximation of degree days from Allen JC. 1976.  A Modified Sine Wave Method for Calculating Degree Days. Environmental Entomology 5:388–396. 
#' @param T_min Minimum temperature of the day (°C)
#' @param T_max Maximum temperature of the day (°C)
#' @param LDT lower developmental threshold (°C).
#' @param UDT upper developmental threshold (°C).
#' @param method type of method being used. Current choices: "single.sine","double.sine", "single.triangulation" or "double.triangulation".
#' @return degree days (°C)
#' @keywords degree days
#' @family microclimate functions
#' @export
#' @examples
#' \dontrun{
#' degree_days(T_min=7, T_max=14, LDT=12, UDT=33, method="single.sine")
#' degree_days(T_min=7, T_max=14, LDT=12, UDT=33, method="single.triangulation")
#' }
#' 

degree_days=function(T_min,T_max,LDT=NA,UDT=NA, method="single.sine"){
  
  stopifnot(T_max>=T_min, method %in% c("single.sine","double.sine", "single.triangulation", "double.triangulation"))
  
  #amplitude
  alpha=(T_max-T_min)/2 
  dd = 0
  #Single sine calculation
  if (method == "single.sine") {
    
    
    if (T_min >= UDT && T_max > UDT) { # entirely above both thresholds
      dd = (T_max - T_min)
    } else  if ( T_min > LDT  && T_max > UDT) { #Intercepted by upper threshold 
      theta2=asin((UDT-(T_max+T_min)/2)/alpha)
      dd = 1 / pi * (((T_max + T_min) / 2 - LDT) * (theta2 + pi / 2) + (UDT - LDT) *
                       (pi / 2 - theta2) - alpha * cos(theta2))
    } else  if (T_min < LDT &&  T_max > UDT ) {  #Intercepted by both thresholds
      theta2=asin((UDT-(T_max+T_min)/2)/alpha)
      theta1=asin((LDT-(T_max+T_min)/2)/alpha)
      dd = 1 / pi * ((  ((T_max + T_min) / 2) - LDT) * (theta2 - theta1) + alpha * (cos(theta1) -
                                                                                      cos(theta2)) + (UDT - LDT) * (pi / 2 - theta2))
    } else if (T_min > LDT &&  T_max < UDT ) { #Entirely between both thresholds
      dd = ((T_max + T_min) / 2) - LDT
    } else if (T_min < LDT && T_max > LDT) {  # intercepted by LDT  
      #theta1=asin((LDT-(T_max+T_min)/2)/alpha)
      # credit - http://stackoverflow.com/questions/20998460/unexpected-behavior-in-asin
      #It's a floating point issue. The way floating point numbers work is that all 
      #numbers need to be mapped to the nearest one which can be expressed as a finite 
      #sum of powers of two and this may lead to small inaccuracies in the expected output 
      #and can be dependent upon how the numbers are calculated
      theta1= asin(pmax(-1,pmin(1,(LDT-(T_max+T_min)/2)/alpha)))
      dd = 1 / pi * (( ((T_max + T_min) / 2) - LDT) * ( (pi / 2) - theta1) + alpha * cos(theta1))
    } else if (T_min < LDT && T_max <= LDT) { # entirely below both thresholds
      dd = 0
    }
  } #end single sine method
  
  #double sine calculation
  if (method == "double.sine") {
    
    if (T_min >= LDT && T_max > UDT) { # entirely above both thresholds
      dd = (UDT - LDT) / 2
    } else if (T_min > LDT  && T_max > UDT) { #Intercepted by upper threshold
      theta2=asin((UDT-(T_max+T_min)/2)/alpha)
      dd = 1 / (2 * pi) * (((T_max + T_min) / 2 - LDT) * (theta2 + pi / 2) + (UDT -
                                                                                LDT) * (pi / 2 - theta2) - alpha * cos(theta2))
    } else if (T_min < LDT &&  T_max > UDT) { #Intercepted by both thresholds
      theta2=asin((UDT-(T_max+T_min)/2)/alpha)
      theta1=asin((LDT-(T_max+T_min)/2)/alpha)
      dd = 1 / (2 * pi) * (((T_max + T_min) / 2 - LDT) * (theta2 - theta1) + alpha *
                             (cos(theta1) - cos(theta2)) + (UDT - LDT) * (pi / 2 - theta2))
    } else if (T_min > LDT &&  T_max < UDT) { #Entirely between both thresholds
      dd = 0.5 * ((T_max + T_min) / 2 - LDT)
    } else if (T_min < LDT && T_max > LDT) { # intercepted by LDT
      theta1= asin(pmax(-1,pmin(1,(LDT-(T_max+T_min)/2)/alpha)))
      dd = 1 / (2 * pi) * (((T_max + T_min) / 2 - LDT) * (pi / 2 - theta1) + alpha *
                             cos(theta1))
    } else if (T_min < LDT && T_max <= LDT) { # entirely below both thresholds
      dd = 0
    }
    dd= dd*2
  } #end double sine method
  
  #Single triangulation - with simplified formula
  if (method == "single.triangulation") {
    
    MT = (T_max+T_min)/2
    if (T_min >= UDT && T_max > UDT) { # entirely above both thresholds
      dd = (UDT - LDT)
    } else  if ( T_min > LDT  && T_max > UDT) { #Intercepted by upper threshold 
      dd = (MT-LDT)-((T_max-UDT)^2/((T_max-T_min)*2))
    } else  if (T_min < LDT &&  T_max > UDT ) {  #Intercepted by both thresholds
      dd = ((T_max-LDT)^2-(T_max-UDT)^2)/((T_max-T_min)*2)
    } else if (T_min > LDT &&  T_max < UDT ) { #Entirely between both thresholds
      dd = MT-LDT
    } else if (T_min < LDT && T_max > LDT) {  # intercepted by LDT  
      dd = (T_max-LDT)^2/((T_max-T_min)*2)
    } else if (T_min < LDT && T_max <= LDT) { # entirely below both thresholds
      dd = 0
    }
  } #end single triangulation method
  
  #Double triangulation - with simplified formula
  if (method == "double.triangulation") {
    
    MT = (T_max+T_min)/2
    if (T_min >= UDT && T_max > UDT) { # entirely above both thresholds
      dd = (UDT - LDT)/2
    } else  if ( T_min > LDT  && T_max > UDT) { #Intercepted by upper threshold 
      dd = (MT-LDT)-((T_max-UDT)^2/((T_max-T_min)*4))
    } else  if (T_min < LDT &&  T_max > UDT ) {  #Intercepted by both thresholds
      dd = ((T_max-LDT)^2-(T_max-UDT)^2)/((T_max-T_min)*4)
    } else if (T_min > LDT &&  T_max < UDT ) { #Entirely between both thresholds
      dd = (MT/4)-(LDT/2)
    } else if (T_min < LDT && T_max > LDT) {  # intercepted by LDT  
      dd = (T_max-LDT)^2/((T_max-T_min)*4)
    } else if (T_min < LDT && T_max <= LDT) { # entirely below both thresholds
      dd = 0
    }
    dd= dd*2 
    
  } #end double triangulation method
  
  return(round(dd,2))
}

