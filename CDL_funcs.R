########                       BEGINNING function definitions                      #########
#### if then else raster function [sim. to GRASS r.mapcalc if(x,a,b)]:
Cond=function(condition, trueValue, falseValue){
  return(condition * trueValue + (!condition)*falseValue)
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

