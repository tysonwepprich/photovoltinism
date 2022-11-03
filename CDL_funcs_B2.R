# Function definitions----

'%!in%' <- function(x,y)!('%in%'(x,y))

isFALSE <- function (x) {
  identical(x, FALSE)
}

logit <- function(p){log(p/(1-p))}

#### if then else raster function [sim. to GRASS r.mapcalc if(x,a,b)]:
Cond<- function(condition, trueValue, falseValue){
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

# Simple Mean Temp DD Calc method: ((tmean > LDT) * (tmean - LDT))
# same result as  max((tmax + tmin)/2 - LDT,0)  
# so no need for tmean PRISM data. 
SimpDD <- function(tmax,tmin,LDT) {
  return(max((tmax + tmin)/2 - LDT,0))
}

# Averaging DD Calc method (max + min/2 - tlow) but with horizontal (substitution) upper threshold:
AvgDD <- function(tmax, tmin, LDT, UDT) {
  return(Cond(tmax < LDT, 0, Cond(tmin > UDT, 0, Cond(tmax > UDT, (UDT + tmin)/2 - LDT, Cond((tmax + tmin)/2-LDT < 0,0,(tmax + tmin)/2 - LDT)))))
}

# Single triangle with upper threshold (Sevachurian et al. 1977) - also a good substitution for 
# single sine method
TriDD <- function(tmax, tmin, LDT, UDT) {
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
                   "NW_SMALL"      = extent(-125.1, -117.9, 38.5, 49.2),
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
  pfile <- str_split(string = prismfiles, pattern = coll("/"), numsplits) %>% purrr::map(numsplits)
  qa <- str_split(string = pfile, pattern = coll("_"), 6) %>% purrr::map(3) %>% unlist()
  
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
    dplyr::mutate(datarank = 1:n()) %>% 
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

#### BRITTANY'S FUNCTIONS ####
# Save raster - simply the "writeRaster" function from raster library but prints progress in log file
SaveRaster <- function(r, d, sim, outnam) {
  writeRaster(r, file = paste0(outnam, "_", sim, "_", d), format="GTiff", overwrite=TRUE)
  print(noquote(paste0("Saving raster: ", outnam,"_", d,".tif")))
}

# Convert raster (TIF file) into a data frame for plotting
ConvDF <- function(rast) {
  spdf <- as(rast, "SpatialPixelsDataFrame")
  df <- as.data.frame(spdf)
  colnames(df) <- c("value", "x", "y")
  return(df)
}

# Classify data into bins so it can be visualized in categories (e.g., 1-10, 11-20, 21-30) 
Cut_bins <- function(f) {
  f$value_orig <- f$value # keep old value so can sort factors against it
  f$value <- ifelse(f$value < 0.1, round_any(f$value, 0.01),  # round to hundredth dec
                    ifelse(f$value < 1 & f$value > 0.1, round_any(f$value, .1), # round to tenth
                           round_any(f$value, 1)))
  #f$value <- ifelse(f$value < 1, round_any(f$value, .1), round_any(f$value, 1)) # round decimal places appropriately
  f2 <- f %>% mutate(value = cut(value, breaks = 10, dig.lab=4)) # cut values into values
  f2$value <- gsub(pattern = "[()]|\\[|\\]", replacement = "", f2$value) # remove brackets and parentheses from output value
  f2$value <- gsub(pattern = "[,]", replacement = "-", f2$value) # replace comma with a dash
  f2$value <- factor(f2$value, levels = unique(f2$value[order(f2$value_orig)])) # Order by values so plots in numerical order
  return(f2)
}

# Plot theme to be used for all summary (PNG) maps
#mytheme <- theme(plot.title=element_text(size=10, hjust=0.5, vjust=-3),  # Specify title attributes
                # plot.subtitle = element_text(size=6.5, hjust=0.5, vjust=-3, lineheight = 1.25),
                # plot.margin = unit(c(0,0.5,0,0.5),"line"),
                # legend.text=element_text(size = 6), # Specify legend attributes
                # legend.position="right", 
                # legend.justification="left", 
                # legend.margin=margin(0,0,0,0),
                # legend.key.width = unit(0.75,"line"), legend.key.height = unit(0.4,"line"))

# Base features used for all summary (PNG) maps in "PlotMap" function
Base_map <- function(df) {
  p <- ggplot2::ggplot(states, ggplot2::aes(x=long, y=lat)) + 
    ggplot2::geom_raster(data=df, ggplot2::aes(x=x, y=y, fill=value)) + # plot raster (raster values are now in a data frame)
    ggplot2::geom_path(ggplot2::aes(group=group), color="gray20", lwd = 0.25) +
    ggthemes::theme_map(base_size = 7)
}

# Creates a gradient color ramp consisting of 4 colors
Colfunc <- function(x,y) {
  colorRampPalette(c(x,y))(4)
}

# Deal with 0 vs non-zero values when plotting chill and heat stress unit rasters
Stress_Val_Conv <- function(x) {
  if(all(x$value == 0)){
    x2 <- x
  } else {
    x2 <- Cut_bins(x)
  }
  x2$value <- factor(x2$value)
  return(x2)
}

PlotMap2 <- function(r,d,titl,lgd,outfl) {
  df <- ConvDF(r) # convert raster to data frame
  nam <- deparse(substitute(r)) # Get name of raster, capitalize first letter
 if (nam == "Lifestage") {
   df <- dplyr::mutate(df, value = ifelse(value==1,"OA",ifelse(value==2,"egg",ifelse(value==3,"larva",ifelse(value==4,"pupa","adult")))))
    df$value <- factor(df$value, levels = c("OA","egg","larva", "pupa", "adult"))   
    col_key <- cbind(data.frame("cols" = c("gray70", "#E41A1C","#377EB8","#4DAF4A","#984EA3"), stringsAsFactors = FALSE),
                     data.frame("value" = c("OA","egg","larva", "pupa", "adult")))
    col_df <- suppressWarnings(suppressMessages(dplyr::semi_join(col_key,df,by="value"))) # join the colors to the data
    cols <- setNames(as.character(col_df$cols),col_df$value) # convert col_df to a named vector
    # make the plot
    p <- Base_map(df) +
      ggplot2::scale_fill_manual(values=cols, name = paste0(lgd)) +
      #labs(title = paste(sp,titl), subtitle = str_wrap(subtitl, width = 75)) +
      mytheme
    # number of generations
  } else if (nam == "NumGen") {
    df$value <- factor(df$value) 
    # make the color ramp - enough for 15 generations
    col_key <- cbind(data.frame("cols" = c(RColorBrewer::brewer.pal(8,"Set1"),RColorBrewer::brewer.pal(8,"Pastel2")), stringsAsFactors = FALSE),
                     data.frame("value" = c(as.character(0:15)), stringsAsFactors = FALSE))
    col_df <- suppressWarnings(suppressMessages(dplyr::semi_join(col_key,df,by="value"))) # join the colors to the data
    cols <- setNames(as.character(col_df$cols),col_df$value) # convert col_df to a named vector
    # make the plot
    p <- Base_map(df) +       
      ggplot2::scale_fill_manual(values=cols, name = paste0(lgd)) +
      #labs(title = paste(sp,titl), subtitle = str_wrap(subtitl, width = 75)) +
      mytheme
  }
    # save the plot, or else report that there was an error and skip - see "rmessages.txt" for error report
    ggplot2::ggsave(p,file=paste0(outfl,"_",d,"_",".png"), width = 6, height = 3, units = c('in'), dpi = 200) # save the plot
    print(noquote(paste0("Saving summary map: ",outfl,"_",d,".png"))) # print progress in log file
}


