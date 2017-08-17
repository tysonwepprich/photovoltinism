#!/usr/bin/Rscript
options(echo=FALSE)
args <- commandArgs(trailingOnly = TRUE)
require(sp)
require(rgdal)
require(raster)
require(foreach) # for parallelized loops
require(doMC)    # parallel backend for foreach
registerDoMC(cores=24)

######## BEGINNING of Header Documentation Section #########
 # ddcalc.R: a degree-day calculator in R which produces grids
 #                   By Dan Upper and Len Coop
 #
 #  Uses grid input and output code from DDRPv24
 #  Uses degree-day calculation methods from our ddcalc.pl, which was
 #     derived from an awk ddcalc, which....

######## END of Header Documentation Section #########

######## BEGINNING function definitions #########
#### Degree-day calculation methods.
source("/usr/local/dds/DDRP/docalcs.R")
#### Raster file selection functions

numQual = new.env(hash=TRUE, parent=emptyenv(), size=10L)
numQual[["10yr"]]=1
numQual[["nmme"]]=10
numQual[["early"]]=21
numQual[["provisional"]]=22
numQual[["stable"]]=23

quality=function(filename){
  extra = 0
  qual = strsplit(filename, split="_")[[1]][3]
  if (grepl("^10yr", qual)) {
      # If it's an average layer ...

      # Extract the most recent of the averaged year and add a
      # fraction of it to the return value so that newer averages look
      # better than older averages.

      # R's regex ops don't give you the chars captured in parens in
      # anything even vaguely resembling a reasonable way.  What I
      # want to do here is $qual =~ /10yr\\d\\d(\\d\\d)/; $extra =
      # $1/100; Doing the closest thing to a straightforward
      # equivalent I can find looks like this:
      #
      # rout = regexpr("10yr\d\d(\d\d)", qual, perl=TRUE)
      # dd = substr(qual, attr(rout, "capture.start"),
      #             attr(rout, "capture.start"> attr(rout, "capture.length")-1)
      # extra = as.integer(dd)/100
      #
      # Instead I'm doing the equivalent of this:
      # ($dd = $qual) =~ s/10yr\d\d(\d\d)//; $extra = $dd/100;

      dd = sub("10yr\\d\\d(\\d\\d)", "\\1", qual, perl=TRUE)
      extra = as.integer(dd)/100
      qual="10yr"
  }
  return (numQual[[qual]] + extra)
}

getBestRaster = function(param, datestr) {
    if (! exists(datestr, envir=filesMap))  # if we don't have that date
        stop(paste("datestr '", datestr, "' not anticipated"))
    # get all the files for the date in question
    files = filesMap[[datestr]]
    # get the ones with the right param
    pfiles = grep(param, files, value=TRUE)  # nothing for that param, date
    if (length(pfiles) == 0)
        stop(paste("parameter '", param,
                   "' not available for datestr '", datestr, "'."))
    # find the best one, as defined by the quality function

    best = pfiles[order(sapply(pfiles, quality, USE.NAMES=FALSE),
                        decreasing=TRUE)[1]]
    #cat(datestr, " / ", param, " -> ", best, "\n")

    # load it and crop it.
    return(crop(raster(best),REGION))
}

#### Write out maps as .tif GIS data files, use GRASS to produce PNG copies for viewing as maps
ConvPNG=function(tif_file){
cat("CONVERSION:",tif_file,'Degree-Days','scrap',regionParam,outputDir,fahrenheit,startDate,endDate,"\n")
        system2('/usr/local/dds/DDRP/dograss_NW_NB_54',args=c("/usr/local/dds/DDRP/geotif2png.sh",tif_file,'Degree-Days','scrap',regionParam,outputDir,fahrenheit,startDate,endDate),stdout="",stderr="",stdin="")
}


writeMap=function(ras){
    prevDir = getwd()
    setwd(outputDir)

    writeRaster(ras, outfilename, format="GTiff", overwrite=TRUE)
    ConvPNG(outfilename)

    setwd(prevDir)
}

######## END function definitions #########

######## BEGINNING of Param Handling Section #########
#### Default values for command line params ####
tlo           = 10          # in C (50F)
thi           = 54.4        # in C (130F)
calcType      = "T1"        # S1, S2, T1, T2, I1, I2, A, G, HD, CD.. endless
startDate     = as.Date("2017-01-01", format="%Y-%m-%d")
                            # start date -- year and biofix             
endDate       = as.Date("2017-12-31", format="%Y-%m-%d")
                            # end date (12-31 for voltinism)
regionParam   = "MIDWEST"   # default REGION to use
normalFlag    = FALSE       # if true, use normals not observations
outdir        = "pid"         # output dir for finding results; NULL means construct the default
outputDir     = "/tmp"      # directory in which outfile will be written
fahrenheit    = FALSE       # convert output to degrees F if true
downscale     = FALSE       # ignored for now

#### Process command line args ####
if (length(args) < 10 || length(args) > 10) {
   cat('
DDRP: Degree-Day, establishment Risk, and Pest Event Mapping System
Required params: tlo          (low threshold, degrees C)
                 thi          (high threshold, degrees C)
                 calcType     (calculation method: e.g., T1 or S2)
                 startDate    (YYYY-MM-DD)
                 endDate      (YYYY-MM-DD)
                 REGION       (CONUS,regions,48 states currently)
                 normalFlag   (use AVG data not observations: true/TRUE/false/FALSE; default false)
                 outdir       (output dir with optional file)
                 fahrenheit   (true/TRUE/false/FALSE; convert to degrees F)
                 downscale    (true/TRUE/false/FALSE; ignore for now)
                 ... (future optional downscaling parameters which are
                      only expected if downscale is true/TRUE).

The default outfile is /tmp/dds{region}_{tlo}_{thi}_{calcType}.png, and
it is ued if the outfile argument is the empty string.  If you specify
an outfile, ".png" will be added to it.  outfile may only contain
[A-Za-z0-9_./].

WARNING: outfile MUST NOT come from an untrusted source.  The directory
portion of outfile is NOT checked for sanity.

Run examples:

ddcalc.R 10 54.4 T1 2016-03-15 2016-12-31 MIDWEST false "" false false
    This creates /tmp/ddsMIDWEST_10_54_T1.png in degrees C

ddcalc.R 10 54.4 S2 2015-08-01 2016-05-31 OR false /tmp/myname true false
    This creates /tmp/myname.png in degrees F

', "\n")  
   q()
} else {
	cat("PARAM LIST:  ")
   print(args)
}


#### Read in command line args ####
tlo           = as.double(args[1])
thi           = as.double(args[2])
calcType     = args[3]
startDate    = as.Date(args[4], format="%Y-%m-%d")
endDate      = as.Date(args[5], format="%Y-%m-%d")
regionParam  = args[6]
if (args[7] == 'true' || args[7] == 'TRUE')    normalFlag = TRUE
#outfile = args[8]
outdir  = args[8]
if (args[9] == 'true' || args[9] == 'TRUE')    fahrenheit = TRUE
# downscale -- if flag true read whatever options are needed

#### Use gsub to make secure cgi-bin args
#: equiv. of perl spp =~ s/[^A-Za-z0-9]//g; ####
# I'm not doing the numerics on the assumption that as.integer() is enough
calcType     = gsub('[^A-Za-z0-9."]', "", calcType, perl=TRUE)
regionParam  = gsub('[^A-Za-z0-9."]', "", regionParam, perl=TRUE)
#outfile      = gsub('[^A-Za-z0-9._/"]', "", outfile, perl=TRUE)
outdir       = gsub('[^A-Za-z0-9."]', "", outdir, perl=TRUE)

#if (FALSE) {
if (TRUE) {
    print(paste("tlo:", tlo))
    print(paste("thi:", thi))
    print(paste("calcType:", calcType))
    print(paste("dates:",startDate, "-", endDate))
    print(paste("region:", regionParam))
    print(paste("normalFlag:", normalFlag))
}

# set the output filename.  Note that we need the .tif ending here;
# conversion to .png happens later.
#if ( is.na(outfile) || outfile == "") {
if ( is.na(outdir) || outdir == "") { # EMPTY - NOT USING FOR WEB USE
    # construct the default output filename
    outfilename = sprintf("dds%s_%02.0f_%s.tif", regionParam, tlo, calcType)
    # outputDir remains /tmp
} else {
    # we were given an output dir
    # does it contain path components?
    if ( grepl("/", outdir, perl=TRUE) ) { # NOT USING for WEB USE
        outputDir = sub("[^/]+$", "", outdir, perl=TRUE)
        outfilename = sprintf("dds%s_%02.0f_%s.tif", regionParam, tlo, calcType)
    } else { # THIS IS THE ONE FOR WEB USE by /dd/mapper program
		  base_dir  <- "/web/tmp2/"
        outputDir <- sprintf("%s%s%s",base_dir,outdir,"/")
        outfilename = sprintf("dds%s_%02.0f_%s.tif", regionParam, tlo, calcType)
    }
    #outfilename = sprintf("%s.tif", outfilename)
}

cat("DIRS INVOLVED:, outdir=", outdir,
    "    outputDir=", outputDir,
    "    outfilename=", outfilename, "\n")

#### Set up regions - use switch() (works like a single use hash) ####
#"MIDWEST"      = extent(-92,-90,30,49),
REGION = switch(regionParam,
"CONUS"        = extent(-125.0,-66.5,24.0,50.0),
"MIDWEST"      = extent(-104.2,-87,30,49),
"NORTHWEST"    = extent(-125.1,-103.8,40.6,49.2),
"SOUTHWEST"    = extent(-124.6,-101.5,31.2,42.2),
"SOUTHCENTRAL" = extent(-83.6,-78.3,31.8,35.3),
"NORTHCENTRAL" = extent(-104.3,-80.2,35.7,49.7),
"SOUTHEAST"    = extent(-107.1,-75.0,24.1,39.6),
"NORTHEAST"    = extent(-84.2,-64.3,36.9,48.1),
"AL"           = extent(-88.5294,-84.7506,30.1186,35.1911),
"AR"           = extent(-94.8878,-89.5094,32.8189,36.6936),
"AZ"           = extent(-115, -108.98, 31.2, 37),
"CA"           = extent(-124.6211, -113.7428, 32.2978, 42.2931),
"CO"           = extent(-109.2625, -101.8625, 36.7461, 41.2214),
"CT"           = extent(-73.7700, -71.7870, 40.9529, 42.0355),
"DL"           = extent(-76.1392, -74.1761, 38.3508, 39.9919),
"FL"           = extent(-87.8064, -79.9003, 24.4494, 31.1214),
"GA"           = extent(-85.7850, -80.5917, 30.1767, 35.1594),
"IA"           = extent(-96.8617, -89.9697, 40.1147, 43.7353),
"ID"           = extent(-117.3917, -110.6167, 41.4500, 49.3583),
"IL"           = extent(-91.5897, -87.0461, 36.8903, 42.6375),
"IN"           = extent(-88.1686, -84.4686, 37.7836, 41.9794),
"KS"           = extent(-102.3342, -94.1756, 36.6369, 40.2836),
"KT"           = extent(-89.3581, -81.8425, 36.4208, 39.3347),
"LA"           = extent(-94.3019, -88.7758, 28.8333, 33.2994),
"MA"           = extent(-73.5639, -69.7961, 41.1689, 42.9525),
"MD"           = extent(-79.7014, -74.8833, 37.0631, 39.9075),
"ME"           = extent(-71.4056, -66.6667, 42.9525, 47.5228),
"MI"           = extent(-90.5542, -82.3047, 41.6311, 47.5739),
"MN"           = extent(-97.4000, -89.3786, 43.2550, 49.3506),
"MO"           = extent(-95.8803, -88.9883, 35.8822, 40.7058),
"MS"           = extent(-91.7475, -87.8522, 29.9842, 35.2631),
"MT"           = extent(-116.3667, -103.8250, 44.0667, 49.3917),
"NC"           = extent(-84.440918, -75.300293, 33.682906, 36.646092),
"ND"           = extent(-104.2708, -96.3075, 45.6403, 49.1817),
"NE"           = extent(-104.3553, -95.0464, 39.7506, 43.2022),
"NH"           = extent(-72.6617, -70.6142, 42.6256, 45.4700),
"NJ"           = extent(-75.9175, -73.1892, 38.8944, 41.5806),
"NM"           = extent(-109.2942, -102.6383, 31.1892, 37.2000),
"NV"           = extent(-120.3358, -113.6803, 34.7356, 42.2981),
"NY"           = extent(-80.0867, -71.7381, 40.4828, 45.1692),
"OH"           = extent(-85.0439, -80.2464, 38.2797, 42.0217),
"OK"           = extent(-103.2850, -94.1964, 33.3839, 37.2850),
"OR"           = extent(-124.7294, -116.2949, 41.7150, 46.4612),
"PA"           = extent(-80.7672, -74.5033, 39.4694, 42.5094),
"SC"           = extent(-83.6422, -78.3275, 31.8814, 35.3811),
"SD"           = extent(-104.3553, -96.0806, 42.3050, 46.2050),
"TN"           = extent(-90.3239, -81.5047, 34.5578, 37.1125),
"TX"           = extent(-107.1592, -93.2411, 25.8614, 36.7200),
"UT"           = extent(-114.2925, -108.7450, 36.7778, 42.2347),
"VA"           = extent(-83.8322, -75.6200, 36.3892, 39.7886),
"VT"           = extent(-73.6747, -71.4108, 42.5886, 45.1956),
"WA"           = extent(-124.9585, -116.8364, 45.4554, 49.1664),
"WI"           = extent(-93.1572, -86.6822, 42.2733, 46.9914),
"WV"           = extent(-82.8783, -77.5114, 37.1158, 40.7836),
"WY"           = extent(-111.6167, -103.7333, 40.6667, 45.4833))


docalcs = switch(calcType,
                 "T1" = docalcsT1,
                 "T2" = docalcsT2,
                 "S1" = docalcsS1,
                 "S2" = docalcsS2,
                 "S1L" = docalcsS1L,
                 "S2L" = docalcsS2L,
                 "I1" = docalcsI1,
                 "I2" = docalcsI2,
                 "A"  = docalcsA,
                 "G"  = docalcsG,
                 "HD" = docalcsHD,
                 "CD" = docalcsCD,
                 "XD" = docalcsXD,
                 "Simple" = docalcsSimple
                 )

# unlike standard ddcalc, we're using the next day's min for the T2,
# S2, I2, not the previous day's max.  At least we'll try it.  It
# seems to me that it more accurately reflects the requested date
# interval, because the min happens is closer to midnight than the
# max.  But the difference is probably pretty small.
needNextMin = switch(calcType,
                     "T1" = FALSE,
                     "T2" = TRUE,
                     "S1" = FALSE,
                     "S2" = TRUE,
                     "S1L" = FALSE,
                     "S2L" = TRUE,
                     "I1" = FALSE,
                     "I2" = TRUE,
                     "A"  = FALSE,
                     "G"  = FALSE,
                     "HD" = FALSE,
                     "CD" = FALSE)
neededParams = c('tmax', 'tmin')
paramRegex = "tmax|tmin"

######## END of Param Handling Section #########

######## WEATHER INPUTS & OUTPUTS directory init #########
#for PRISM climate data w/subdirs 4-digit year

dates = seq(startDate, to=endDate, by='days')
years = unique(format(dates, "%Y"))
if (normalFlag) {
    pattern = paste("PRISM_(", paramRegex, ")_10yr.*_bil.bil$", sep="")
} else {
    pattern = paste("PRISM_(", paramRegex, ").*_bil.bil$", sep="")
}

baseDir = "/data/PRISM/symlinks/"
cat("BASE DIR: ",baseDir,"\n")
setwd(baseDir)

# create the map from dates to vectors of filenames and initialize all
# of the vectors to empty vectors.
# I assume there's an Rish way to do this, but envs won't take vector params
filesMap = new.env(hash=TRUE, parent=emptyenv(), size=2*length(dates))
for (datestr in format(dates, "%Y%m%d")) {
    filesMap[[datestr]] = c()
}

for (year in years){
    #yearDir = sprintf("%s%s",baseDir,year)
    #cat("WORKING DIR: ",yearDir,"\n")
    #setwd(yearDir)

    files = list.files(path=year,
                       pattern=pattern,
                       all.files=FALSE,
                       full.names=TRUE)

    for (file in files){
        datestr = strsplit(file,split="_")[[1]][5]
        filesMap[[datestr]] = append(filesMap[[datestr]], file)
    }
}

plus = function(a, b) a+b   # for use in .combine
final = foreach(date=dates, .combine=plus, .inorder=FALSE) %dopar% {
#final = foreach(date=dates, .combine=plus) %do% {
    datestr = format(date, "%Y%m%d")
    tmin = getBestRaster('tmin', datestr)
    tmax = getBestRaster('tmax', datestr)
    if (needNextMin) {
        #print("n6l6, ")
        tryCatch({
            nextDatestr = format(date+1, "%Y%m%d")
            nextmin = getBestRaster('tmin', nextDatestr)
        }, error=function(e) {
            # tomorrow's minimum not available; pretend today's is enough
            nextmin = min
        })
    } else
        nextmin = NA

    # We have all the layers we need; now do the calculations
    # Pass nextmin whether we need it or not (if not it's NA and gets ignored)
    dds = docalcs(tmin, tmax, nextmin, tlo, thi)
}

# keep us 'murkins happy
if (fahrenheit) final = 1.8*final

final
writeMap(final)
stop("END of DDCALC.R")
##

