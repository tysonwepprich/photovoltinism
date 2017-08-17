#!/usr/bin/Rscript
options(echo=FALSE)
## phen_model.R
# this version for ENTO adapt for CGI
args <- commandArgs(trailingOnly = TRUE)
 require(sp)
 require(rgdal)
 require(raster)
########           BEGINNING of Header Documentation Section                       #########
 #	              DDRP: Degree-Day, establishment Risk, and Pest event mapping system
 #	              By Len Coop, Gericke Cook, and Dan Upper for APHIS PPQ and IPM needs
 # alternative names to consider:
 #	  PestEventAndRiskMappingSystem PEARMS                                                
 #	  RiskDDAndPestEventMappingSystem RADDAPEMS                                                
 #	  DDRiskPOPS: 
 #	  LifeStageMaps
 #   LBPERMS R-based Lifestage-Based Phenology Events and Risks Mapping System
 #
 # Features in program provided by Gericke/PPQ (208 lines code):
 # 1. DD model plus exclusions (climatic limits) for combined risk maps predicting both potential
 #    number of generations and failure to establish
 # 2. Allow for multi-stage developmental and exclusion params (e.g. separate thresholds and limits
 #    for eggs larvae pupae adults)
 # 3. Step through daily PRISM 4km climate maps (Tmax,Tmin,Tmean,Precip) for calculations
 #
 # Features added/pending by Len/Dan:
 # 1. WORKING:v5: Additional DD calc formulas allowing better precision and use of upper thresholds (UDT)
 # 2. WORKING:v7: Addition of pest event map (PEM) creation
 # 3. WORKING v8: Option to switch off exclusion and PEM processing
 # 4. WORKING:v15: Addition of stage overlap factor e.g. eggs do not all hatch to larvae at same time; with
 #    use of an overlap parameter and ramp functions
 # 5. WORKING:v16: Define a subregion, set up numerous subregions for various needs
 # 6. WORKING:v18: Add OW stage initialization and flexible stage ordering e.g.1: OW adult->e,l,p,a,repeat
 #    e.g.2: OW egg->l,p,a,e,repeat
 # 7. WORKING: v20: 1) Simple parameter passing suitable for command line and cgi-bin scripts (658 lines code):
 #                  2) Stage overlap param (0.0 to 0.50) now working 
 #      Params passed in: 
 #               SPP (3letterabbrev) 
 #               YEAR (4 digit) 
 #               START_DOY (dayofyear) 
 #               REGION (CONUS,OR,MIDWEST currently)
 #               STAGE_OVERLAP (0-0.5) 
 #               EXCLUSIONS (Y/N)
 #               PEMS (Y/N)
 #      ex: ./DDRiskPEMSv21.R FCM 2010 39 OR 0.35 N N
 #   showing correct DD calcs, NumGens, PEMs for 3 stages and 2 generations as examples
 # 8. WORKING v22: 1) Read pest/insect params from external library/database file (736 lines code)
 #                 2) Add data output dir different from data input dir
 #            ex: ./DDRiskPEMSv22.R FCM 2010 39 OR 0.35 N N results3
 # 9. WORKING: v23 file handling (only need tmax, tmin, ppt): 
 #              1) a) past years - use approp year daily 4km PRISM data (stable) - yes 
 #                    (800m daily PRISM is available to us on word of PRISM group but our servers probably not
 #                    ready to handle such high resolution, hopefully new server will be)
 #					   b) for current year (2015) use PRISM in order of avail: 
 #                   b1) stable b2) provisional b3) early  
 #                c) forecast (NDFD 2.5km is avail to try), need to work on naming conventions and resolution
 #		            d) extended forecast starting w/10-year averages (800m avail), need naming conventions and resolution
 #                e)   by 2016 should have 90+ day forecasts via ARDP grant work
 #                f) FUTURE PRISM data: possibly Dewpoint, RHmean/VPD or other moisture params
 #             3) Control over outputs:
 #                a) PEMS should be automatically on if PEMS is on
 #                b) More control over stages and times of PEMS w/species param files - working
 #                c) need params passed in for freq of outputs: daily, weekly, 10 days, 14 days, monthly, just the
 #                    final maps.
 #                    normally would set a last date and only output a map on that date
 #                    also a way to specify stage to focus on?
 #                    try starting with categories:
 #                    0=no extra output (defaults in spp file?)
 #                    1=outputs last day of month
 #                d) Control over which stages to map
 #                e) Consider degrees F on input (harder?) as well as on output (easier)
 #             4) Web/CGI-BIN interface         
 #                 a) all command line params available as form elements
 #                 b) have a download link for all files
 #                 c) make .png files for viewing that includes legend, title, statelines
 #                 d) use PID or random number (web) vs user specified output directory
 #                 e) full set of CONUS states and regions
 #                 f) have full metadata files available
 # 10. IN DEVEL: v24 expand exclusions to allow:
 #                 a) accumulate chill units 
 #                 b) separate tracking for end-of-year masking in addition to realtime exclusions as in previous versions
########                END of Header Documentation Section                        #########

########                BEGINNING of Param Handling Section                        #########
#### Default values for command line params ####
spp           <- "FCM"       # default species to use
start_year    <- 2011        # only one 4 digit year currently supported 
start_doy     <- 39          # day of year to use for biofix             
end_doy       <- 365         # day of year to use for stop date - need 365 if voltinism map 
region_param  <- "MIDWEST"   # default REGION to use
ovlp          <- 0.30        # proportion overlap stage to stage; must be < 0.5
exclusions    <- 0           # ability to turn on/off climate limits/exclusions 1=on 2=also apply exclusions in realtime 3 = stress unit exclusions
pems          <- 0           # turn on/off pest event maps
out_dir       <- "results"   # output directory (currently appended to /data/PRISM)
out_option    <- 1           # output option category (**to be defined**)
mapA          <- 1           # make maps for adult stage
mapE          <- 0           # make maps for egg stage
mapL          <- 0           # make maps for larval stage
mapP          <- 0           # make maps for pupal stage
do_ovlp       <- 0           # added 6/1/16 to only calc ovlp if we make maps on that day
do_maps       <- 0           # added 6/1/16 to only make maps on days selected
                             # could add: ovlp differ by stage, incr each generation
#### Process command line args ####
if (length(args) < 12) {
   cat("  ","\n")
   cat("DDRP: Degree-Day, establishment Risk, and Pest Event Mapping System","\n")
   cat("Required params: SPP (3letterabbrev)       ","\n") 
   cat("                 YEAR (4 digit)            ","\n") 
   cat("                 START_DOY (dayofyear)     ","\n") 
   cat("                 END_DOY (dayofyear)       ","\n") 
   cat("                 REGION (CONUS,regions,48 states currently)","\n")
   cat("                 OVERLAP (0-0.5)           ","\n") 
   cat("                 EXCLUSIONS (0=off,1=on,2=apply in realtime,3=stressunit-based exclusions)","\n")
   cat("                 Pest Event Maps (Y/N)     ","\n")
   cat("                 Name of output dir (results)","\n")
   cat("                 Output option 1-8 (1=default)","\n")
   cat("                 Make Maps for Adults (1=yes)","\n")
   cat("                 Make Maps for Eggs   (1=yes)","\n")
   cat("                 Make Maps for Larvae (1=yes)","\n")
   cat("                 Make Maps for Pupae  (1=yes)","\n")
   cat("Run example:     DDRPv24.R FCM 2010 39 365 MIDWEST 0.25 0 N results 1 1 0 0 0","\n")
   cat("  ","\n")
   q()
} else {
	cat("PARAM LIST:  ")
   print(args)
}

#### Read in command line args ####
spp           <- args[1] 
start_year    <- as.integer(args[2])
start_doy     <- as.integer(args[3])
end_doy       <- as.integer(args[4])
region_param  <- args[5]
ovlp          <- as.double(args[6])
exclusions    <- args[7]
pems          <- args[8]
out_dir       <- args[9]
out_option    <- args[10]
mapA          <- args[11]
mapE          <- args[12]
mapL          <- args[13]
mapP          <- args[14]

#### Use gsub to make secure cgi-bin args : equiv. of perl spp =~ s/[^A-Za-z0-9]//g; ####
spp           <- gsub('[^A-Za-z0-9]',"",spp,perl=TRUE)
start_year    <- gsub('[^0-9]',"",start_year,perl=TRUE)
start_doy     <- gsub('[^0-9.]',"",start_doy,perl=TRUE)
end_doy       <- gsub('[^0-9.]',"",end_doy,perl=TRUE)
region_param  <- gsub('[^A-Za-z0-9."]',"",region_param,perl=TRUE)
ovlp          <- gsub('[^0-9.]',"",ovlp,perl=TRUE)
#exclusions    <- gsub('[^A-Z]',"",exclusions,perl=TRUE)
exclusions    <- gsub('[^0-9]',"",exclusions,perl=TRUE)
pems          <- gsub('[^A-Z]',"",pems,perl=TRUE)
out_dir       <- gsub('[^A-Za-z0-9./]',"",out_dir,perl=TRUE)
out_option    <- gsub('[^0-9.]',"",out_option,perl=TRUE)
mapA          <- gsub('[^0-9.]',"",mapA,perl=TRUE)
mapE          <- gsub('[^0-9.]',"",mapE,perl=TRUE)
mapL          <- gsub('[^0-9.]',"",mapL,perl=TRUE)
mapP          <- gsub('[^0-9.]',"",mapP,perl=TRUE)

#### Convert var types for some vars ####
start_year    <- as.integer(start_year)
start_doy     <- as.integer(start_doy)
end_doy       <- as.integer(end_doy)
ovlp          <- as.double(ovlp)
out_option    <- as.integer(out_option)
mapA          <- as.integer(mapA)
mapE          <- as.integer(mapE)
mapL          <- as.integer(mapL)
mapP          <- as.integer(mapP)


#### Set up regions - use switch() (works like a single use hash) ####
#"MIDWEST"      = extent(-92,-90,30,49),
REGION <- switch(region_param,
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
	
#if (exclusions == "Y") { 
#	exclusions <- 1 
#} else {
#	exclusions <- 0 
#}

if (exclusions == "2") { # OK both use climate limiting exclusions and apply realtime masks affecting all output 
	exclusions <- 1 
	exclusionsrealtime <- 1 
   exclusions_stressunits <- 0
} else if (exclusions == "1") { # use exclusions but dont affect phenology/other maps until we apply them at end of model run 
	exclusions <- 1 
	exclusionsrealtime <- 0 
   exclusions_stressunits <- 0
} else if (exclusions == "3") {  #  == "3" NEW type stress unit based exclusions
	exclusions <- 0 
	exclusionsrealtime <- 0 
   exclusions_stressunits <- 1
} else {  #  == "0" dont use any exclusions
	exclusions <- 0 
	exclusionsrealtime <- 0 
   exclusions_stressunits <- 0
}


if (pems == "Y") { 
	pems <- 1 
} else { 
	pems <- 0 
}

#### Init. mapping output freq. options ####
finalmaps    <- 0     # minimal output only when last day of model run is done
monthlymaps  <- 0
biweeklymaps <- 0
dekadmaps    <- 0    # a dekad is every 10 days FYI e.g. Aug 1, 11, 21, 31
weeklymaps   <- 0
evendaymaps  <- 0
dailymaps    <- 0

if (out_option == "0") { 
  finalmaps = 1
} else if (out_option == "1") { 
  monthlymaps  <- 1
  finalmaps    <- 1
} else if (out_option == "2") { 
  biweeklymaps <- 1
  finalmaps    <- 1
} else if (out_option == "3") { 
  dekadmaps    <- 1
  finalmaps    <- 1
} else if (out_option == "4") { 
  weeklymaps   <- 1
  finalmaps    <- 1
} else if (out_option == "5") { 
  evendaymaps  <- 1
  finalmaps    <- 1
} else if (out_option == "6") { 
  dailymaps    <- 1
  finalmaps    <- 1
}
cat("CLEAN PARAMS: spp startyr startdoy enddoy region_param ovlp exclusions pems out_dir output_option mapA mapE mapL mapP: \n              ",spp,start_year,start_doy,end_doy,region_param,ovlp,exclusions,pems,out_dir,out_option,mapA,mapE,mapL,mapP,"\n")
########                      END of Param Handling Section                        #########

########                   BEGINNING of Directory Init Section                     #########
#### PARAM INPUTS - this is location for species params; thresholds etc. ####
params_dir <- "/usr/local/dds/DDRP/spp_params/"

#### FOR WEATHER INPUTS & OUTPUTS - PRISM climate data w/subdirs 4-digit year ####
if (grepl("16",start_year,perl=TRUE)) {  # if outdir has 2 consec. numbers, assume webuser
  base_dir <- "/mnt/ssd1/PRISM/"
} else {  # otherwise just use base dir
  base_dir <- "/data/PRISM/"
}
cat("BASE DIR: ",base_dir,"\n")
prism_dir <- sprintf("%s%s",base_dir,start_year)
cat("WORKING DIR: ",prism_dir,"\n")

# MAP OUTPUTS
if (grepl("[0-9]{3,5}",out_dir,perl=TRUE)) {  # if outdir has 3 to 5 consec. numbers, assume webuser
  sink("/tmp/rlogging.txt")
  output_dir <- sprintf("%s%s","/home/httpd/html/tmp/",out_dir)
} else {  # otherwise just use base dir
  output_dir <- sprintf("%s%s",base_dir,out_dir)
}

if (file.exists(output_dir)) {
  cat("EXISTING OUTPUT DIR: ",output_dir,"\n")
} else {
  dir.create(output_dir)
  cat("NEW OUTPUT DIR: ",output_dir,"\n")
}
setwd(prism_dir)
########                        END of Directory Init Section                      #########

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


#### Write out maps as .tif GIS data files, use GRASS to produce PNG copies for viewing as maps
ConvPNG=function(tif_file){
        #system2('/usr/local/dds/DDRP/dograss_NW_NB_54',args=c("/usr/local/dds/DDRP/dograss_png.sh",tif_file,spp),stdout="",stderr="",stdin="")
        cat("2spp fullname: ",spp,fullname,region_param,"\n")
        system2('/usr/local/dds/DDRP/dograss_NW_NB_54',args=c("/usr/local/dds/DDRP/dograss_png.sh",tif_file,spp,fullname,region_param),stdout="",stderr="",stdin="")
}
cat("6WORKING DIR: ",prism_dir,"\n")

WriteMaps=function(d){   
			  setwd(output_dir)
			  # always show lifestages, NumGen, StageCount
			  if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays main maps
             rast_name <- sprintf("%s%s%s","All_Stg_Excl_",d,".tif")
             writeRaster(allminEXCL,paste("All_Stg_Excl_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
				 Lifestage <- Cond(allminEXCL,Lifestage,0)
             rast_name <- sprintf("%s%s%s","Lifestage_",d,".tif")
             writeRaster(Lifestage,paste("Lifestage_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
				 NumGen <- Cond(allminEXCL,NumGen,0)
             rast_name <- sprintf("%s%s%s","NumGen_",d,".tif")
             writeRaster(NumGen,paste("NumGen_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
				 StageCount <- Cond(allminEXCL,StageCount,0)
             rast_name <- sprintf("%s%s%s","StageCount_",d,".tif")
             writeRaster(StageCount,paste("StageCount_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
				 # dont apply to ddtotal
             rast_name <- sprintf("%s%s%s","ddtotal_",d,".tif")
             writeRaster(ddtotal,paste("ddtotal_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
			  }
			  else { # write regular old maps
             rast_name <- sprintf("%s%s%s","Lifestage_",d,".tif")
             writeRaster(Lifestage,paste("Lifestage_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             rast_name <- sprintf("%s%s%s","NumGen_",d,".tif")
             writeRaster(NumGen,paste("NumGen_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             rast_name <- sprintf("%s%s%s","StageCount_",d,".tif")
             writeRaster(StageCount,paste("StageCount_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             rast_name <- sprintf("%s%s%s","ddtotal_",d,".tif")
             writeRaster(ddtotal,paste("ddtotal_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
			  }
			  if (exclusions_stressunits == 1) { # make additional maps for "stageless" stress maps
             rast_name <- sprintf("%s%s%s","Chill_Stress_Units_",d,".tif")
             writeRaster(chillunitsCUM,paste("Chill_Stress_Units_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             rast_name <- sprintf("%s%s%s","Chill_Stress_Excl_",d,".tif")
             writeRaster(chillEXCL,paste("Chill_Stress_Excl_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             rast_name <- sprintf("%s%s%s","Heat_Stress_Units_",d,".tif")
             writeRaster(heatunitsCUM,paste("Heat_Stress_Units_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             rast_name <- sprintf("%s%s%s","Heat_Stress_Excl_",d,".tif")
             writeRaster(heatEXCL,paste("Heat_Stress_Excl_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             rast_name <- sprintf("%s%s%s","All_Stress_Excl_",d,".tif")
             writeRaster(AllEXCL,paste("All_Stress_Excl_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
				 # Exclusions included: -2=severe, -1=moderate stress
				 #Lifestage here:
             LifestageEXCL1 <- Cond((AllEXCL > -2),Lifestage,-2)
             rast_name <- sprintf("%s%s%s","LifestageExcl1_",d,".tif")
             writeRaster(LifestageEXCL1,paste("LifestageExcl1_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             LifestageEXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,Lifestage))
             rast_name <- sprintf("%s%s%s","LifestageExcl2_",d,".tif")
             writeRaster(LifestageEXCL2,paste("LifestageExcl2_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
				 #NumGen here:
             NumGenEXCL1 <- Cond((AllEXCL > -2),NumGen,-2)
             rast_name <- sprintf("%s%s%s","NumGenExcl1_",d,".tif")
             writeRaster(NumGenEXCL1,paste("NumGenExcl1_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             NumGenEXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,NumGen))
             rast_name <- sprintf("%s%s%s","NumGenExcl2_",d,".tif")
             writeRaster(NumGenEXCL2,paste("NumGenExcl2_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
				 #StageCount here:
             StageCountEXCL1 <- Cond((AllEXCL > -2),StageCount,-2)
             rast_name <- sprintf("%s%s%s","StageCountExcl1_",d,".tif")
             writeRaster(StageCountEXCL1,paste("StageCountExcl1_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             StageCountEXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,StageCount))
             rast_name <- sprintf("%s%s%s","StageCountExcl2_",d,".tif")
             writeRaster(StageCountEXCL2,paste("StageCountExcl2_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
			  }
             # rast_label not in use delete all except this example:
             ##rast_label <- sprintf("%s%s%s",spp," Percent Egg Devel. on ",d)
           if (mapE == 1) { #show egg stage
             #rast_name <- sprintf("%s%s%s","EggDev_",d,".tif")
             #writeRaster(eggDev,paste("EggDev_",d,sep=""), format="GTiff", overwrite=TRUE)
             #ConvPNG(rast_name)
			    if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays DDs
               rast_name <- sprintf("%s%s%s","EggP_Exclusions_",d,".tif")
					eggP <- Cond(allminEXCL,eggP,0)
               writeRaster(eggP,paste("EggP_Exclusions",d,sep=""), format="GTiff", overwrite=TRUE)
               ConvPNG(rast_name)
				 }
				 else {    # no exclusions used this run
               rast_name <- sprintf("%s%s%s","EggP_",d,".tif")
               writeRaster(eggP,paste("EggP_",d,sep=""), format="GTiff", overwrite=TRUE)
               ConvPNG(rast_name)
				 }
				 if (exclusions == 1) {   # display eggmindays map
               rast_name <- sprintf("%s%s%s","Egg_chill_days",d,".tif")
               writeRaster(eggmindays,paste("Egg_chill_days",d,sep=""), format="GTiff", overwrite=TRUE)
               ConvPNG(rast_name)
               rast_name <- sprintf("%s%s%s","Egg_exclusion_areas",d,".tif")
               writeRaster(eggminEXCL,paste("Egg_exclusion_areas",d,sep=""), format="GTiff", overwrite=TRUE)
               ConvPNG(rast_name)
				 }
			  }
			  if (mapL == 1) {
			    if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays DDs
               rast_name <- sprintf("%s%s%s","LarvP_Exclusions_",d,".tif")
					larvP <- Cond(allminEXCL,larvP,0)
               writeRaster(larvP,paste("LarvP_Exclusions",d,sep=""), format="GTiff", overwrite=TRUE)
               ConvPNG(rast_name)
				 }
				 else {    # no exclusions used this run
               rast_name <- sprintf("%s%s%s","LarvP_",d,".tif")
               writeRaster(larvP,paste("LarvP_",d,sep=""), format="GTiff", overwrite=TRUE)
               ConvPNG(rast_name)
				 }
				 if (exclusions == 1) {   # display larvaemindays map
               rast_name <- sprintf("%s%s%s","Larv_chill_days",d,".tif")
               writeRaster(larvaemindays,paste("Larv_chill_days",d,sep=""), format="GTiff", overwrite=TRUE)
               ConvPNG(rast_name)
               rast_name <- sprintf("%s%s%s","Larv_exclusion_areas",d,".tif")
               writeRaster(larvaeminEXCL,paste("Larv_exclusion_areas",d,sep=""), format="GTiff", overwrite=TRUE)
               ConvPNG(rast_name)
				 }
			  }
# FIXME - Pupal exclusions should be fully fleshed out
			  if (mapP == 1) {
             rast_name <- sprintf("%s%s%s","PupP_",d,".tif")
             writeRaster(pupP,paste("PupP_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
			  }
			  if (mapA == 1) {
			    if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays DDs
               rast_name <- sprintf("%s%s%s","AdultP_",d,".tif")
					adultP <- Cond(allminEXCL,adultP,0)
               writeRaster(adultP,paste("AdultP_",d,sep=""), format="GTiff", overwrite=TRUE)
               ConvPNG(rast_name)
				 }
				 else {    # no exclusions used this run
               rast_name <- sprintf("%s%s%s","AdultP_",d,".tif")
               writeRaster(adultP,paste("AdultP_",d,sep=""), format="GTiff", overwrite=TRUE)
               ConvPNG(rast_name)
				 }
				 if (exclusions == 1) {   # display adultmindays map
               rast_name <- sprintf("%s%s%s","Adult_chill_days",d,".tif")
               writeRaster(adultmindays,paste("Adult_chill_days",d,sep=""), format="GTiff", overwrite=TRUE)
               ConvPNG(rast_name)
               rast_name <- sprintf("%s%s%s","Adult_exclusion_areas",d,".tif")
               writeRaster(adultminEXCL,paste("Adult_exclusion_areas",d,sep=""), format="GTiff", overwrite=TRUE)
               ConvPNG(rast_name)
				 }
           }
			  setwd(prism_dir)
}
########                        END function definitions                           #########

########                   BEGINNING of Initialization Section                     #########
#### Abbrevs used in model params:
  #LDT = lower development threshold, temp at which growth = 0 (using PRISM tmean)
  #DD = degree days, number of cumulative heat units to complete that lifestage
  #OW = overwintering (stage, DD reqs., etc)

#### Pest Specific, Multiple Life Stage Phenology Model Parameters: ####
## now read from source param files in ./spp_params/SPP.params
  param_file <- sprintf("%s%s",spp,".params")
  species_params <- sprintf("%s%s",params_dir,param_file)
  if (file.exists(species_params)) {
    cat("SPP PARAMS: ",species_params,"\n")
    source(species_params)
    cat("Reading Params for SPP: ",spp," Fullname: ",fullname,"\n")
  } else {
    cat("PARAM FILE: ",species_params,"...Not found; exiting Program...\n")
    q()  # no reason to keep going without any params
  }

#### Push out a metadata file with all inputs used in model   ####
    cat("output dir: ",output_dir,"\n")
cat("made it here line xxx \n")
setwd(output_dir)
metadata <- sprintf("%s%s","./","metadata.txt")
if (file.exists(metadata)) {
  cat("Metadata for DDRP - Degree-Day, Risk, and Pest event Maps \n",file=metadata)
} else {
  file.create(metadata)
  cat("Metadata for DDRP - Degree-Day, Risk, and Pest event Maps \n",file=metadata)
}
########                      END of Param Handling Section                        #########


########                      BEGIN Start Metadata Output File                     #########

#cat("CLEAN PARAMS: spp startyr startdoy enddoy region_param ovlp exclusions pems out_dir output_option mapA mapE mapL mapP: \n            ",spp,start_year,start_doy,end_doy,region_param,ovlp,exclusions,pems,out_dir,out_option,mapA,mapE,mapL,mapP,"\n")

cat("\n Model Species Params:\n Species Abbrev:",spp,"\n Full Name:",fullname,"\n Pest of:",pestof,file=metadata,append=TRUE)
cat("\n Overwintering Stage:",owstage,"\n Egg Lower Devel Threshold:",eggLDT,"\n Egg Upper Devel Threshold:",eggUDT,file=metadata,append=TRUE)
cat("\n Larvae Lower Devel Threshold:",larvaeLDT,"\n Larvae Upper Devel Threshold:",larvaeUDT,"\n Pupae Lower Devel Threshold:",pupaeLDT,file=metadata,append=TRUE)
cat("\n Pupae Upper Devel Threshold:",pupaeUDT,"\n Adult Lower Devel Threshold:",adultLDT,"\n Adult Upper Devel Threshold:",adultUDT,file=metadata,append=TRUE)

if (exclusions_stressunits) {
   if (chillstress_threshold) {
   cat("\n Lower Chill Threshold:",chillstress_threshold,"\n Upper Heat Threshold:",heatstress_threshold,file=metadata,append=TRUE)
   cat("\n Max Chill Units (lower bound):",chillstress_units_max1,"\n Max Chill Units (upper bound):",chillstress_units_max2,file=metadata,append=TRUE)
   cat("\n Max Heat Stress Units (lower bound):",heatstress_units_max1,"\n Max Heat Stress Units (upper bound):",heatstress_units_max2,file=metadata,append=TRUE)
   }
}

if (exclusions) {
   if (eggLLT) {
   cat("\n Egg Lower Lethal Threshold:",eggLLT,"\n Egg Upper Lethal Threshold:",eggULT,file=metadata,append=TRUE)
   }
   if (eggLLDAYS)    {
   cat("\n Egg Lower Lethal Days Needed:",eggLLDAYS,file=metadata,append=TRUE)
   }
   cat("\n Larvae Lower Lethal Threshold:",larvaeLLT,"\n Larvae Upper Lethal Threshold:",larvaeULT,file=metadata,append=TRUE)
   if (larvaeLLDAYS) {
   cat("\n Larvae Lower Lethal Days Needed:",larvaeLLDAYS,file=metadata,append=TRUE)
   }
   cat("\n Adult Lower Lethal Threshold:",adultLLT,file=metadata,append=TRUE)
   if (adultLLDAYS) {
   cat("\n Adult Lower Lethal Days Needed:",adultLLDAYS,file=metadata,append=TRUE)
   }
}

if (pems) {
cat("\n Number of generations to make Pest Event Maps (PEMs):",PEMnumgens,file=metadata,append=TRUE)
cat("\n Egg Event DDs and Label:",eggEventDD," ",eggEventLabel,file=metadata,append=TRUE)
cat("\n Larvae Event DDs and Label:",larvaeEventDD," ",larvaeEventLabel,file=metadata,append=TRUE)
cat("\n Pupae Event DDs and Label:",pupaeEventDD," ",pupaeEventLabel,file=metadata,append=TRUE)
cat("\n Adult Event DDs and Label:",adultEventDD," ",adultEventLabel,file=metadata,append=TRUE)
#  adultEventDD <- 30 # PEMs for adult stage (1st oviposition) is 22 DDs into adult stage
#  adultEventLabel <- "Approx. 1st egglaying by females" # Label for PEM adult stage
}

cat("\n\n Model Input Params:\n Start Year:",start_year,"\n Start Day-of-year:",start_doy,file=metadata,append=TRUE)
cat("\n End Day-of-Year:",end_doy,"\n Region:",region_param,"\n Stage overlap:",ovlp,file=metadata,append=TRUE)
cat("\n Exclusion Maps:",exclusions,"\n PestEvent Maps:",pems,"\n Output_Dir:",out_dir,file=metadata,append=TRUE)
cat("\n Output Option:",out_option,"\n Map Adults:",mapA,"\n Map Eggs:",mapE,file=metadata,append=TRUE)
cat("\n Map Larvae:",mapL,"\n Map Pupae:",mapP,file=metadata,append=TRUE)

#spp startyr startdoy enddoy region_param ovlp exclusions pems out_dir output_option mapA mapE mapL mapP: \n            ",spp,start_year,start_doy,end_doy,region_param,ovlp,exclusions,pems,out_dir,out_option,mapA,mapE,mapL,mapP,"\n")

#cat(species_params,file=metadata,append=TRUE)

setwd(prism_dir)
########                        END Start Metadata Output File                     #########


#### Example params spp=FCM in case you dont have any param files: ####
  #fullname   <- "False Codling Moth"
  #stgorder   <- c("OA","E","L","P","A","F")
  #owstage    <- "OA"
  #eggLDT     <- 11.93
  #eggUDT     <- 40  # Daiber: 40C as upper dev. threshold
  #larvaeLDT  <- 11.6
  #larvaeUDT  <- 34 #upper dev. threshold-need to verify 
  #pupaeLDT   <- 11.9
  #pupaeUDT   <- 40
  #adultLDT   <- 12.2 #for oviposition
  #adultUDT   <- 40
  #eggDD      <- 69  # round from 69.3
  #larvaeDD   <- 156
  #pupDD      <- 174 #females
  #OWpupDD    <- 86  # text OW stage dev 39 DD "post diapause"
  #adultDD    <- 79 # round from 79.2 time to 50% eggs laid
  #OWadultDD  <- 86  # text OW stage dev 39 DD "post diapause"
  #calctype   <-"average"
  #LLT = lower lethal temperature (PRISM tmin), ULT = upper lethal temperature (PRISM tmax)
  #eggLLT     <- -3
  #eggLLDAYS  <- 7  # NEW v24 add # days to accumulate for local extinction
  #eggULT     <- 41
  #larvaeLLT  <- -12
  #larvaeLLDAYS <- 2
  #larvaeULT  <- 40
  #adultLLT   <- 0.5
  #adultLLDAYS <- 5
  # Pest Event Maps (PEMs) must be turned on for these to get used:
#if (pems) {  # init as spp_param files may not have these (but should if PEMS turned on)
#  PEMnumgens <- 2  # create PEMS for up to this many generations (max is 4)
#  eggEventDD <- 5 # PEMs for egg stage is 5 DDs into egg stage
#  eggEventLabel <- "Beginning of egg hatch" # Label for PEM egg stage
#  larvaeEventDD <- 10 # PEMs for larvae stage is 78 DDs into larval stage
#  larvaeEventLabel <- "Early larval development" # Label for PEM larval stage
#  pupaeEventDD <- 10 # PEMs for pupal stage is 58 DDs into pupal stage
#  pupaeEventLabel <- "Early pupal development" # Label for PEM pupal stage
#  adultEventDD <- 30 # PEMs for adult stage (1st oviposition) is 22 DDs into adult stage
#  adultEventLabel <- "Approx. 1st egglaying by females" # Label for PEM adult stage
#}
#### END Pest Specific, Multiple Life Stage Phenology Model Parameters ####

#### Search pattern for PRISM/other daily temperature grids. Load them for processing. ####
pattern = paste("(PRISM_tmax_)(.*)(_bil.bil)$", sep="")
#pattern = paste("(*_tmax_)(.*)(_bil.bil)$", sep="")
files <- list.files(pattern=pattern, all.files=FALSE, full.names=TRUE)
#Check that there are enough files for the year
#length(files)
numlist <- vector()

#### New code (V23) to add to DDRPx to pick most useful/recent file types needed for current year data ####
# stable (4)  > provisional (3) > early (2)  > 7day forecast (NA) > 90day fc (NA) > 10yr Avg data (new!) (1) > none (0)
# just do for tmax, assume tmin is same (add error flagging somehow??)
filelist <- vector()
# set up R version of hash tables (new.env)
maxfiles <- new.env(hash=TRUE, parent=emptyenv(), size=100L)
minfiles <- new.env(hash=TRUE, parent=emptyenv(), size=100L)
quality  <- new.env(hash=TRUE, parent=emptyenv(), size=100L)

# init file type hash keys
for (mon in 1:12) {
   for (day in 1:31) {
      skey = sprintf("%04s%02d%02d",start_year,mon,day)
      quality[[skey]] <- 0
   }
}

fkey <- 0
########                      END of Initialization Section                        #########


########                 Begin of actual stepping through the model                #########
#### Go through all files, store best set in a hash using these 0 to 4 priority numbers ####
for (file in sort(files)) {
    status = strsplit(file,split="_")[[1]][3]  # type of data: early provisional stable etc
    num = strsplit(file,split="_")[[1]][5]     # date e.g. 20150917
    fkey <- num

    if ((status == "stable") && (quality[[fkey]] < 5)) {
      filelist <-c(filelist,file)
      numlist <-c(numlist,num)
      quality[[fkey]] <- 5
      maxfiles[[fkey]] <- file  # tmax stored now check tmin
      minfile  <- gsub("tmax","tmin",file,perl=TRUE)
      if (file.exists(minfile)) {
        minfiles[[fkey]] <- minfile  # tmin exists assume its good
		}
		else {
        cat("no tmin file found equiv to avail tmax file: ",minfile,"\n")  # FIXME no smarts to use lower grade tmin
		}
    } else if ((status == "provisional") && (quality[[fkey]] < 4)) {
      filelist <-c(filelist,file)
      numlist <-c(numlist,num)
      quality[[fkey]] <- 4
      maxfiles[[fkey]] <- file
      minfile  <- gsub("tmax","tmin",file,perl=TRUE)
      if (file.exists(minfile)) {
        minfiles[[fkey]] <- minfile
		}
		else {
        cat("no tmin file found equiv to avail tmax file: ",minfile,"\n")
		}
    } else if ((status == "early") && (quality[[fkey]] < 3)) {
      filelist <-c(filelist,file)
      numlist <-c(numlist,num)
      quality[[fkey]] <- 3
      maxfiles[[fkey]] <- file
      minfile  <- gsub("tmax","tmin",file,perl=TRUE)
      if (file.exists(minfile)) {
        minfiles[[fkey]] <- minfile
		}
		else {
        cat("no tmin file found equiv to avail tmax file: ",minfile,"\n")
		}
	#
	# future short and extended forecast options/filetypes get inserted here
	#
    #} else if ((status == "10yr0716") && (quality[[fkey]] < 2)) {  # 10 yr avg 2006-2015 PRISM data
    } else if (status == "10yr0716")  {  # 10 yr avg 2006-2015 PRISM data
      filelist <-c(filelist,file)
      numlist <-c(numlist,num)
      quality[[fkey]] <- 2
      maxfiles[[fkey]] <- file
      minfile  <- gsub("tmax","tmin",file,perl=TRUE)
      if (file.exists(minfile)) {
        cat("found 10yr minfile: ",minfile,"\n")
        minfiles[[fkey]] <- minfile
		}
		else {
        cat("no tmin file found equiv to avail tmax file: ",minfile,"\n")
		}
   } else if ((status == "10yr0615") && (quality[[fkey]] < 1)) {
      filelist <-c(filelist,file)
      numlist <-c(numlist,num)
      quality[[fkey]] <- 1
      maxfiles[[fkey]] <- file
      minfile  <- gsub("tmax","tmin",file,perl=TRUE)
      if (file.exists(minfile)) {
        #cat("found nmme minfile: ",minfile,"\n")
        minfiles[[fkey]] <- minfile
      }
      else {
        cat("no tmin file found equiv to avail tmax file: ",minfile,"\n")
      }

    } else {
        cat("no match for this file: ",file," ",num,"\n")
    }
}
# cat("Done1. Full numlist: ",sort(numlist),"\n")
sortedlist <- sort(numlist)
#### END code to distinguish file type/quality  ####

#### Read in first raster as a template, crop to REGIONS setting and intialize tracking rasters ####
# SET Region Here:
  #  cat("LINE 313:",sortedlist,"\n")
template <- crop(raster(files[1]),REGION)
#template <- crop(raster(maxfiles[[1]]),REGION)
template[!is.na(template)] <- 0

#### Initialize all tracking rasters as zero with the template ####
ddtotal <- template
DDaccum <- template
Lifestage <- template
Lifestage <- Lifestage -1  # Need for OW designated as stage -1
NumGen <- template
StageCount <- template

if (exclusions) {
 eggminEXCL        <- template  # final map with regions with mindays >= mindaysMAX
 eggmindays        <- template  # track no. days Tlow < eggLLT
 eggmindaysMAX     <- template
 eggmindaysMAX     <- eggLLDAYS # max no. days Tlow < eggLLT = exclusion
 larvaeminEXCL     <- template  # final map with regions with mindays >= mindaysMAX
 larvaemindays     <- template
 larvaemindaysMAX  <- template
 larvaemindaysMAX  <- larvaeLLDAYS # max no. days Tlow < larvaeLLT = exclusion
 #pupaeminEXCL      <- template  # final map with regions with mindays >= mindaysMAX
 #pupaemindays      <- template
 #pupaemindaysMAX   <- template
 #pupaemindaysMAX   <- pupaeLLDAYS # max no. days Tlow < pupaeLLT = exclusion
 adultminEXCL      <- template  # final map with regions with mindays >= mindaysMAX
 adultmindays      <- template
 adultmindaysMAX   <- template
 adultmindaysMAX   <- adultLLDAYS # max no. days Tlow < adultLLT = exclusion
 allminEXCL        <- template  # product of all EXCL maps
} else if (exclusions_stressunits) {
	# NEW approach chill/heat stress units
 chillmask         <- template  # binary mask for daily chill units
 chillstress       <- template  # count of daily chill units
 chillstressTHRESH  <- template  # mask for chillstress units threshold
 chillstressTHRESH  <- chillstress_threshold # mask for chillstress units threshold
 chillunitsCUM     <- template  # cumulative chill units
 chillstressMAX1    <- template  # use for max chill before mortality??
 chillstressMAX1    <- chillstress_units_max1 # use for max chill before mortality??
 chillstressMAX2    <- template  # use for max chill before mortality??
 chillstressMAX2    <- chillstress_units_max2 # use for uncertainty zone max chill before mortality??
 chillEXCL         <- template  # EXCL map for chilling
 heatmask          <- template  # binary mask for daily heat stress units
 heatstress        <- template  # use for heat stress units
 heatstressTHRESH  <- template  # mask for heatstress units threshold
 heatstressTHRESH  <- heatstress_threshold # mask for heatstress units threshold
 heatunitsCUM      <- template  # cumulative heat stress units
 heatstressMAX1    <- template  # use for max chill before mortality??
 heatstressMAX1    <- heatstress_units_max1 # use for max chill before mortality??
 heatstressMAX2    <- template  # use for max chill before mortality??
 heatstressMAX2    <- heatstress_units_max2 # use for max chill before mortality??
 heatEXCL          <- template  # EXCL map for heat stress             
 AllEXCL           <- template  # EXCL map for combined stresses (chill,heat,later: moisture)            
  # end NEW
}
cat("made it here line endxxxx \n")
  #  cat("LINE 323: \n")

#### If Pest Event Maps (PEMS) wanted then init PEM rasters ####
if (pems) { 
  if (eggEventDD) {  # must be specified in spp.params file
  eggEvent <- template # raster to track DDs for egg hatch PEMs
  eggEvent <- eggEvent + eggEventDD  # allow extra DDs for event to be triggered
    if (PEMnumgens > 0) {
      PEMe1 <- template  # egg DOYs for when cumDDs > eggEvent threshold  1st Gen
    } 
    if (PEMnumgens > 1) {
      PEMe2 <- template  # egg DOYs for when cumDDs > eggEvent threshold  2nd Gen
    } 
    if (PEMnumgens > 2) {
      PEMe3 <- template  # egg DOYs for when cumDDs > eggEvent threshold  3rd Gen
    } 
    if (PEMnumgens > 3) {
      PEMe4 <- template  # egg DOYs for when cumDDs > eggEvent threshold  4th Gen
    }
  }
  if (larvaeEventDD) {
  larvaeEvent <- template # raster to track DDs for larvae PEMs
  larvaeEvent <- larvaeEvent + larvaeEventDD  # mid-larval dev after ca. xx DDs 
    if (PEMnumgens > 0) {
      PEMl1 <- template  # larval DOYs for when cumDDs > larvaeEvent threshold  1st Gen
    } 
    if (PEMnumgens > 1) {
      PEMl2 <- template  # larval DOYs for when cumDDs > larvaeEvent threshold  2nd Gen
    } 
    if (PEMnumgens > 2) {
      PEMl3 <- template  # larval DOYs for when cumDDs > larvaeEvent threshold  3rd Gen
    } 
    if (PEMnumgens > 3) {
      PEMl4 <- template  # larval DOYs for when cumDDs > larvaeEvent threshold  4th Gen
    } 
  }
  if (pupaeEventDD) {
  pupaeEvent <- template # raster to track DDs for pupal PEMs
  pupaeEvent <- pupaeEvent + pupaeEventDD  # mid-pupal dev after ca. xx DDs 
    if (PEMnumgens > 0) {
      PEMp1 <- template  # pupal DOYs for when cumDDs > pupalEvent threshold  1st Gen
    } 
    if (PEMnumgens > 1) {
      PEMp2 <- template  # pupal DOYs for when cumDDs > pupalEvent threshold  2nd Gen
    } 
    if (PEMnumgens > 2) {
      PEMp3 <- template  # pupal DOYs for when cumDDs > pupalEvent threshold  3rd Gen
    } 
    if (PEMnumgens > 3) {
      PEMp4 <- template  # pupal DOYs for when cumDDs > pupalEvent threshold  4th Gen
    } 
  }
  if (adultEventDD) {
  adultEvent <- template # raster of DOY for adult emerge PEMs
  adultEvent <- adultEvent + adultEventDD # adult emerge occurs after xx DDs during adult stage
    if (PEMnumgens > 0) {
      PEMa1 <- template  # adult DOYs for when cumDDs > adultEvent threshold  1st Gen
    } 
    if (PEMnumgens > 1) {
      PEMa2 <- template  # adult DOYs for when cumDDs > adultEvent threshold  2nd Gen
    } 
    if (PEMnumgens > 2) {
      PEMa3 <- template  # adult DOYs for when cumDDs > adultEvent threshold  3rd Gen
    } 
    if (PEMnumgens > 3) {
      PEMa4 <- template  # adult DOYs for when cumDDs > adultEvent threshold  4th Gen
    } 
  }
}
#### END Pest Event Maps (PEMS) wanted then init PEM rasters ####

#### Using OVLP feature; init maps that track this #### 
if (ovlp > 0.0) {
  eggDev <- template # proportion egg has developed thus far
  eggP <- template # proportion in egg stage
  larvDev <- template # proportion larvae has developed thus far
  larvP <- template # proportion in larval stage
  pupDev <- template # proportion pupae have developed thus far
  pupP <- template # proportion in pupal stage
  adultDev <- template # proportion adult has developed thus far
  adultP <- template # proportion in adult stage
}

#### Init OWSTAGE maps #### 
if (owstage == "OE") {
  OWeggP <- template # proportion in OW egg stage
} else if (owstage == "OL") {
  OWlarvP <- template # proportion in OW larval stage
} else if (owstage == "OP") {
  OWpupP <- template # proportion in OW pup stage
} else if (owstage == "OA") {
  OWadultP <- template # proportion in OW adult stage
} 

#rm(template)
#   does this change with OW stage change ??
#Lifestage: [-1] = OW stage [0] = egg, [1] = larvae, [2] = pupae, [3] = adult

#### Accumulate degree days and reclass cells to NA with temperature exclusion ####
# NA for exclusion means that DD accum cannot occur anymore with incomplete
# generation development and no oviposition to next generation. This may
# be tested for sensitivity by allowing exlusion cells to reset to zero instead.

#### Set up start and stop dates for loop ####
doy <- start_doy
#start <- doy + 1
uniq_list <- unique(sortedlist)
sublist <- uniq_list[doy:end_doy]

#### Initialize stage specific lifestage binary rasters ####
# Limits operations to a mask for cells that are in that lifestage
# This is what allows for pixel by pixel tracking of what lifestage
# that cell is in - move updates for these to end of day

#### Init Lifestage tracking
    if (owstage == "OE") {
      LSOW0 <- Lifestage == -1 
	 } else if (owstage == "OL") {
      LSOW1 <- Lifestage == -1 
	 } else if (owstage == "OP") {
      LSOW2 <- Lifestage == -1 
	 } else if (owstage == "OA") {
      LSOW3 <- Lifestage == -1 
	 }
    LS0 <- Lifestage == 0
    LS1 <- Lifestage == 1
    LS2 <- Lifestage == 2
    LS3 <- Lifestage == 3

# find todays date since we may want to use it for detn freq of output
     today <- Sys.Date()
     showtoday <- format(today, format="%Y%m%d")
     yesterday <- format(today-1, format="%Y%m%d")
     todayplus1 <- format(today+1, format="%Y%m%d")
     todayplus2 <- format(today+2, format="%Y%m%d")
     todayplus3 <- format(today+3, format="%Y%m%d")
     todayplus7 <- format(today+7, format="%Y%m%d")
cat("today and tomorrow etc =",showtoday,todayplus1,todayplus2,todayplus7,"\n")  # *** 
########                         END initializations                               #########

########                 Begin of actual stepping through the model                #########
########                                                                           #########
for (d in sublist) {
	 doy <- doy + 1
    #cat("d and doy =",d,doy,"\n")
#cat("d,doy, filetype and file =",d," ",doy," ",quality[[d]]," ",maxfiles[[d]],"\n")  # *** 

if ("d" == "showtoday") {
   cat("d, doy, showtoday =",d,doy,showtoday,"\n")  # *** 
} else {
   #cat("d, doy, NOT showtoday =",d,doy,showtoday,"\n")  # *** 
}

    if (.Platform$OS.type == "windows") flush.console()
    Sys.sleep(1)
    #Read in that day's PRISM raster files ## DEPRECATED PLAN TO REMOVE USE of PRISM_tmean
    if(calctype=="simple") { # only need these for simple DD calcs
      pattern = paste("(PRISM_tmean_)(.*)(",d,")(_bil.bil)$", sep="")
      temp <- list.files(pattern=pattern,all.files=FALSE, full.names=TRUE)
      #tmean <- raster(temp)
      tmean <- crop(raster(temp),REGION)
    }

#### Prepare & crop to current region: tmax and tmin raster files ####
    temp <- maxfiles[[d]]
    tmax <- crop(raster(temp),REGION)
			  #setwd(output_dir)
           #writeRaster(tmax,paste("Tmax_",d,sep=""), format="GTiff", overwrite=TRUE)
			  #setwd(prism_dir)
	 #dataType(tmax) <- "FLT4S"
    temp <- minfiles[[d]]
    tmin <- crop(raster(temp),REGION)
    tmean <- crop(raster(temp),REGION) # NEW do I need this??
	 if(calctype=="simple") { # only need these for simple DD calcs
	    tmean <- (tmax+tmin)/2
    }
    else {  # always need for NEW exclusion mapping
	    tmean <- (tmax+tmin)/2
    }

#### Loop through Stages; order of stages now read from SPP param file ####
#stgorder <- c("OE","L","P","A","E","F")
#stgorder <- c("OL","P","A","E","L","F")

    for (i in stgorder) {  # Handle stages in the model
    ####  MAIN STEPS FOR EGG STAGE ####
		if (i == "E" | i == "OE") {   # Egg Stage
           if(calctype=="average") { #devel DDs (zero values for temps below LDT)
             dd0tmp <- AvgDD(tmax,tmin,eggLDT,eggUDT)
			    #setwd(output_dir)
             #writeRaster(dd0tmp,paste("dd0tmp_",d,sep=""), format="GTiff", overwrite=TRUE)
             #rast_name <- sprintf("%s%s%s","dd0tmp_",d,".tif")
             #rast_label <- sprintf("%s%s%s",spp," Daily Degree-Days on ",d)
				 #Sys.sleep(0.3)
			    #ConvPNG(rast_name) 
			    #setwd(prism_dir)
           } else if(calctype=="triangle") {
		       dd0tmp <- TriDD(tmax,tmin,eggLDT,eggUDT)
			  } else { # assume (calctype=="simple") 
             dd0tmp <- SimpDD(tmean,eggLDT)
           }

			if(exclusions) { #Calculate lower lethal threshold and exclusion mask
			  # reminder: 0 means exclude; 1 means dont exclude
           #eggmin <- tmin > eggLLT  # old version eggmin
			  # new section to accumulate days below threshold - how to increment
			  eggmindays <- Cond(tmin <= eggLLT, eggmindays+1, eggmindays)
			  eggminEXCL <- Cond(eggmindays >= eggLLDAYS, 0, 1) 
			  # end new section to accumulate days below thresh.
           eggminEXCL[eggminEXCL==0] <- NA  # why does this need to be here??
           #writeRaster(eggmin,paste("EggMin_",d,sep=""), format="GTiff",overwrite=TRUE)
           #Calculate upper lethal threshold and exclusion mask (keep old method for now?)
           #eggmax <- tmax < eggULT
           #eggmax[eggmax==0] <- NA
	
		}
			#Apply exclusions and mask to daily DDs and limit to correct stage
			# NEED TO rethink how exclusions get used; should not be used here in fact
			#if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays DDs
			#  if (i == "OE") { dd0 <- dd0tmp * eggmin * eggmax * LSOW0 }
			#  else if (i == "E") { dd0 <- dd0tmp * eggmin * eggmax * LS0 }
			#}
			#else { # just use lifestage mask
			  if (i == "OE") { dd0 <- dd0tmp * LSOW0 }
			  else if (i == "E") { dd0 <- dd0tmp * LS0 }
         #}
           #Apply exclusions and lifestage mask to daily degree day surface
			  #if (i == "OE") { dd0 <- dd0tmp * eggmin * eggmax * LSOW0
			  #  } else if (i == "E") { dd0 <- dd0tmp * eggmin * eggmax * LS0
			  #}
			 #}
			 #else { # just use lifestage mask
			 # if (i == "OE") { dd0 <- dd0tmp * LSOW0
			 ##   } else if (i == "E") { dd0 <- dd0tmp * LS0
			 # }
          #}

         #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
         #count cells with dd>0
         dd.stat <- cellStats(dd0,stat='max',na.rm=TRUE)
         if (dd.stat > 0) {
           DDaccum <- DDaccum + dd0
			  if (pems) {
			  # set DOY in PEM if curr gen, not set yet, DDs > eventDDs, and Stage is eggs, else keep same value
             if (PEMnumgens > 0) {
			      PEMe1 <- Cond(NumGen==1 & PEMe1 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe1) 
             }
				 if (PEMnumgens > 1) {
			      PEMe2 <- Cond(NumGen==2 & PEMe2 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe2) 
             }
				 if (PEMnumgens > 2) {
			      PEMe3 <- Cond(NumGen==3 & PEMe3 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe3) 
             }
				 if (PEMnumgens > 3) {
			      PEMe4 <- Cond(NumGen==4 & PEMe4 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe4) 
             }
			  }
           #Calculate lifestage progression: for dd accum in correct lifestage, is it >= egg DD threshold?
			  if (i == "OE") {
    #cat("date and doy, i =",d,doy,i,"\n")
              progressOW0 <- (DDaccum * LSOW0) >= OWeggDD
              #writeRaster(progressOW0,paste("ProgressOW0_",d,sep=""), format="GTiff",overwrite=TRUE)
              Lifestage <- Cond(LSOW0 == 1 & progressOW0 == 1, 1, Lifestage)
              #Reset the DDaccum cells to zero for cells that progressed to next lifestage
              DDaccum <- Cond(progressOW0 == 1,(DDaccum - OWeggDD) * LSOW0, DDaccum)
			  } else if (i == "E") {
    #cat("date and doy, i =",d,doy,i,"\n")
              progress0 <- (DDaccum * LS0) >= eggDD
              #writeRaster(progress0,paste("Progress0_",d,sep=""), format="GTiff",overwrite=TRUE)
              Lifestage <- Cond(LS0 == 1 & progress0 == 1, 1, Lifestage)
              #Reset the DDaccum cells to zero for cells that progressed to next lifestage
              DDaccum <- Cond(progress0 == 1,(DDaccum - eggDD) * LS0, DDaccum)
			  }
         }

    ####  MAIN STEPS FOR LARVAL STAGE ####
		} else if (i == "L" | i == "OL") {  # Larval Stage
         #developmental degree days
         if(calctype=="average") {
           dd1tmp <- AvgDD(tmax,tmin,larvaeLDT,larvaeUDT)
         } else if(calctype=="triangle") {
		     dd1tmp <- TriDD(tmax,tmin,larvaeLDT,larvaeUDT)
			} else { # assume (calctype=="simple") 
           dd1tmp <- SimpDD(tmean,larvaeLDT)
         }
         ddtotal <- ddtotal + dd1tmp #LEN : Accumulate total degree days for the year for larvae

			if(exclusions) { # Calculate lower lethal threshold and exclusion mask
           #larvaemin <- tmin > larvaeLLT
			  # new section to accumulate days below threshold - how to increment
			  larvaemindays <- Cond(tmin <= larvaeLLT, larvaemindays+1, larvaemindays)
			  larvaeminEXCL <- Cond(larvaemindays >= larvaeLLDAYS, 0, 1) 
			  # end new section to accumulate days below thresh.
           larvaeminEXCL[larvaeminEXCL==0] <- NA
           #Calculate upper lethal threshold and exclusion mask
           #larvaemax <- tmax < larvaeULT
           #larvaemax[larvaemax == 0] <- NA
           #writeRaster(larvaemax,paste("Larvaemax_",d,sep=""), format="GTiff",overwrite=TRUE)
          } 
          else if(exclusions_stressunits) {
	 # NEW version exclusions - applies to all stages
    chillmask <- tmin < chillstressTHRESH  # make todays chill mask
    chillstress <- chillmask * abs(chillstressTHRESH - tmin) # compute todays Chill stress DDs
    chillunitsCUM <- chillunitsCUM + chillstress
	 # ASSUME NEW -2=severe -1=mod 0=none throughout
    chillEXCL <- Cond(chillunitsCUM >= chillstressMAX2,-2,Cond(chillunitsCUM >= chillstressMAX1,-1,0))
    ##-- Heat Stress Accumulation--------------------
    heatmask <- tmax > heatstressTHRESH  # make todays heat mask
    heatstress <- heatmask * abs(tmax - heatstressTHRESH) # compute todays heat stress DDs
    heatunitsCUM <- heatunitsCUM + heatstress
    heatEXCL <- Cond(heatunitsCUM >= heatstressMAX2,-2,Cond(heatunitsCUM >= heatstressMAX1,-1,0))
	 AllEXCL <- Cond((chillEXCL == 0) & (heatEXCL == 0),0,
	            Cond((chillEXCL == -1) & (heatEXCL >= -1),-1,
	            Cond((chillEXCL >= -1) & (heatEXCL == -1),-1,-2)))
	 # end NEW exclusions
			}
			#Apply exclusions and mask to daily DDs and limit to correct stage
			#if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays DDs
			#  if (i == "OL") { dd1 <- dd1tmp * larvaemin * larvaemax * LSOW1 }
			#  else if (i == "L") { dd1 <- dd1tmp * larvaemin * larvaemax * LS1 }
			#}
			#else { # just use lifestage mask
			  if (i == "OL") { dd1 <- dd1tmp * LSOW1 }
			  else if (i == "L") { dd1 <- dd1tmp * LS1 }
         #}

         #Accumulate degree days, if dd1 > 0 otherwise exclusion masks get applied to this generation.
         dd.stat <- cellStats(dd1,stat='max',na.rm=TRUE)
         if (dd.stat > 0) {
           DDaccum <- DDaccum + dd1
			  if (pems & larvaeEventDD) {
			  # set DOY in PEM if curr. gen, not set yet, DDs > eventDDs, Stage is larvae, else keep same value
             if (PEMnumgens > 0) {
			      PEMl1 <- Cond(NumGen==1 & PEMl1 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl1) 
             }
				 if (PEMnumgens > 1) {
			      PEMl2 <- Cond(NumGen==2 & PEMl2 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl2) 
             }
				 if (PEMnumgens > 2) {
			      PEMl3 <- Cond(NumGen==3 & PEMl3 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl3) 
             }
				 if (PEMnumgens > 3) {
			      PEMl4 <- Cond(NumGen==4 & PEMl4 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl4) 
             }
           }
           #Calculate lifestage progression: for dd accum in correct lifestage, is it >= egg DD threshold?
			  if (i == "OL") {
              progressOW1 <- (DDaccum * LSOW1) >= OWlarvaeDD
              #writeRaster(progressOW1,paste("ProgressOW1_",d,sep=""), format="GTiff",overwrite=TRUE)
              Lifestage <- Cond(LSOW1 == 1 & progressOW1 == 1, 2, Lifestage)
              #Reset the DDaccum cells to zero for cells that progressed to next lifestage
              DDaccum <- Cond(progressOW1 == 1,(DDaccum - OWlarvaeDD) * LSOW1, DDaccum)
			  } else if (i == "L") {
              progress1 <- (DDaccum * LS1) >= larvaeDD
              #writeRaster(progress1,paste("Progress1_",d,sep=""), format="GTiff",overwrite=TRUE)
              Lifestage <- Cond(LS1 == 1 & progress1 == 1, 2, Lifestage)
              DDaccum <- Cond(progress1 == 1,(DDaccum - larvaeDD) * LS1, DDaccum)
           }
			}
    ####  MAIN STEPS FOR PUPAL STAGE ####
		} else if (i == "P" | i == "OP") {   # Pupal Stage
         #developmental degree days
         if(calctype=="average") {
           dd2tmp <- AvgDD(tmax,tmin,pupaeLDT,pupaeUDT)
         } else if(calctype=="triangle") {
		     dd2tmp <- TriDD(tmax,tmin,pupaeLDT,pupaeUDT)
			} else { # assume (calctype=="simple") 
           dd2tmp <- SimpDD(tmean,pupaeLDT)
         }
         #Apply exclus (none for pupae stage) & lifestage mask to daily DD surface and limit to correct stage
			  # new section to accumulate days below threshold - how to increment
			  #pupaemindays <- Cond(tmin <= pupaeLLT, pupaemindays+1, pupaemindays)
			  #pupaeminEXCL <- Cond(pupaemindays >= pupaeLLDAYS, 0, 1) 
			  # end new section to accumulate days below thresh.
           #pupaemin[pupsemin==0] <- NA
           #Calculate upper lethal threshold and exclusion mask
           #pupaemax <- tmax < pupaeULT
           #pupaemax[pupaemax == 0] <- NA
			if (i == "OP") { dd2 <- dd2tmp * LSOW2
			} else if (i == "P") { dd2 <- dd2tmp * LS2
			}

         #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
         dd.stat <- cellStats(dd2,stat='max',na.rm=TRUE)
         if (dd.stat > 0) {
           DDaccum <- DDaccum + dd2
			  if (pems & pupaeEventDD) {
			  # set DOY in PEM if curr. gen, not set yet, DDs > eventDDs, Stage is pupa, else keep same value
             if (PEMnumgens > 0) {
			      PEMp1 <- Cond(NumGen==1 & PEMp1 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp1) 
				 }
             if (PEMnumgens > 1) {
			      PEMp2 <- Cond(NumGen==2 & PEMp2 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp2) 
				 }
             if (PEMnumgens > 2) {
			      PEMp3 <- Cond(NumGen==3 & PEMp3 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp3) 
				 }
             if (PEMnumgens > 3) {
			      PEMp4 <- Cond(NumGen==4 & PEMp4 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp4) 
				 }
           }
           #Calculate lifestage progression: for dd accum in correct lifestage, is it >= pupae DD threshold?
			  if (i == "OP") {
              progressOW2 <- (DDaccum * LSOW2) >= OWpupDD
              #writeRaster(progressOW2,paste("ProgressOW2_",d,sep=""), format="GTiff",overwrite=TRUE)
              Lifestage <- Cond(LSOW2 == 1 & progressOW2 == 1, 3, Lifestage)
              #Reset the DDaccum cells to zero for cells that progressed to next lifestage
              DDaccum <- Cond(progressOW2 == 1,(DDaccum - OWpupDD) * LSOW2, DDaccum)
			  } else if (i == "P") {
              progress2 <- (DDaccum * LS2) >= pupDD
              #writeRaster(progress2,paste("Progress2_",d,sep=""), format="GTiff",overwrite=TRUE)
              Lifestage <- Cond(LS2 == 1 & progress2 == 1, 3, Lifestage)
              DDaccum <- Cond(progress2 == 1,(DDaccum - pupDD) * LS2, DDaccum)
           }
         }
    ####  MAIN STEPS FOR ADULT STAGE ####
		} else if (i == "A" | i == "OA") {  # Adult stage, or time to 50% oviposition
         #developmental degree days
         if(calctype=="average") {
           dd3tmp <- AvgDD(tmax,tmin,adultLDT,adultUDT)
         } else if(calctype=="triangle") {
		     dd3tmp <- TriDD(tmax,tmin,adultLDT,adultUDT)
			} else { # assume (calctype==simple) 
           dd3tmp <- SimpDD(tmean,adultLDT)
         }

			if(exclusions) { # Calculate lower lethal threshold and exclusion mask
           adultmin <- tmin > adultLLT
			  # new section to accumulate days below threshold - how to increment
			  adultmindays <- Cond(tmin <= adultLLT, adultmindays+1, adultmindays)
			  adultminEXCL <- Cond(adultmindays >= adultLLDAYS, 0, 1) 
           adultminEXCL[adultminEXCL==0] <- NA
    #      writeRaster(adultmin,paste("Admin_",d,sep=""), format="GTiff",overwrite=TRUE)
			  # NOW sum up all EXCL maps 
			  allminEXCL <- (eggminEXCL * larvaeminEXCL * adultminEXCL)
			  # end new section to accumulate days below thresh.
			  #Apply exclusions and mask to daily DDs and limit to correct stage
			}
			#Apply exclusions and mask to daily DDs and limit to correct stage
			#if(exclusions && exclusionsrealtime) { # apply the exclusion maps to todays DDs
			#  if (i == "OA") { dd3 <- dd3tmp * adultmin * adultmax * LSOW3 }
			#  else if (i == "A") { dd3 <- dd3tmp * adultmin * adultmax * LS3 }
			#}
			#else { # just use lifestage mask
			  if (i == "OA") { dd3 <- dd3tmp * LSOW3 }
			  else if (i == "A") { dd3 <- dd3tmp * LS3 }
         #}
			#  if (i == "OA") { dd3 <- dd3tmp * adultmin * LSOW3 }
			#  else if (i == "A") { dd3 <- dd3tmp * adultmin * LS3 }
			# }
			# else { # just use lifestage mask
			#  if (i == "OA") { dd3 <- dd3tmp * LSOW3 }
			#  else if (i == "A") { dd3 <- dd3tmp * LS3 }
         # }

         #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
         dd.stat <- cellStats(dd3,stat='max',na.rm=TRUE)
			print(dd.stat)
         if (dd.stat > 0) {
           DDaccum <- DDaccum + dd3
			  if(pems & adultEventDD) { # record DOYs in PEMS for adult events
             if (PEMnumgens > 0) {
			      PEMa1 <- Cond(NumGen==1 & PEMa1 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa1) 
             }
				 if (PEMnumgens > 1) {
			      PEMa2 <- Cond(NumGen==2 & PEMa2 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa2) 
             }
				 if (PEMnumgens > 2) {
			      PEMa3 <- Cond(NumGen==3 & PEMa3 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa3) 
             }
				 if (PEMnumgens > 3) {
			      PEMa4 <- Cond(NumGen==4 & PEMa4 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa4) 
			    }
			  }

           #Calculate lifestage progression: for dd accum in correct lifestage, is it >= adult DD threshold?
			  if (i == "OA") {
              progressOW3 <- (DDaccum * LSOW3) >= OWadultDD
              #writeRaster(progressOW3,paste("ProgressOW3_",d,sep=""), format="GTiff",overwrite=TRUE)
              DDaccum <- Cond(progressOW3 == 1,(DDaccum - OWadultDD) * LSOW3, DDaccum)
              progressOW3[is.na(progressOW3)] <- template[is.na(progressOW3)]
              Lifestage <- Cond(LSOW3 == 1 & progressOW3 == 1, 0, Lifestage)
              Lifestage[is.na(Lifestage)] <- template[is.na(Lifestage)]
              #Increment NumGen + 1
              NumGen <- NumGen + progressOW3
              #writeRaster(NumGen,paste("NumGen",d,sep=""), format="GTiff",overwrite=TRUE)
			  } else if (i == "A") {
              progress3 <- (DDaccum * LS3) >= adultDD
              #writeRaster(progress3,paste("Progress3_",d,sep=""), format="GTiff",overwrite=TRUE)
              #Reset the DDaccum cells to zero for cells that progressed to next lifestage
              DDaccum <- Cond(progress3 == 1,(DDaccum - adultDD) * LS3, DDaccum)
              #Remove masking effect from progress counter so it doesn't perpetuate through NumGen raster
              progress3[is.na(progress3)] <- template[is.na(progress3)]
              Lifestage <- Cond(LS3 == 1 & progress3 == 1,0, Lifestage)
              Lifestage[is.na(Lifestage)] <- template[is.na(Lifestage)]
              #Increment NumGen + 1
              NumGen <- NumGen + progress3
              #writeRaster(NumGen,paste("NumGen",d,sep=""), format="GTiff",overwrite=TRUE)
			  }

         }

    ####  MAIN STEPS FOR END OF EACH DAY ####
      } else if (i == "F") { # end of the day placeholder 

          # update lifestage masks
			  if (owstage == "OE") {
             LSOW0 <- Lifestage == -1
			  } else if (owstage == "OL") {
             LSOW1 <- Lifestage == -1
			  } else if (owstage == "OP") {
             LSOW2 <- Lifestage == -1
			  } else if (owstage == "OA") {
             LSOW3 <- Lifestage == -1
			  }
          LS0 <- Lifestage == 0
          LS1 <- Lifestage == 1
          LS2 <- Lifestage == 2
          LS3 <- Lifestage == 3

        # Need to do anything else (overlap and write maps)?
           if (doy %% 30 == 0 && monthlymaps == 1) { ## monthly
             cat("### Monthly maps output ### d, doy, i =",d,doy,i,"\n")
             do_ovlp <- 1     # OK to calc ovlp
             do_maps <- 1     # OK to write maps
           } else if (doy %% 14 == 0 && biweeklymaps == 1) { ## biweekly
             cat("### Biweekly maps output ### d, doy, i =",d,doy,i,"\n")
             do_ovlp <- 1     # OK to calc ovlp
             do_maps <- 1     # OK to write maps
           } else if (doy %% 10 == 0 && dekadmaps == 1) { ## every 10 days
             cat("### Dekad maps output ### d, doy, i =",d,doy,i,"\n")
             do_ovlp <- 1     # OK to calc ovlp
             do_maps <- 1     # OK to write maps
           } else if (doy %% 7 == 0 && weeklymaps == 1) { ## weekly
             cat("### Weekly maps output ### d, doy, i =",d,doy,i,"\n")
             do_ovlp <- 1     # OK to calc ovlp
             do_maps <- 1     # OK to write maps
           } else if (doy %% 2 == 0 && evendaymaps == 1) { ## every other day
             cat("### Even day maps output ### d, doy, i =",d,doy,i,"\n")
             do_ovlp <- 1     # OK to calc ovlp
             do_maps <- 1     # OK to write maps
           } else if (dailymaps == 1) { ## must be every day
             cat("### Daily maps output ### d, doy, i =",d,doy,i,"\n")
             do_ovlp <- 1     # OK to calc ovlp
             do_maps <- 1     # OK to write maps
           } else { ## no need to do ovlp or write maps
             do_ovlp <- 0     # NO to calc ovlp
             do_maps <- 0     # NO to write maps
           }

        if (ovlp > 0.0 && do_ovlp == 1) {  # use overlap feature which just ramps a transition
           OWCurr_Overlap=function(d,c){
           	 return(Cond(d <= c,1,
           	 Cond(d <= (1-c),1,
           	 0.5*(1+(1-d)/c))))
           }
           Curr_Overlap=function(d,c){
           	 return(Cond(d <= c,0.5*(1+d/c),
           	 Cond(d <= (1-c),1,
           	 0.5*(1+(1-d)/c))))
           }
           Prev_Overlap=function(d,c){
           	 return(Cond(d <= c,0.5*(1-d/c),0))
           }
           Next_Overlap=function(d,c){
           	 return(Cond(d <= (1-c),0,
	          0.5*(1-(1-d)/c)))
           }
           #progress1 <- (DDaccum * LS1) >= larvaeDD 
			  if (owstage == "OE") {
             #cat("date and doy,owstage, i =",d,doy,owstage,i,"\n")
             OWeggP <- template # proportion in OW egg stage
			    OWeggDev <- (DDaccum * LSOW0)/OWeggDD
             OWeggP <- OWeggP + LSOW0 * OWCurr_Overlap(OWeggDev,ovlp)
			    # no prev stage overlap for OW stage
             larvP <- larvP + LSOW0 * Next_Overlap(OWeggDev,ovlp)
				 #setwd(output_dir)
             #writeRaster(OWeggP,paste("OWEggP_",d,sep=""), format="GTiff",overwrite=TRUE)
             #writeRaster(OWeggDev,paste("OWEggDev_",d,sep=""), format="GTiff",overwrite=TRUE)
				 #setwd(prism_dir)
			  } else if (owstage == "OL") {
             OWlarvP <- template # proportion in OW larval stage
			    OWlarvDev <- (DDaccum * LSOW1)/OWlarvaeDD
             OWlarvP <- OWlarvP + LSOW1 * OWCurr_Overlap(OWlarvDev,ovlp)
             pupP <- pupP + LSOW1 * Next_Overlap(OWlarvDev,ovlp)
			  } else if (owstage == "OP") {
             OWpupP <- template # proportion in OW pup stage
			    OWpupDev <- (DDaccum * LSOW2)/OWpupDD
             OWpupP <- OWpupP + LSOW2 * OWCurr_Overlap(OWpupDev,ovlp)
             adultP <- adultP + LSOW2 * Next_Overlap(OWpupDev,ovlp)
			  } else if (owstage == "OA") {
             OWadultP <- template # proportion in OW adult stage
			    OWadultDev <- (DDaccum * LSOW3)/OWadultDD
             OWadultP <- OWadultP + LSOW3 * OWCurr_Overlap(OWadultDev,ovlp)
             eggP <- eggP + LSOW3 * Next_Overlap(OWadultDev,ovlp)
			  }
			  eggP <- template
			  larvP <- template
			  pupP <- template
			  adultP <- template

			  eggDev <- (DDaccum * LS0)/eggDD
           eggP <- eggP + LS0 * Curr_Overlap(eggDev,ovlp)
           adultP <- adultP + LS0 * Prev_Overlap(eggDev,ovlp)
           larvP <- larvP + LS0 * Next_Overlap(eggDev,ovlp)

			  larvDev <- (DDaccum * LS1)/larvaeDD
           larvP <- larvP + LS1 * Curr_Overlap(larvDev,ovlp)
           eggP <- eggP + LS1 * Prev_Overlap(larvDev,ovlp)
           pupP <- pupP + LS1 * Next_Overlap(larvDev,ovlp)

			  pupDev <- (DDaccum * LS2)/pupDD
           pupP <- pupP + LS2 * Curr_Overlap(pupDev,ovlp)
           larvP <- larvP + LS2 * Prev_Overlap(pupDev,ovlp)
           adultP <-  adultP + LS2 * Next_Overlap(pupDev,ovlp)

			  adultDev <- (DDaccum * LS3)/adultDD
           adultP <- adultP + LS3 * Curr_Overlap(adultDev,ovlp)
           pupP <- pupP + LS3 * Prev_Overlap(adultDev,ovlp)
           eggP <- eggP + LS3 * Next_Overlap(adultDev,ovlp)
			  
           # add OW stageP to regular (curr and next) stageP
			  if (owstage == "OE") {
           eggP <- eggP + LSOW0 * OWeggP
           larvP <- larvP + LSOW0 * (1 - OWeggP)
			  } else if (owstage == "OL") {
           larvP <- larvP + LSOW1 * OWlarvP
           pupP <- pupP + LSOW1 * (1 - OWlarvP)
			  } else if (owstage == "OP") {
           pupP <- pupP + LSOW2 * OWpupP
           adultP <- adultP + LSOW2 * (1 - OWpupP)
			  } else if (owstage == "OA") {
           adultP <- adultP + LSOW3 * OWadultP
           eggP <- eggP + LSOW3 * (1 - OWadultP)
			  }

           # new Nov 2015 tally stages plus gens
           StageCount <- Lifestage + (NumGen * 4)

           #cat("today and showtoday,d,doy,i =",today,showtoday,d,doy,i,"\n")

             if (do_maps == 1) { ## today
               cat("### its today - write maps today ### d, doy, i =",d,doy,i,"\n")
	            WriteMaps(d)
			    } 
           }  # ovlp > 0 && do_ovlp == 1
           else if (do_maps == 1) { # no overlap but maybe write maps today
             cat("### its today - write maps today ### d, doy, i =",d,doy,i,"\n")
             WriteMaps(d)
           }
	    }  # lifestage F or 5 (end of day calcs)
   }  # lifestage for loop
}  # daily loop

#### Model Done - final map production ####
setwd(output_dir)
rast_name <- sprintf("%s%s%s%s",spp,"_DD_",d,".tif")
writeRaster(ddtotal,paste(spp,"_DD_",d,sep=""), format="GTiff",overwrite=TRUE)
				 Sys.sleep(0.3)
ConvPNG(rast_name)
rast_name <- sprintf("%s%s%s%s",spp,"_NumGen_",d,".tif")
writeRaster(NumGen,paste(spp,"_NumGen_",d,sep=""), format="GTiff",overwrite=TRUE)
				 Sys.sleep(0.3)
ConvPNG(rast_name)

if (exclusions_stressunits) {
	# Final Stress Unit Maps:
             rast_name <- sprintf("%s%s%s","Chill_Stress_Units_",d,".tif")
             writeRaster(chillunitsCUM,paste("Chill_Stress_Units_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             rast_name <- sprintf("%s%s%s","Chill_Stress_Excl_",d,".tif")
             writeRaster(chillEXCL,paste("Chill_Stress_Excl_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             rast_name <- sprintf("%s%s%s","Heat_Stress_Units_",d,".tif")
             writeRaster(heatunitsCUM,paste("Heat_Stress_Units_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             rast_name <- sprintf("%s%s%s","Heat_Stress_Excl_",d,".tif")
             writeRaster(heatEXCL,paste("Heat_Stress_Excl_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             rast_name <- sprintf("%s%s%s","All_Stress_Excl_",d,".tif")
             writeRaster(AllEXCL,paste("All_Stress_Excl_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
	#
				 # Exclusions included: -2=severe, -1=moderate stress
             NumGenEXCL1 <- Cond((AllEXCL > -2),NumGen,-2)
             rast_name <- sprintf("%s%s%s","NumGenExcl1_",d,".tif")
             writeRaster(NumGenEXCL1,paste("NumGenExcl1_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
             NumGenEXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,NumGen))
             rast_name <- sprintf("%s%s%s","NumGenExcl2_",d,".tif")
             writeRaster(NumGenEXCL2,paste("NumGenExcl2_",d,sep=""), format="GTiff", overwrite=TRUE)
             ConvPNG(rast_name)
#THIS WAS DONE FOR LBAM AND IS NOT GENERIC SO DROP FOR NOW
#EggHeatMort <- 100*(1-exp(-0.0609*MaxHeat))
#rast_name <- sprintf("%s%s%s%s",spp,"_EggHeatMort_",d,".tif")
#writeRaster(EggHeatMort,paste(spp,"_EggHeatMort_",d,sep=""), format="GTiff",overwrite=TRUE)
#				 Sys.sleep(0.3)
#ConvPNG(rast_name)
#
#LarvHeatMort <- 100*(1-exp(-0.1046*MaxHeat))
#rast_name <- sprintf("%s%s%s%s",spp,"_LarvHeatMort_",d,".tif")
#writeRaster(LarvHeatMort,paste(spp,"_LarvHeatMort_",d,sep=""), format="GTiff",overwrite=TRUE)
#				 Sys.sleep(0.3)
#ConvPNG(rast_name)
#
#PupHeatMort <- 100*(1-exp(-0.133*MaxHeat))
#rast_name <- sprintf("%s%s%s%s",spp,"_PupHeatMort_",d,".tif")
#writeRaster(PupHeatMort,paste(spp,"_PupHeatMort_",d,sep=""), format="GTiff",overwrite=TRUE)
#				 Sys.sleep(0.3)
#ConvPNG(rast_name)
# 
#AdHeatMort <- 100*(1-exp(-0.0473*MaxHeat))
#rast_name <- sprintf("%s%s%s%s",spp,"_AdHeatMort_",d,".tif")
#writeRaster(AdHeatMort,paste(spp,"_AdHeatMort_",d,sep=""), format="GTiff",overwrite=TRUE)
#				 Sys.sleep(0.3)
#ConvPNG(rast_name)
#
#LarvChillMort <- 100*(1-exp(-0.024*MaxChill))
#rast_name <- sprintf("%s%s%s%s",spp,"_LarvChillMort_",d,".tif")
#writeRaster(LarvChillMort,paste(spp,"_LarvChillMort_",d,sep=""), format="GTiff",overwrite=TRUE)
#				 Sys.sleep(0.3)
#ConvPNG(rast_name)
#
#ClimateSuit <- 100 * (1-EggHeatMort/100) * (1-LarvHeatMort/100) * (1-PupHeatMort/100) * (1-AdHeatMort/100) * (1-LarvChillMort/100)
#rast_name <- sprintf("%s%s%s%s",spp,"_ClimateSuit_",d,".tif")
#writeRaster(ClimateSuit,paste(spp,"_ClimateSuit_",d,sep=""), format="GTiff",overwrite=TRUE)
#				 Sys.sleep(0.3)
#ConvPNG(rast_name)
##writeRaster(ClimateSuit,"LBAM_ClimateSuitability_2006-2015.tif",format="GTiff",overwrite=TRUE)
#end NEW Apply mortality functions
}



setwd(prism_dir)
if(pems) {  # should PEM maps only at end of year; may not be correct if not run full year
  setwd(output_dir)
  if (mapE == 1) {
    if (PEMnumgens > 0) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMe1",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",eggEventLabel," Gen. 1")
      writeRaster(PEMe1,paste(spp,"_PEMe1",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 1st Gen eggstage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMe1EXCL1 <- Cond((AllEXCL > -2),PEMe1,-2)
         rast_name <- sprintf("%s%s%s","PEMe1Excl1_",d,".tif")
         writeRaster(PEMe1EXCL1,paste("PEMe1Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMe1EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMe1))
         rast_name <- sprintf("%s%s%s","PEMe1Excl2_",d,".tif")
         writeRaster(PEMe1EXCL2,paste("PEMe1Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
	 }
    if (PEMnumgens > 1) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMe2",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",eggEventLabel," Gen. 2")
      writeRaster(PEMe2,paste(spp,"_PEMe2",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 2nd Gen eggstage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMe2EXCL1 <- Cond((AllEXCL > -2),PEMe2,-2)
         rast_name <- sprintf("%s%s%s","PEMe2Excl1_",d,".tif")
         writeRaster(PEMe2EXCL1,paste("PEMe2Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMe2EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMe2))
         rast_name <- sprintf("%s%s%s","PEMe2Excl2_",d,".tif")
         writeRaster(PEMe2EXCL2,paste("PEMe2Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
	 }
    if (PEMnumgens > 2) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMe3",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",eggEventLabel," Gen. 3")
      writeRaster(PEMe3,paste(spp,"_PEMe3",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 3rd Gen eggstage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMe3EXCL1 <- Cond((AllEXCL > -2),PEMe3,-2)
         rast_name <- sprintf("%s%s%s","PEMe3Excl1_",d,".tif")
         writeRaster(PEMe3EXCL1,paste("PEMe3Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMe3EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMe3))
         rast_name <- sprintf("%s%s%s","PEMe3Excl2_",d,".tif")
         writeRaster(PEMe3EXCL2,paste("PEMe3Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
	 }
    if (PEMnumgens > 3) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMe4",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",eggEventLabel," Gen. 4")
      writeRaster(PEMe4,paste(spp,"_PEMe4",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 4th Gen eggstage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMe4EXCL1 <- Cond((AllEXCL > -2),PEMe4,-2)
         rast_name <- sprintf("%s%s%s","PEMe4Excl1_",d,".tif")
         writeRaster(PEMe4EXCL1,paste("PEMe4Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMe4EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMe4))
         rast_name <- sprintf("%s%s%s","PEMe4Excl2_",d,".tif")
         writeRaster(PEMe4EXCL2,paste("PEMe4Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
	 }
  }
  if (mapL == 1) {
    if (PEMnumgens > 0) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMl1",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",larvaeEventLabel," Gen. 1")
      #rast_label <- sprintf("%s%s",fullname," - Date of Larval Stage Gen. 1")
      writeRaster(PEMl1,paste(spp,"_PEMl1",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 1st Gen larval stage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMl1EXCL1 <- Cond((AllEXCL > -2),PEMl1,-2)
         rast_name <- sprintf("%s%s%s","PEMl1Excl1_",d,".tif")
         writeRaster(PEMl1EXCL1,paste("PEMl1Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMl1EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMl1))
         rast_name <- sprintf("%s%s%s","PEMl1Excl2_",d,".tif")
         writeRaster(PEMl1EXCL2,paste("PEMl1Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
    }
	 if (PEMnumgens > 1) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMl2",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",larvaeEventLabel," Gen. 2")
      writeRaster(PEMl2,paste(spp,"_PEMl2",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 2nd Gen larval stage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMl2EXCL1 <- Cond((AllEXCL > -2),PEMl2,-2)
         rast_name <- sprintf("%s%s%s","PEMl2Excl1_",d,".tif")
         writeRaster(PEMl2EXCL1,paste("PEMl2Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMl2EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMl2))
         rast_name <- sprintf("%s%s%s","PEMl2Excl2_",d,".tif")
         writeRaster(PEMl2EXCL2,paste("PEMl2Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
    }
	 if (PEMnumgens > 2) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMl3",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",larvaeEventLabel," Gen. 3")
      writeRaster(PEMl3,paste(spp,"_PEMl3",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 3rd Gen larval stage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMl3EXCL1 <- Cond((AllEXCL > -2),PEMl3,-2)
         rast_name <- sprintf("%s%s%s","PEMl3Excl1_",d,".tif")
         writeRaster(PEMl3EXCL1,paste("PEMl3Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMl3EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMl3))
         rast_name <- sprintf("%s%s%s","PEMl3Excl2_",d,".tif")
         writeRaster(PEMl3EXCL2,paste("PEMl3Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
    }
	 if (PEMnumgens > 3) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMl4",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",larvaeEventLabel," Gen. 4")
      writeRaster(PEMl4,paste(spp,"_PEMl4",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 4th Gen larval stage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMl4EXCL1 <- Cond((AllEXCL > -2),PEMl4,-2)
         rast_name <- sprintf("%s%s%s","PEMl4Excl1_",d,".tif")
         writeRaster(PEMl4EXCL1,paste("PEMl4Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMl4EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMl4))
         rast_name <- sprintf("%s%s%s","PEMl4Excl2_",d,".tif")
         writeRaster(PEMl4EXCL2,paste("PEMl4Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
    }
  }
  if (mapP == 1) {
    if (PEMnumgens > 0) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMp1",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",pupaeEventLabel," Gen. 1")
      writeRaster(PEMp1,paste(spp,"_PEMp1",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 1st Gen pupal stage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMp1EXCL1 <- Cond((AllEXCL > -2),PEMp1,-2)
         rast_name <- sprintf("%s%s%s","PEMp1Excl1_",d,".tif")
         writeRaster(PEMp1EXCL1,paste("PEMp1Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMp1EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMp1))
         rast_name <- sprintf("%s%s%s","PEMp1Excl2_",d,".tif")
         writeRaster(PEMp1EXCL2,paste("PEMp1Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}

    }
	 if (PEMnumgens > 1) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMp2",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",pupaeEventLabel," Gen. 2")
      writeRaster(PEMp2,paste(spp,"_PEMp2",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 2nd Gen pupal stage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMp2EXCL1 <- Cond((AllEXCL > -2),PEMp2,-2)
         rast_name <- sprintf("%s%s%s","PEMp2Excl1_",d,".tif")
         writeRaster(PEMp2EXCL1,paste("PEMp2Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMp2EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMp2))
         rast_name <- sprintf("%s%s%s","PEMp2Excl2_",d,".tif")
         writeRaster(PEMp2EXCL2,paste("PEMp2Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
    }
	 if (PEMnumgens > 2) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMp3",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",pupaeEventLabel," Gen. 3")
      writeRaster(PEMp3,paste(spp,"_PEMp3",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 3rd Gen pupal stage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMp3EXCL1 <- Cond((AllEXCL > -2),PEMp3,-2)
         rast_name <- sprintf("%s%s%s","PEMp3Excl1_",d,".tif")
         writeRaster(PEMp3EXCL1,paste("PEMp3Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMp3EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMp3))
         rast_name <- sprintf("%s%s%s","PEMp3Excl2_",d,".tif")
         writeRaster(PEMp3EXCL2,paste("PEMp3Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
    }
	 if (PEMnumgens > 3) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMp4",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",pupaeEventLabel," Gen. 4")
      writeRaster(PEMp4,paste(spp,"_PEMp4",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 4th Gen pupal stage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMp4EXCL1 <- Cond((AllEXCL > -2),PEMp4,-2)
         rast_name <- sprintf("%s%s%s","PEMp4Excl1_",d,".tif")
         writeRaster(PEMp4EXCL1,paste("PEMp4Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMp4EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMp4))
         rast_name <- sprintf("%s%s%s","PEMp4Excl2_",d,".tif")
         writeRaster(PEMp4EXCL2,paste("PEMp4Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
    }
  }
  if (mapA == 1) {
    if (PEMnumgens > 0) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMa1",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",adultEventLabel," Gen. 1")
	   writeRaster(PEMa1,paste(spp,"_PEMa1",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 1st Gen adult stage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMa1EXCL1 <- Cond((AllEXCL > -2),PEMa1,-2)
         rast_name <- sprintf("%s%s%s","PEMa1Excl1_",d,".tif")
         writeRaster(PEMa1EXCL1,paste("PEMa1Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMa1EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMa1))
         rast_name <- sprintf("%s%s%s","PEMa1Excl2_",d,".tif")
         writeRaster(PEMa1EXCL2,paste("PEMa1Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
    }
    if (PEMnumgens > 1) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMa2",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",adultEventLabel," Gen. 2")
      writeRaster(PEMa2,paste(spp,"_PEMa2",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 2nd Gen adult stage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMa2EXCL1 <- Cond((AllEXCL > -2),PEMa2,-2)
         rast_name <- sprintf("%s%s%s","PEMa2Excl1_",d,".tif")
         writeRaster(PEMa2EXCL1,paste("PEMa2Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMa2EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMa2))
         rast_name <- sprintf("%s%s%s","PEMa2Excl2_",d,".tif")
         writeRaster(PEMa2EXCL2,paste("PEMa2Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
    }
    if (PEMnumgens > 2) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMa3",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",adultEventLabel," Gen. 3")
      writeRaster(PEMa3,paste(spp,"_PEMa3",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 3rd Gen adult stage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMa3EXCL1 <- Cond((AllEXCL > -2),PEMa3,-2)
         rast_name <- sprintf("%s%s%s","PEMa3Excl1_",d,".tif")
         writeRaster(PEMa3EXCL1,paste("PEMa3Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMa3EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMa3))
         rast_name <- sprintf("%s%s%s","PEMa3Excl2_",d,".tif")
         writeRaster(PEMa3EXCL2,paste("PEMa3Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
    }
    if (PEMnumgens > 3) {
      rast_name <- sprintf("%s%s%s",spp,"_PEMa4",".tif")
      #rast_label <- sprintf("%s%s%s%s",fullname," - Date of ",adultEventLabel," Gen. 4")
      writeRaster(PEMa4,paste(spp,"_PEMa4",sep=""), format="GTiff",overwrite=TRUE)
      ConvPNG(rast_name)
      if (exclusions_stressunits) { # 4th Gen adult stage includ exclusions
		 # Exclusions included: -2=severe, -1=moderate stress
         PEMa4EXCL1 <- Cond((AllEXCL > -2),PEMa4,-2)
         rast_name <- sprintf("%s%s%s","PEMa4Excl1_",d,".tif")
         writeRaster(PEMa4EXCL1,paste("PEMa4Excl1_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
         PEMa4EXCL2 <- Cond((AllEXCL == -2),-2, Cond((AllEXCL == -1),-1,PEMa4))
         rast_name <- sprintf("%s%s%s","PEMa4Excl2_",d,".tif")
         writeRaster(PEMa4EXCL2,paste("PEMa4Excl2_",d,sep=""), format="GTiff", overwrite=TRUE)
         ConvPNG(rast_name)
		}
    }
  }
} #### if PEMS
setwd(prism_dir)
#### END Model Done - final map production ####

warnings()
#Possibility to calculate any number of outputs. This example was for 2014
#data only, but will want to look at multi-year calculations and how we
#can express uncertainty (annual variability) for more static risk maps.
