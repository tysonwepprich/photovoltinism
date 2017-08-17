 require(sp)
 require(rgdal)
 require(raster)


setwd("F:/PRISM/2014")
Con=function(condition, trueValue, falseValue){
    return(condition * trueValue + (!condition)*falseValue)
    }

#Daily PRISM Statistics (degree days and thresholds) for pest development models

#### The old NAPPFAST model parameters:
#Lower Base temp = 12 degrees C, Upper Base Temp = 40 degrees C  (tmean)
#growing degree days threshold = 450 DD
#Exclusion, cold temperature: 5+ days of tmean -3 degrees C or less

####Pest Specific, Multiple Life Stage Phenology Model Parameters:
#LDT = lower development threshold, temp at which growth = 0 (using PRISM tmean)
#DD = degree days, number of cumulative heat units to complete that lifestage
eggLDT = 11.93
larvaeLDT = 11.6
pupaeLDT = 11.9
adultLDT = 12.2 #for oviposition
eggDD = 69.3
larvaeDD = 156
pupaeDD = 174 #females
adultDD = 79.2 #time to 50% eggs laid
#Temperature exclusion threshold parameters
#LLT = lower lethal temperature (PRISM tmin), ULT = upper lethal temperature (PRISM tmax)
eggLLT = -3
eggULT = 41
larvaeLLT = -12
larvaeULT = 40
adultLLT = 0.5

#Search pattern for PRISM daily temperature grids. Load them for processing.
pattern = paste("(PRISM_tmean_)(.*)(_bil.bil)$", sep="")
files <- list.files(pattern=pattern, all.files=FALSE, full.names=TRUE)
#Check that there are enough files for the year
length(files)
numlist <- vector()
for (file in files) {
    num = strsplit(file,split="_")[[1]][5]
    numlist <-c(numlist,num)
    }
#Order by date starting with first downloaded date. Up to user to download date range
#in line with species' biofix. Assumes Jan 01  biofix
#Sorting is necessary in most recent year with provisional data and different filename structure
sortedlist <- sort(numlist)
filelist <- vector()
#Initialize Tracking Rasters
#Read in first raster as a template and intialize tracking rasters
template <- raster(files[1])
template[!is.na(template)] <- 0
#Initialize all tracking rasters as zero with the template
DDaccum <- template
Lifestage <- template
NumGen <- template
#rm(template)
#Lifestage: [0] = egg, [1] = larvae, [2] = pupae, [3] = adult

# Accumulate degree days and reclass cells to NA with temperature exclusion
# NA for exclusion means that DD accum cannot occur anymore with incomplete
# generation development and no oviposition to next generation. This may
# be tested for sensitivity by allowing exlusion cells to reset to zero instead.

###For testing only
## sublist <- sortedlist[190:220]
## for (d in sublist) {

for (d in sublist) {
    print(d)
    if (.Platform$OS.type == "windows") flush.console()
    Sys.sleep(1)
    #Read in that day's PRISM raster files
    pattern = paste("(PRISM_tmean_)(.*)(",d,")(_bil.bil)$", sep="")
    temp <- list.files(pattern=pattern,all.files=FALSE, full.names=TRUE)
    tmean <- raster(temp)
    pattern = paste("(PRISM_tmin_)(.*)(",d,")(_bil.bil)$", sep="")
    temp <- list.files(pattern=pattern,all.files=FALSE, full.names=TRUE)
    tmin <- raster(temp)
    pattern = paste("(PRISM_tmax_)(.*)(",d,")(_bil.bil)$", sep="")
    temp <- list.files(pattern=pattern,all.files=FALSE, full.names=TRUE)
    tmax <- raster(temp)
    rm(pattern,temp)

    #Create stage specific lifestage binary rasters
    #Limits operations to a mask for cells that are in that lifestage
    #This is what allows for pixel by pixel tracking of what lifestage
    # that cell is in

    LS0 <- Lifestage == 0
    LS1 <- Lifestage == 1
    LS2 <- Lifestage == 2
    LS3 <- Lifestage == 3
    writeRaster(Lifestage,paste("Lifestage_",d,sep=""), format="GTiff",overwrite=TRUE)

    for (i in 1:4) {
      if (i==1){
         #developmental degree days  (zero values for temps below LDT)
         dd0 <- ((tmean > eggLDT) * (tmean - eggLDT))
         #Calculate lower lethal threshold and exclusion mask
         eggmin <- tmin > eggLLT
         eggmin[eggmin==0] <- NA
         writeRaster(eggmin,paste("EggMin_",d,sep=""), format="GTiff",overwrite=TRUE)
         #Calculate upper lethal threshold and exclusion mask
         eggmax <- tmax < eggULT
         eggmax[eggmax==0] <- NA
         writeRaster(eggmax,paste("EggMax_",d,sep=""), format="GTiff",overwrite=TRUE)
         #Apply exclusions and lifestage mask to daily degree day surface
         dd0 <- dd0 * eggmin * eggmax * LS0

         #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
         #count cells with dd>0
         dd.stat <- cellStats(dd0,stat='max',na.rm=TRUE)
         if (dd.stat > 0) {
           DDaccum <- DDaccum + dd0
           #Calculate lifestage progression: for dd accum in correct lifestage, is it >= egg DD threshold?
           progress0 <- (DDaccum * LS0) >= eggDD
           writeRaster(progress0,paste("Progress0_",d,sep=""), format="GTiff",overwrite=TRUE)
           Lifestage <- Con(LS0 == 1 & progress0 == 1, 1, Lifestage)
           #Reset the DDaccum cells to zero for cells that progressed to next lifestage
           DDaccum <- Con(progress0 == 0, DDaccum, 0)
         }

      } else if (i == 2) {
         #developmental degree days
         dd1 <- ((tmean > larvaeLDT) * (tmean - larvaeLDT))
         #Calculate lower lethal threshold and exclusion mask
         larmin <- tmin > larvaeLLT
         larmin[larmin==0] <- NA
         writeRaster(larmin,paste("Larmin_",d,sep=""), format="GTiff",overwrite=TRUE)
         #Calculate upper lethal threshold and exclusion mask
         larmax <- tmax < larvaeULT
         larmax[larmax == 0] <- NA
         writeRaster(larmax,paste("Larmax_",d,sep=""), format="GTiff",overwrite=TRUE)
         #Apply exclusions and lifestage mask to daily degree day surface and limit to correct stage
         dd1 <- dd1 * larmin * larmax * LS1

         #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
         dd.stat <- cellStats(dd1,stat='max',na.rm=TRUE)
         if (dd.stat > 0) {
           DDaccum <- DDaccum + dd1
           #Calculate lifestage progression: for dd accum in correct lifestage, is it >= larvae DD threshold?
           progress1 <- (DDaccum * LS1) >= larvaeDD
           writeRaster(progress1,paste("Progress1_",d,sep=""), format="GTiff",overwrite=TRUE)
           Lifestage <- Con(LS1 == 1 & progress1 == 1, 2, Lifestage)
           #Reset the DDaccum cells to zero for cells that progressed to next lifestage
           DDaccum <- Con(progress1 == 0, DDaccum, 0)
         }

      } else if (i == 3) {
         #developmental degree days
         dd2 <- ((tmean > pupaeLDT) *(tmean - pupaeLDT))
         #Apply exclusions (none for pupae stage) and lifestage mask to daily degree day surface and limit to correct stage
         dd2 <- dd2 * LS2

         #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
         dd.stat <- cellStats(dd2,stat='max',na.rm=TRUE)
         if (dd.stat > 0) {
           DDaccum <- DDaccum + dd2
           #Calculate lifestage progression: for dd accum in correct lifestage, is it >= pupae DD threshold?
           progress2 <- (DDaccum * LS2) >= pupaeDD
           writeRaster(progress2,paste("Progress2_",d,sep=""), format="GTiff",overwrite=TRUE)
           Lifestage <- Con(LS2 == 1 & progress2 == 1, 3, Lifestage)
           #Reset the DDaccum cells to zero for cells that progressed to next lifestage
           DDaccum <- Con(progress2 == 0, DDaccum, 0)
         }

      } else { #adult stage, or time to 50% oviposition
         #developmental degree days
         dd3 <- ((tmean > adultLDT) * (tmean - adultLDT))
         #Calculate lower lethal threshold and exclusion mask
         admin <- tmin > adultLLT
         admin[admin==0] <- NA
         writeRaster(admin,paste("Admin_",d,sep=""), format="GTiff",overwrite=TRUE)
         #Apply exclusions and lifestage mask to daily degree day surface
         dd3 <- dd3 * admin * LS3

         #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
         dd.stat <- cellStats(dd3,stat='max',na.rm=TRUE)
         if (dd.stat > 0) {
           DDaccum <- DDaccum + dd3
           #Calculate lifestage progression: for dd accum in correct lifestage, is it >= adult DD threshold?
           progress3 <- (DDaccum * LS3) >= adultDD
           writeRaster(progress3,paste("Progress3_",d,sep=""), format="GTiff",overwrite=TRUE)
           #Reset the DDaccum cells to zero for cells that progressed to next lifestage
           DDaccum <- Con(progress3 == 1,0, DDaccum)
           #Remove masking effect from progress counter so it doesn't perpetuate through NumGen raster
           progress3[is.na(progress3)] <- template[is.na(progress3)]
           #Reset Lifestage raster with new generation. By replacing NA with template zero values, we remove generation dependency
           Lifestage <- Con(LS3 == 1 & progress3 == 1,0, Lifestage)
           Lifestage[is.na(Lifestage)] <- template[is.na(Lifestage)]


           #Increment NumGen + 1
           NumGen <- NumGen + progress3
           writeRaster(NumGen,paste("NumGen",d,sep=""), format="GTiff",overwrite=TRUE)
         }
      }
    }
}

writeRaster(NumGen,"FCM_NumGenerations.tif", format="GTiff",overwrite=TRUE)
#Possibility to calculate any number of outputs. This example was for 2014
#data only, but will want to look at multi-year calculations and how we
#can express uncertainty (annual variability) for more static risk maps.
