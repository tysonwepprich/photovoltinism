# Take Gerricke's FCM code and convert to Galerucella parameters
# Test out photoperiod and substages with constrained extext (Oregon)


library(sp)
library(rgdal)
library(raster)

# for parallel simulations with control over seed for reproducibility
# library(doRNG)


# setwd("F:/PRISM/2014") # remove this
prism_path <- "/data/PRISM/2014"

source('CDL_funcs.R')
#Daily PRISM Statistics (degree days and thresholds) for pest development models


####Pest Specific, Multiple Life Stage Phenology Model Parameters:
#LDT = lower development threshold, temp at which growth = 0 (using PRISM tmean)
#DD = degree days, number of cumulative heat units to complete that lifestage
start_doy  <- 100
end_doy    <- 300
stgorder   <- c("OA","E","L","P","A","F")
photo_sens <- 3 #c(-1, 3) # integer life stages for now
CDL        <- 14.25
model_CDL  <- 1 # if 1, model photoperiod decision
owstage    <- "OA"
OWadultDD  <- 108 # text OW stage dev 39 DD "post diapause"
calctype   <-"triangle"
eggLDT     <- 10
eggUDT     <- 35  #
larvaeLDT  <- 10
larvaeUDT  <- 35  #upper dev. threshold-need to verify
pupaeLDT   <- 10
pupaeUDT   <- 35
adultLDT   <- 13.5 #for oviposition
adultUDT   <- 35
eggDD = 93.3
larvaeDD = 136.4 # 46.9 + 45.8 + 43.7 instars
pupDD = 137.7 
adultDD = 125.9 #time to oviposition
OWadultDD = 100
region_param <- "OR"

#Search pattern for PRISM daily temperature grids. Load them for processing.
pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
files <- list.files(path = prism_path, pattern=pattern, all.files=FALSE, full.names=TRUE)
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
REGION <- switch(region_param,
                 "CONUS"        = extent(-125.0,-66.5,24.0,50.0),
                 "NORTHWEST"    = extent(-125.1,-103.8,40.6,49.2),
                 "OR"           = extent(-124.7294, -116.2949, 41.7150, 46.4612))
template <- crop(raster(files[1]),REGION)

template[!is.na(template)] <- 0
#Initialize all tracking rasters as zero with the template
DDaccum <- template
ddtotal <- template
Lifestage <- template
Lifestage <- Lifestage -1  # Need for OW designated as stage -1
NumGen <- template
#rm(template)
#Lifestage: [0] = egg, [1] = larvae, [2] = pupae, [3] = adult

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
# Accumulate degree days and reclass cells to NA with temperature exclusion
# NA for exclusion means that DD accum cannot occur anymore with incomplete
# generation development and no oviposition to next generation. This may
# be tested for sensitivity by allowing exlusion cells to reset to zero instead.

###For testing only
sublist <- sortedlist[start_doy:end_doy]
## for (d in sublist) {

for (d in sublist) {
  print(d)
  if (.Platform$OS.type == "windows") flush.console()
  Sys.sleep(1)
  #Read in that day's PRISM raster files
  # pattern = paste("(PRISM_tmean_)(.*)(",d,")(_bil.bil)$", sep="")
  # temp <- list.files(pattern=pattern,all.files=FALSE, full.names=TRUE)
  # tmean <- raster(temp)
  pattern = paste("(PRISM_tmin_)(.*)(",d,")(_bil.bil)$", sep="")
  temp <- list.files(path = prism_path, pattern=pattern,all.files=FALSE, full.names=TRUE)
  tmin <- crop(raster(temp), REGION)
  pattern = paste("(PRISM_tmax_)(.*)(",d,")(_bil.bil)$", sep="")
  temp <- list.files(path = prism_path, pattern=pattern,all.files=FALSE, full.names=TRUE)
  tmax <- crop(raster(temp), REGION)
  rm(pattern,temp)
  
  # tmean not in GRUB PRISM data, use approximation
  tmean <- (tmax + tmin) / 2

  # photoperiod for this day across raster
  if (model_CDL == 1){
  doy <- lubridate::yday(lubridate::ymd(d))
  photo <- RasterPhoto(tmin, doy)
  }
  
  #### Loop through Stages; order of stages now read from SPP param file ####
  
  for (i in stgorder) {  # Handle stages in the model
    
    
    ####  MAIN STEPS FOR EGG STAGE ####
    if (i == "E" | i == "OE") {   # Egg Stage
      if(calctype=="average") { #devel DDs (zero values for temps below LDT)
        dd0tmp <- AvgDD(tmax,tmin,eggLDT,eggUDT)
      } else if(calctype=="triangle") {
        dd0tmp <- TriDD(tmax,tmin,eggLDT,eggUDT)
      } else { # assume (calctype=="simple") 
        dd0tmp <- SimpDD(tmean,eggLDT)
      }
      
      if (i == "OE") { dd0 <- dd0tmp * LSOW0 
      } else if (i == "E") { 
        dd0 <- dd0tmp * LS0 
      }
      
      
      #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
      #count cells with dd>0
      dd.stat <- cellStats(dd0,stat='max',na.rm=TRUE)
      if (dd.stat > 0) {
        DDaccum <- DDaccum + dd0
        # if (pems) {
        #   # set DOY in PEM if curr gen, not set yet, DDs > eventDDs, and Stage is eggs, else keep same value
        #   if (PEMnumgens > 0) {
        #     PEMe1 <- Cond(NumGen==1 & PEMe1 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe1) 
        #   }
        #   if (PEMnumgens > 1) {
        #     PEMe2 <- Cond(NumGen==2 & PEMe2 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe2) 
        #   }
        #   if (PEMnumgens > 2) {
        #     PEMe3 <- Cond(NumGen==3 & PEMe3 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe3) 
        #   }
        #   if (PEMnumgens > 3) {
        #     PEMe4 <- Cond(NumGen==4 & PEMe4 == 0 & DDaccum >= eggEvent, doy * LS0, PEMe4) 
        #   }
        # }
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
      
      #else { # just use lifestage mask
      if (i == "OL") { dd1 <- dd1tmp * LSOW1 }
      else if (i == "L") { dd1 <- dd1tmp * LS1 }
      #}
      
      #Accumulate degree days, if dd1 > 0 otherwise exclusion masks get applied to this generation.
      dd.stat <- cellStats(dd1,stat='max',na.rm=TRUE)
      if (dd.stat > 0) {
        DDaccum <- DDaccum + dd1
        # if (pems & larvaeEventDD) {
        #   # set DOY in PEM if curr. gen, not set yet, DDs > eventDDs, Stage is larvae, else keep same value
        #   if (PEMnumgens > 0) {
        #     PEMl1 <- Cond(NumGen==1 & PEMl1 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl1) 
        #   }
        #   if (PEMnumgens > 1) {
        #     PEMl2 <- Cond(NumGen==2 & PEMl2 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl2) 
        #   }
        #   if (PEMnumgens > 2) {
        #     PEMl3 <- Cond(NumGen==3 & PEMl3 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl3) 
        #   }
        #   if (PEMnumgens > 3) {
        #     PEMl4 <- Cond(NumGen==4 & PEMl4 == 0 & DDaccum >= larvaeEvent, doy * LS1, PEMl4) 
        #   }
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
      
      if (i == "OP") { dd2 <- dd2tmp * LSOW2
      } else if (i == "P") { dd2 <- dd2tmp * LS2
      }
      
      #Accumulate degree days, if dd0 > 0 otherwise exclusion masks get applied to this generation.
      dd.stat <- cellStats(dd2,stat='max',na.rm=TRUE)
      if (dd.stat > 0) {
        DDaccum <- DDaccum + dd2
        # if (pems & pupaeEventDD) {
        #   # set DOY in PEM if curr. gen, not set yet, DDs > eventDDs, Stage is pupa, else keep same value
        #   if (PEMnumgens > 0) {
        #     PEMp1 <- Cond(NumGen==1 & PEMp1 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp1) 
        #   }
        #   if (PEMnumgens > 1) {
        #     PEMp2 <- Cond(NumGen==2 & PEMp2 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp2) 
        #   }
        #   if (PEMnumgens > 2) {
        #     PEMp3 <- Cond(NumGen==3 & PEMp3 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp3) 
        #   }
        #   if (PEMnumgens > 3) {
        #     PEMp4 <- Cond(NumGen==4 & PEMp4 == 0 & DDaccum >= pupaeEvent, doy * LS2, PEMp4) 
        #   }
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
      if (dd.stat > 0) {
        DDaccum <- DDaccum + dd3
        # if(pems & adultEventDD) { # record DOYs in PEMS for adult events
        #   if (PEMnumgens > 0) {
        #     PEMa1 <- Cond(NumGen==1 & PEMa1 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa1) 
        #   }
        #   if (PEMnumgens > 1) {
        #     PEMa2 <- Cond(NumGen==2 & PEMa2 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa2) 
        #   }
        #   if (PEMnumgens > 2) {
        #     PEMa3 <- Cond(NumGen==3 & PEMa3 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa3) 
        #   }
        #   if (PEMnumgens > 3) {
        #     PEMa4 <- Cond(NumGen==4 & PEMa4 == 0 & DDaccum >= adultEvent, doy * LS3, PEMa4) 
        #   }
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
      
      
      ####  MAIN STEPS FOR END OF EACH DAY ####
    } else if (i == "F") { # end of the day placeholder 
      
      if (model_CDL == 1){
        # diapause decision
        # add one more integer to Lifestage, 4 == diapause
        sens_mask <- Cond(Lifestage %in% photo_sens, 1, 0)
        diap_mask <- Cond(photo < CDL, sens_mask, 0)
        Lifestage <- Cond(diap_mask == 1, 4, Lifestage)
        LS4 <- Lifestage == 4
      }
      
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
      
      
      # # new Nov 2015 tally stages plus gens
      # StageCount <- Lifestage + (NumGen * 4)
      
    }  # lifestage F or 5 (end of day calcs)
  }  # lifestage for loop
}  # daily loop

# writeRaster(NumGen,"FCM_NumGenerations.tif", format="GTiff",overwrite=TRUE)
plot(NumGen)

#Possibility to calculate any number of outputs. This example was for 2014
#data only, but will want to look at multi-year calculations and how we
#can express uncertainty (annual variability) for more static risk maps.
