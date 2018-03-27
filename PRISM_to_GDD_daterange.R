# This code will:
# 1. Download PRISM data if needed (unfortunately, does entire CONUS for each year)
# 2. Calculate degree-days for each years at multiple sites
# 3. Forecast future events based on average temperature in previous years
# 4. Plot date range when degree-day accumulation has occurred or is forecasted

# Note, I have 2018 PRISM data with 10-year-average forecasts from Len's server (so not reproducible)
# 5-year-average forecasts are made with the code below.

library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(lubridate)
library(sp)
library(rgdal)
library(raster)
library(purrr)
library(stringr)
library(prism)

# Input development thresholds
LDT <- 10
UDT <- 37

# Input where your PRISM downloaded data are
# base_path <- "/data/PRISM/"
base_path <- "prismDL"

# Input years to analyze
years <- c(2013:2018)

# Input sites with ID, x, y columns
sites <- read.csv("data/GCA_modeling_sites.csv", header = TRUE)


# Utility functions

# Take .bil files from PRISM yearly directories
# Return best data for each day
# Remove leap year if not needed
ExtractBestPRISM <- function(prismfiles, yr){
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
  
  # still has issue with leap year
  if (yr %% 4 != 0) {
    df2 <- df2[!grepl(pattern = paste0(yr, "0229"), x = df2$dates), ]
  }
  
  best <- prismfiles[df2$rownum]
  dates <- regexpr(pattern = "[0-9]{8}", text = best)
  
  fileorder <- order(regmatches(best, dates))
  prismfiles <- best[fileorder]
  return(prismfiles)
}

# Ifelse function used for rasters/vectors
Cond=function(condition, trueValue, falseValue){
  return(condition * trueValue + (!condition)*falseValue)
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



# PRISM download
# using ropensci 'prism' package to access webservice
# downloads entire CONUS, so files are large

startdate <- "2018-02-28"
enddate <- "2018-03-14"
yr <- year(ymd(startdate))

# need to set download directory
options(prism.path = paste("prismDL", yr, sep = "/"))

get_prism_dailys("tmin", minDate = startdate, 
                 maxDate = enddate, keepZip = FALSE)
get_prism_dailys("tmax", minDate = startdate, 
                 maxDate = enddate, keepZip = FALSE)



# Site based growing degree-days from PRISM
# Also can do this on rasters to make maps
outlist <- list()
for (yr in years){
  prism_path <- paste(base_path, yr, sep = "/")
  
  # Search pattern for PRISM daily temperature grids. Load them for processing.
  # TMIN
  pattern = paste("(PRISM_tmin_)(.*)(_bil.bil)$", sep="") # changed this to min, mean not on GRUB?
  tminfiles <- list.files(path = prism_path, pattern=pattern, 
                          all.files=FALSE, full.names=TRUE, recursive = TRUE)
  
  # PRISM has different names for files with different quality
  # this function selects according to hierarchy
  tminfiles <- ExtractBestPRISM(tminfiles, yr)
  
  r <- raster(tminfiles[1])
  tminstack <- stack(r)
  tminstack@layers <- sapply(tminfiles, function(x) { r@file@name=x; r } ) 
  
  # TMAX
  pattern = paste("(PRISM_tmax_)(.*)(_bil.bil)$", sep="")
  tmaxfiles <- list.files(path = prism_path, pattern=pattern, 
                          all.files=FALSE, full.names=TRUE, recursive = TRUE)
  tmaxfiles <- ExtractBestPRISM(tmaxfiles, yr)
  
  r <- raster(tmaxfiles[1])
  tmaxstack <- stack(r)
  tmaxstack@layers <- sapply(tmaxfiles, function(x) { r@file@name=x; r } ) 
  
  # by site
  tmin <- extract(tminstack, sites[, c("x", "y")])
  tmax <- extract(tmaxstack, sites[, c("x", "y")])
  
  GDD <- data.frame(TriDD(tmax, tmin, LDT, UDT))
  names(GDD) <- paste("yday", 1:ncol(GDD), sep = "_")
  gdd <- sites %>% 
    bind_cols(GDD) %>% 
    tidyr::gather(key = yday, value = DD, -ID, -x, -y) %>% 
    mutate(yday = unlist(map(stringr::str_split(string = yday, pattern = coll("_"), 2), 2)),
           year = yr)
  outlist[[length(outlist)+1]] <- gdd 
}

gdddf <- bind_rows(outlist)

saveRDS(gdddf, "gdddf.rds")

gdddf <- readRDS("gdddf.rds")


# Map range of GDD at different sites, 100-150DD dates
cutoffday <- 77 # day of year when PRISM observations switch to forecasts

# Add new predictions for current year
preddf <- gdddf %>% 
  filter(complete.cases(.),
         year != 2018) %>% 
  mutate(yday = as.numeric(yday)) %>% 
  arrange(yday) %>% 
  group_by(ID, yday) %>%   
  summarise(DD = mean(DD)) %>% 
  filter(yday > cutoffday) %>% 
  mutate(year = "2018+5yrmean")

# replace predictions after cutoff day
obs <- gdddf %>% 
  dplyr::select(-x, -y) %>% 
  mutate(yday = as.numeric(yday)) %>% 
  filter(complete.cases(.),
         year == 2018, 
         yday <= cutoffday) %>% 
  mutate(year = "2018+5yrmean")

preds <- bind_rows(preddf, obs)

gdddf$year <- as.factor(as.character(gdddf$year))
levels(gdddf$year)[6] <- "2018+10yrmean"

resdf <- gdddf %>% 
  dplyr::select(-x, -y) %>% 
  filter(complete.cases(.)) %>% 
  mutate(yday = as.numeric(yday)) %>% 
  bind_rows(preds) %>% 
  arrange(yday) %>% 
  group_by(year, ID) %>% 
  mutate(AccumDD = cumsum(DD)) %>% 
  filter(AccumDD >= 100) %>% 
  filter(AccumDD <= 150) %>% 
  filter(row_number() == 1L | row_number() == n()) %>% 
  mutate(date = as.Date(yday, origin=as.Date("2000-12-31")),
         bound = c("start", "end")) %>% 
  dplyr::select(ID, year, date, bound) %>% 
  tidyr::spread(bound, date) %>% 
  filter(ID %in% levels(gdddf$ID)[c(1, 3, 4, 5, 10, 12, 15, 16, 17, 21, 22, 25, 26, 28, 30)]) %>%  # East and West
  # filter(ID %in% levels(gdddf$ID)[c(1, 3, 5, 10, 15, 17, 21, 22, 25, 26, 30)]) %>%  # Just NW sites
  group_by(ID) %>% 
  mutate(meanstart = mean(start))




# Plot results


theme_set(theme_bw(base_size = 16)) 


plt <- ggplot(resdf, aes(x = reorder(ID, meanstart), ymin = start, ymax = end, group = as.factor(year), color = as.factor(year))) +
  # geom_jitter(size = 3, width = .15) +
  geom_linerange(position = position_dodge(width = .6), size = 2) +
  scale_color_viridis(discrete = TRUE, name = "Year", option = "D") +
  scale_y_date(date_breaks = "1 week", date_labels = "%b-%d") +
  coord_flip() +
  ggtitle("Date range for predicted Galerucella\noverwintering adult emergence (100-150 degree-days)") +
  xlab("Sites ordered by\nmean start of emergence") +
  guides(color = guide_legend(reverse=TRUE)) 

plt

ggsave(filename = "Galerucella_daterange3.png", plot = plt,
       device = "png", width = 16, height = 10, units = "in")


# output data to csv for import into Google Calendar
outdf <- resdf %>% 
  filter(year == "2018+5yrmean") %>% 
  dplyr::select(ID, start, end) %>% 
  mutate(start = start + duration(17, units = "years"),
         end = end + duration(17, units = "years"))
names(outdf) <- c("Subject", "Start Date", "End Date")
write.csv(outdf, file = "GoogleEvents.csv", row.names = FALSE)



