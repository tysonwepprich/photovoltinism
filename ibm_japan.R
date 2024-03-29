# Japan sites optimal CP

source("CDL_funcs.R")
library(GSODR)
library(dplyr)
library(ggplot2)
library(ggridges)


# look at other nearby stations to Inunaki pass (33.675, 130.55)
# Fritzi collected in Kumamoto prefecture at 747-838m elev (center ~ 32.73, 130.79)
# Northern population near Lake Toya, Hokkaido (center ~ 42.602851, 140.851949)
# New biotype for Europe from Murakami (38.2, 139.5)
tbar_stations <- nearest_stations(LAT = 38.25,
                                  LON = 140.35,
                                  distance = 15)
load(system.file("extdata", "isd_history.rda", package = "GSODR"))

stations <- isd_history %>% 
  mutate(enddate = lubridate::ymd(END),
         startdate = lubridate::ymd(BEGIN)) %>% 
  filter(STNID %in% tbar_stations,
         startdate < as.Date("1989-12-31"),
         enddate > as.Date("2015-12-31")) %>% 
  dplyr::pull(STNID)




# just get IIZUKA for now, maybe adjust temperature with lapse rate for "International Standard Atmosephere"
# nothing nearby Inunaki pass at similar elevation for comparison
years <- 1990:1991
# stations <- c("478030-99999", "478070-99999", "478080-99999", "478090-99999")
outlist <- list()
safe_GSOD <- purrr::safely(get_GSOD)
for (yr in years){
  for (st in stations){
    tbar <- safe_GSOD(years = yr, station = st)
    
    # raw temperature without elev adjustment, no focal site
    if(nrow(tbar$result) > 0){
      elev_stat <- tbar$result$ELEV_M[1]
      
      temp <- tbar$result %>% 
        dplyr::select(YDAY, MAX, MIN) %>% 
        arrange(YDAY) %>% 
        mutate(
          DD = TriDD(MAX, MIN, 6.9, 30),
          AccumDD = cumsum(DD),
          YEAR = yr,
          STNID = st,
          ELEV_M = elev_stat)
      outlist[[length(outlist)+1]] <- temp
    }else{
      next
    }
    
    # # uses lapse rates
    # elev_stat <- tbar$ELEV_M[1]
    # elev_site <- 415 
    # 
    # temp <- tbar %>% 
    #   dplyr::select(YDAY, MAX, MIN) %>% 
    #   arrange(YDAY) %>% 
    #   mutate(MAX_adj = MAX - 6.49 * (elev_site - elev_stat) / 1000,
    #          MIN_adj = MIN - 6.49 * (elev_site - elev_stat) / 1000,
    #          DD = TriDD(MAX, MIN, 6.1, 35),
    #          DD_adj = TriDD(MAX_adj, MIN_adj, 6.1, 35),
    #          AccumDD = cumsum(DD),
    #          AccumDDadj = cumsum(DD_adj),
    #          YEAR = yr,
    #          STNID = st)
    
    # outlist[[length(outlist)+1]] <- temp
  }
}
outdf <- bind_rows(outlist) %>% 
  # left_join(stations[, c("STNID", "ELEV_M")]) %>% 
  mutate(Site = "Murakami") %>% 
  group_by(YEAR, YDAY, Site) %>% 
  summarise(tmax = mean(MAX, na.rm = TRUE),
            tmin = mean(MIN, na.rm = TRUE),
            degday = TriDD(tmax, tmin, 10, 37.8))



saveRDS(outdf, "sgermany.rds")


df1 <- readRDS("sjapan.rds")
df2 <- readRDS("njapan.rds")

df <- bind_rows(df1, df2) %>% 
  filter(complete.cases(.)) %>% 
  mutate(date = as.Date(paste(YEAR, YDAY, sep = "-"), "%Y-%j")) %>% 
  mutate(MONTH = lubridate::month(date)) %>% 
  group_by(STNID, YDAY) %>% 
  summarise(tmax = mean(MAX),
            tmin = mean(MIN),
            avgdd = mean(DD),
            accumdd = mean(AccumDD),
            elevation = ELEV_M[1],
            region = REGION[1],
            mon = MONTH[1],
            station = paste(region, elevation, sep = ":"))

theme_set(theme_bw(base_size = 20)) 

# temperature ranges
plt <- ggplot(df, aes(x = YDAY, group = STNID, color = elevation)) +
  geom_linerange(aes(ymin = tmin, ymax = tmax)) +
  facet_wrap(~station, ncol = 2) +
  theme_bw()
plt

# mean by month
sumdf <- df %>% 
  group_by(station, mon) %>% 
  summarise(meantmax = mean(tmax),
            meantmin = mean(tmin),
            Month = month.abb[mon[1]])


plt <- ggplot(sumdf, aes(x = reorder(Month, mon), group = station, color = station)) +
  geom_linerange(aes(ymin = meantmin, ymax = meantmax), size = 2, position=position_dodge(width = .5)) +
  ylab("Average maximum and minimum temperatures") +
  xlab("Month") +
  ggtitle("Weather station temperatures by month\nin North and South Japan at different elevations")
plt


ggsave(filename = "japan_temprange.png", plot = plt,
       device = "png", width = 10, height = 5, units = "in")


# daily degree-days
springdf <- df

plt2 <- ggplot(data = springdf, aes(x = avgdd, y = mon, group = mon)) +
  geom_density_ridges(scale = .9) +
  ylab("Month") +
  xlab("Mean daily degree-days (6.1C LDT)") +
  ggtitle("Daily degree-days across year") +
  facet_wrap(~region, ncol = 2)
plt2

ggsave(filename = "japan_ddrange.png", plot = plt2,
       device = "png", width = 10, height = 8, units = "in")



# simple accum degree day graphs
plt <- ggplot(outdf, aes(x = YDAY, y = AccumDD, group = as.factor(YEAR), color = as.factor(YEAR))) +
  geom_line() +
  facet_wrap(~STNID, ncol = 2)
plt

# what to use based on April 15 emergence observations?
# will use DD adjusted for elevation with lapse rate
outdf %>% filter(YDAY == 105, STNID %in% tbar_stations[c(1,3,4)]) %>% 
  group_by(STNID) %>% 
  summarise(meanDD = mean(AccumDD, na.rm = TRUE), 
            meanDDadj = mean(AccumDDadj, na.rm = TRUE))


