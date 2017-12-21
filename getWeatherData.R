# trying out some weather data packages
# 
# # needs API
# library(rwunderground)
# 
# history_range(set_location(lat_long="33.675,130.55"),
#               date_start = "20140101", date_end = "20140531",
#               limit = 10, no_api = TRUE, use_metric = TRUE,
#               raw = FALSE, message = TRUE)

# Global Summary of the Day seems best

source("CDL_funcs.R")
library(GSODR)
library(dplyr)


# just get IIZUKA for now, maybe adjust temperature with lapse rate for "International Standard Atmosephere"
# nothing nearby Inunaki pass at similar elevation for comparison
years <- 2008:2017
stations <- c("478030-99999", "478070-99999", "478080-99999", "478090-99999")
outlist <- list()
for (yr in years){
  for (st in stations){
    tbar <- get_GSOD(years = yr, station = st)
    elev_stat <- tbar$ELEV_M[1]
    elev_site <- 415 
    
    temp <- tbar %>% 
      dplyr::select(YDAY, MAX, MIN) %>% 
      arrange(YDAY) %>% 
      mutate(MAX_adj = MAX - 6.49 * (elev_site - elev_stat) / 1000,
             MIN_adj = MIN - 6.49 * (elev_site - elev_stat) / 1000,
             DD = TriDD(MAX, MIN, 6.1, 35),
             DD_adj = TriDD(MAX_adj, MIN_adj, 6.1, 35),
             AccumDD = cumsum(DD),
             AccumDDadj = cumsum(DD_adj),
             YEAR = yr,
             STNID = st)
    outlist[[length(outlist)+1]] <- temp
  }
}
outdf <- bind_rows(outlist)

library(ggplot2)
plt <- ggplot(outdf, aes(x = YDAY, y = AccumDD, group = as.factor(YEAR), color = as.factor(YEAR))) +
  geom_line() +
  geom_line(aes(y = AccumDDadj), linetype = 2) +
  facet_wrap(~STNID, ncol = 2)
plt

# what to use based on April 15 emergence observations?
# will use DD adjusted for elevation with lapse rate
outdf %>% filter(YDAY == 105, STNID %in% tbar_stations[c(1,3,4)]) %>% 
  group_by(STNID) %>% 
  summarise(meanDD = mean(AccumDD, na.rm = TRUE), 
            meanDDadj = mean(AccumDDadj, na.rm = TRUE))


# look at other nearby stations to Inunaki pass
tbar_stations <- nearest_stations(LAT = 33.675,
                                  LON = 130.55,
                                  distance = 30)
load(system.file("extdata", "country_list.rda", package = "GSODR"))
load(system.file("extdata", "isd_history.rda", package = "GSODR"))

station_locations <- left_join(isd_history, country_list,
                               by = c("CTRY" = "FIPS"))

Oz <- filter(station_locations, COUNTRY_NAME == "JAPAN")

Oz %>% 
  filter(STNID %in% tbar_stations[c(1,3,4)]) %>%  data.frame()
