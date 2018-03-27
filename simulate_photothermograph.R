# Photothermograph simulation
# for eventual use in Shiny App

# This code will:
# 1. Download Daymet data for a site and year
# 2. Calculate accumulated degree-days from temperature data
# 3. Calculate photoperiod for each day
# 4. Simulate individual variation in critical photoperiod and/or required degree-days.
# 5. Plot photothermographs, optimal voltinism heatmap, voltinism histogram, diapause proportion

# Setup ----
library(ggplot2)
library(viridis)
library(dplyr)
# library(tidyr)
# library(lubridate)
library(purrr)
# library(stringr)
library(daymetr)
library(ggridges)

# Utility functions
source('CDL_funcs.R')

# User input species traits ----
nsim <- 100
lat <- 44.56400
lon <- -123.26300 
startyear <- 2005
endyear <- 2015
# TODO: 3 biocontrol species parameters set for you with biotypes
# custom species params
ldt <- 10
udt <- 37
cdl <- 15.5
emerg_dd <- 100
gen_dd <- 367 # egg to teneral adult
povip_dd <- 126 # pre-ovip period
ovip_dd <- 40 # guessing here, use poisson process
cdl_sd <- .5
emerg_dd_sd <- 15
gen_dd_sd <- 30
povip_dd_sd <- 10
lambda <- 3.5 # population growth rate function of multiple eggs laid per adult, plus random mortality
pars <- data.frame(ldt = ldt, udt = udt, cdl = cdl, emerg_dd = emerg_dd,
                   gen_dd = gen_dd, povip_dd = povip_dd, ovip_dd, cdl_sd = cdl_sd, 
                   emerg_dd_sd = emerg_dd_sd, gen_dd_sd = gen_dd_sd, 
                   povip_dd_sd = povip_dd_sd, lambda = lambda)


# Download weather data ----
temp <- download_daymet(site = "default", lat = lat, lon = lon,
                        start = startyear, end = endyear, internal = TRUE,
                        silent = TRUE, force = FALSE)

# Calculate degree days and photoperiod ----
gdd <- temp$data %>% 
  filter(yday < 350) %>% 
  rowwise() %>% 
  mutate(degday = TriDD(tmax..deg.c., tmin..deg.c., ldt, udt)) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  arrange(yday) %>% 
  mutate(accumdegday = cumsum(degday),
         daylength = photoperiod(lat, yday))



SimLifecycle <- function(pars, nsim, yr, gdd){
  gdd <- gdd %>% filter(year == yr)
  maxgdd <- max(gdd$accumdegday, na.rm = TRUE)
  firstfrost <- min(gdd$accumdegday[which(gdd$tmin..deg.c. <= 0 & gdd$yday > 200)])
  lastfrost <- max(gdd$accumdegday[which(gdd$tmin..deg.c. <= 0 & gdd$yday < 180)])
  
  dflist <- list()
  
  # Simulate overwinter emergence
  emerg_dd <- rnorm(nsim, mean = pars$emerg_dd, sd = pars$emerg_dd_sd)
  gen_dd <- rnorm(nsim, mean = pars$gen_dd, sd = pars$gen_dd_sd)
  povip_dd <- rnorm(nsim, mean = pars$povip_dd, sd = pars$povip_dd_sd)
  ovip_dd <- rpois(nsim, lambda = pars$ovip_dd)
  genlength1 <- emerg_dd + gen_dd
  genlength2 <- gen_dd + povip_dd + ovip_dd
  df <- data.frame(index = seq_len(nsim), sens_dd = emerg_dd, gen = 0, reproduce = TRUE) %>% 
    mutate(rowid = row_number()) %>% 
    group_by(rowid) %>% 
    mutate(chooserow = which.max(sens_dd <= gdd$accumdegday),
           accumdegday = gdd$accumdegday[chooserow],
           yday = gdd$yday[chooserow],
           daylength = gdd$daylength[chooserow],
           cdl = NA,
           diap_dd = genlength1[rowid],
           repro_dd = genlength2[rowid],
           fit_diapause_gdd = maxgdd >= (accumdegday + genlength1[rowid]),
           fit_repro_gdd = maxgdd >= (accumdegday + genlength2[rowid]),
           diapause_before_frost = firstfrost > (accumdegday + genlength1[rowid]),
           emerge_after_frost = lastfrost < sens_dd) %>% 
    ungroup() %>%
    select(-rowid)
  dflist[[1]] <- df
  
  # Maximal potential voltinism without photoperiod
  maxvolt <- vector(length = length(emerg_dd))
  lastsensgdd <- vector(length = length(emerg_dd))
  
  for (i in seq_along(emerg_dd)){
    maxvolt[i] <- length(seq(from = genlength1[i],
                          to = maxgdd, 
                          by = genlength2[i]))
    lastsensgdd[i] <- max(seq(from = genlength1[i],
                                      to = maxgdd, 
                                      by = genlength2[i]))
  }
  
  voltinism <- data.frame(index = seq_len(nsim), emerg_dd = emerg_dd, gen_dd = gen_dd, povip_dd = povip_dd,
                          ovip_dd = ovip_dd, maxvolt = maxvolt, maxgdd = maxgdd, firstfrost = firstfrost, lastfrost = lastfrost,
                          lastsensgdd = lastsensgdd)

  
  eggs <- ceiling(pars$lambda)
  mortality <- (eggs - pars$lambda) / eggs
  
  SimExpand <- function(dfrow){
    outdf <- dfrow[rep(seq_len(nrow(dfrow)), eggs), ]
    outdf$eggs <- seq_len(eggs)
    outdf$sens_dd <- Cond(outdf$gen[1] == 0, outdf$accumdegday + outdf$diap_dd,
                            outdf$accumdegday + outdf$repro_dd)
    outdf <- outdf %>% 
      rowwise() %>% 
      mutate(chooserow = which.max(sens_dd <= gdd$accumdegday),
             accumdegday = ifelse(chooserow == 1, NA, gdd$accumdegday[chooserow]),
             yday = ifelse(chooserow == 1, NA, gdd$yday[chooserow]),
             daylength = ifelse(chooserow == 1, NA, gdd$daylength[chooserow]),
             gen = gen + 1,
             cdl = rnorm(1, mean = pars$cdl, sd = pars$cdl_sd),
             reproduce = daylength >= cdl,
             fit_diapause_gdd = maxgdd >= (accumdegday + diap_dd),
             fit_repro_gdd = maxgdd >= (accumdegday + repro_dd),
             diapause_before_frost = firstfrost > (accumdegday + diap_dd),
             emerge_after_frost = NA)
             
    return(outdf)
  }
  
  while(length(which(df$reproduce == TRUE & df$fit_repro_gdd == TRUE)) > 0){
    df <- df %>% 
      ungroup() %>% 
      filter(reproduce == TRUE & fit_repro_gdd == TRUE) %>% 
      group_by(row_number()) %>% 
      do(SimExpand(.)) %>% 
      ungroup() %>% 
      sample_frac(size = (1 - mortality), replace = FALSE)
    dflist[[length(dflist)+1]] <- df
  }

  simresults <- bind_rows(dflist) %>% 
    left_join(voltinism, by = "index") %>% 
    select(-`row_number()`, -eggs, -chooserow)

  return(simresults)  
}

  
# Simulate ----
sims <- expand.grid(year = unique(gdd$year))
results <- sims %>% 
  group_by(year) %>% 
  do(SimLifecycle(pars, nsim, .$year, gdd))


# Marginal histograms ----
mean_cdl <- mean(results$cdl, na.rm = TRUE)
maxgdd <- max(gdd$accumdegday, na.rm = TRUE)
maxphoto <- max(c(max(gdd$daylength), max(results$cdl, na.rm = TRUE)))
minphoto <- min(gdd$daylength)

cdl_hist <- results %>% 
  ungroup() %>% 
  select(cdl) %>%
  filter(complete.cases(.)) %>% 
  mutate(breaks = cut(cdl, breaks=seq(min(cdl),max(cdl),length.out = 20), 
                      labels=seq(min(cdl),max(cdl), length.out = 20)[-1], 
                      include.lowest=TRUE)) %>% 
  mutate(daylength = as.numeric(as.character(breaks))) %>%
  group_by(daylength) %>% 
  summarise(n = n()) %>% 
  mutate(accumdegday = -50 + 50 * n/max(n))

sens_hist <- results %>% 
  ungroup() %>% 
  select(accumdegday) %>%
  filter(complete.cases(.)) %>% 
  mutate(breaks = cut(accumdegday, breaks = seq(min(accumdegday), max(accumdegday), length.out = 100),
                      labels = seq(min(accumdegday), max(accumdegday), length.out = 100)[-1],
                      include.lowest = TRUE)) %>% 
  mutate(accumdegday = as.numeric(as.character(breaks))) %>% 
  group_by(accumdegday) %>% 
  summarise(n = n()) %>% 
  mutate(daylength = minphoto - .5 + .5 * n/max(n))

# Photothermograph ----
theme_set(theme_bw(base_size = 16)) 

plt <- ggplot(gdd, aes(x = accumdegday, y = daylength, group = as.factor(year), color = as.factor(year))) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(-50, maxgdd + 100), ylim = c(minphoto - .5, maxphoto), expand = FALSE) +
  scale_color_viridis(discrete = TRUE, name = "Year", option = "D") +
  guides(color = guide_legend(reverse=FALSE)) +
  geom_segment(data=cdl_hist, size=2.5, show.legend=FALSE,
               aes(x=-50, xend=accumdegday, y=daylength, yend=daylength), inherit.aes = FALSE) +
  geom_segment(data=sens_hist, size=2.5, show.legend=FALSE,
               aes(x=accumdegday, xend=accumdegday, y=minphoto - .5, yend=daylength), inherit.aes = FALSE) +
  geom_hline(yintercept = mean_cdl, linetype = "dashed", size = 1, alpha = .5) +
  # geom_rug(data = results, aes(x = accumdegday), color = "gray", alpha = .1) +
  # geom_rug(data = results, aes(y = cdl), color = "gray", alpha = .1) +
  ggtitle("Photothermographs with critical photoperiod\nand simulated sensitive stage emergence") +
  xlab("Accumulated degree-days") +
  ylab("Daylength (hours)")
plt


# Voltinism heatmap ----
# NEW ISSUE: these values are from nsim estimates, not including eggs
# so doesn't match up with voltinism cost plots below
# also mortality makes attempted generation wonky
res1 <- results %>% 
  select(year, index, maxvolt, gen) %>% 
  group_by(year, index) %>% 
  summarise(maxvolt = maxvolt[1],
            attemptvolt = max(gen, na.rm = TRUE))
plt1 <- ggplot(res1, aes(x = maxvolt, y = attemptvolt)) +
  # geom_tile(aes(fill = ..count..)) +
  stat_bin2d(aes(fill=..density..), geom="raster", hjust = 0.5, vjust = 0.5, position="identity") +
  scale_fill_viridis() +
  coord_cartesian(xlim = c(0, max(results$maxvolt) + 1), 
                  ylim = c(0, max(results$maxvolt) + 1), expand = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
plt1


# Voltinism cost ----
conseq <- results %>%
  mutate(consequence = ifelse(reproduce == TRUE & fit_diapause_gdd == TRUE & diapause_before_frost == TRUE,
                              "reproduce, F1 fits",
                              # ifelse(reproduce == TRUE & fit_diapause_gdd == TRUE & diapause_before_frost == FALSE,
                              #        "reproduce, F1 dies from frost",
                                     ifelse(reproduce == TRUE & (fit_diapause_gdd == FALSE | diapause_before_frost == FALSE),
                                            "reproduce, F1 doesn't fit",
                                            "diapause")),
         date = as.Date(yday, origin=as.Date("2015-12-31")))




plt2 <- ggplot(conseq, aes(x = accumdegday, y = consequence, 
                            fill = consequence, height = ..count..)) +
  geom_density_ridges(scale = .9, rel_min_height=0.01, stat = "density") +
  # geom_density_ridges(scale = .95, rel_min_height=0.001, stat = "binline", bins = 100) +
  scale_fill_viridis(discrete = TRUE, alpha = .5) +
  scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0.1, .1)) 
  # facet_wrap(~year, nrow = 2)
plt2


plt3 <- ggplot(conseq, aes(x = date, y = consequence, 
                            fill = consequence, height = ..count..)) +
  # geom_density_ridges(scale = .9, rel_min_height=0.01, stat = "density") +
  geom_density_ridges(scale = .95, rel_min_height=0.001, stat = "binline", bins = 50) +
  scale_fill_viridis(discrete = TRUE, alpha = .5) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  # scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0.1, .1)) 
  # facet_wrap(~year, nrow = 2)
plt3





