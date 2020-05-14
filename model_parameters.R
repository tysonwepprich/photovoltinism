# Parameters for optimal CP paper

source("CDL_funcs.R")
library(tidyverse)
library(lubridate)

theme_set(theme_bw(base_size = 18) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())) 
# Diorhabda ----

# Lifestage development ----
# Assuming common threshold

df <- read.csv("data/DCA_herrera.csv", header = TRUE)

ls <- df %>% 
  filter(stage != "total") %>% 
  group_by(temp) %>% 
  summarise(totdays = sum(days),
            devrate = 1 / sum(days))

mod <- lm(devrate ~ temp, data = ls)
(totdd <- 1/coef(mod)[2])
(-coef(mod)[1]/coef(mod)[2])
plot(ls$temp, ls$devrate)
abline(mod)

# Pre-oviposition from Lewis et al 2003
3.9 * (24.1 - 12.0)
# oviposition from Lewis et al 2003
12 * (24.1 - 12.0)


# percent in each stage
perc <- df %>% 
  # filter(temp > 15) %>% 
  filter(stage != "total") %>% 
  group_by(temp) %>% 
  mutate(percdays = days / sum(days))
percmean <- perc %>% 
  group_by(stage) %>% 
  summarise(mean = mean(percdays) * totdd)



ptdf <- df %>%
  filter(stage == "total") %>% 
  group_by(temp) %>%
  do(simdata = rnorm(n = .$n, mean = .$days, sd = .$se * sqrt(.$n))) %>%
  unnest()

ptdf %>% group_by(temp) %>% 
  summarise(sd = sd(x = simdata * (temp - 12.2)), 
            dd = mean(simdata * (temp - 12.2)))

ggplot(ptdf, aes(x = simdata)) +
  geom_density() +
  facet_wrap(~ temp)


# Diapause survival ----
# survival rate over days from Lewis et al. 2003
# Day 0.5 is July 19 == 200 yday
dat <- read.csv("data/DCA_summer_survival.csv", header = TRUE) %>% 
  arrange(Days) %>% 
  mutate(daydiff = c(NA, diff(Days)),
         survdiff = c(NA, diff(Survival)),
         survperc = Survival / 100,
         ydays = round(Days - .5) + 200)

datcut <- dat %>% filter(Days > 8)

mod1 <- lm(survperc ~ ydays, data = datcut)
datcut$pred1 <- predict(mod1, datcut)
plot(datcut$Days, datcut$pred1)
summary(mod1)

# what about survival rate vs degree-days?
library(GSODR)
tbar_stations <- nearest_stations(LAT = 31.057684, 
                                  LON = -97.346432,
                                  distance = 20)

tbar <- get_GSOD(years = 1999, station = tbar_stations) %>% 
  dplyr::select(YDAY, MAX, MIN) %>% 
  mutate(DD = TriDD(tmax = MAX, tmin = MIN, LDT = 12.0, UDT = 40)) %>% 
  arrange(YDAY) %>% 
  mutate(AccumDD = cumsum(DD))

datcut <- left_join(datcut, tbar, by = c("ydays" = "YDAY"))

# datcut <- dat %>% filter(Days > 10)
mod2 <- lm(survperc ~ AccumDD, data = datcut)
datcut$pred2 <- predict(mod2, datcut)
summary(mod2)


newdat <- tbar %>% filter(YDAY > 206) %>% 
  rename(ydays = YDAY) %>% 
  mutate(pred1 = predict(mod1, .),
         pred2 = predict(mod2, .))


dat <- left_join(dat, tbar, by = c("ydays" = "YDAY"))

ggplot(dat, aes(x = ydays, y = survperc)) +
  geom_point() +
  scale_y_continuous(limits = c(0, 1)) +
  geom_line(data = newdat, aes(y = pred1), color = "red") +
  geom_line(data = newdat, aes(y = pred2), color = "blue") +
  annotate("text", x = 350, y = .15, label = "Survival per day", color = "red") +
  annotate("text", x = 350, y = .4, label = "Survival per degree-day", color = "blue") +
  annotate("text", x = 210, y = .95, label = "???", color = "black") +
  xlab("Day of year") +
  ylab("Percent surviving")

ggplot(dat, aes(x = AccumDD, y = survperc)) +
  geom_point() +
  scale_y_continuous(limits = c(0, 1)) +
  geom_line(data = newdat, aes(y = pred1), color = "red") +
  geom_line(data = newdat, aes(y = pred2), color = "blue") +
  annotate("text", x = 3500, y = .25, label = "Survival per day", color = "red") +
  annotate("text", x = 3500, y = .55, label = "Survival per degree-day", color = "blue") +
  annotate("text", x = 2000, y = .95, label = "???", color = "black") +
  xlab("Accumulated degree-days") +
  ylab("Percent surviving")


# Emergence ----
library(daymetr)
tbar_stations <- nearest_stations(LAT = 39.9, 
                                  LON = -105.12,
                                  distance = 30)

tbar <- download_daymet(site = "GJ", lat = 39.9, lon = -105.12, 
                        start = 2012, end = 2016)
tbar <- tbar$data %>% 
  mutate(DD11 = TriDD(tmax = tmax..deg.c., tmin = tmin..deg.c., LDT = 11.1, UDT = 40)) %>% 
  mutate(DD12 = TriDD(tmax = tmax..deg.c., tmin = tmin..deg.c., LDT = 12.0, UDT = 40)) %>% 
  group_by(year) %>% 
  arrange(yday) %>% 
  mutate(AccumDD11 = cumsum(DD11),
         AccumDD12 = cumsum(DD12))

obs <- data.frame(year = c(2012:2016), 
                  yday = lubridate::yday(c("2012-04-19", "2013-04-30", "2014-05-02", "2015-04-14", "2016-04-26")))

obs <- left_join(obs, tbar)

library(fGarch)
x = seq(50, 350, length.out = 1000)
y = dt.ls(x, loc = arg1$mu, sigma2 = arg1$sigma2, shape = arg1$shape, nu = arg1$nu)
plot(x, y)

hist(rsnorm(1000, mean = 100, sd = 25, xi = 2), breaks = 100)


# Galerucella ----

# Lifestage development degree-days ----
# From McAvoy & Kok

# Pre-ovip
df <- data.frame(temp = c(12.5, 15, 20, 25, 27.5),
                 n = c(19, 11, 19, 15, 21),
                 days = c(82.4, 68.4, 13.4, 4.7, 3.8),
                 se = c(11.7, 18.5, 2.4, .2, .2))

# L1-L3 Table 3
df <- data.frame(temp = c(12.5, 15, 20, 25, 27.5),
                 n = c(125, 38, 60, 36, 72),
                 days = c(28.05, 25.44, 15.76, 10.9, 8.23),
                 se = c(0.28, 0.68, .17, .24, .12))

# Pupa
df <- data.frame(temp = c(12.5, 15, 20, 25, 27.5),
                 days = c(38.1, 29.3, 16.95, 9.81, 8.01))

# # Doesn't work with high variances
# ptdf <- df %>% 
#   group_by(temp) %>% 
#   do(simdata = rnorm(n = .$n, mean = .$days, sd = .$se * sqrt(.$n))) %>% 
#   unnest()


mod1 <- lm(1/simdata ~ temp, data = ptdf)
mod2 <- lm(1/days ~ temp, data = df)

(1/coef(mod1)[2])
(1/coef(mod2)[2])
(-coef(mod1)[1]/coef(mod1)[2])
(-coef(mod2)[1]/coef(mod2)[2])
plot(df$temp, 1/df$days)
abline(mod2)

# Assuming common threshold

df <- read.csv("data/GCA_mcavoy.csv", header = TRUE)

ls <- df %>% 
  filter(temp > 13) %>% 
  filter(stage != "ovip") %>% 
  group_by(temp) %>% 
  summarise(totdays = sum(days),
    devrate = 1 / sum(days))

mod <- lm(devrate ~ temp, data = ls)
(totdd <- 1/coef(mod)[2])
(-coef(mod)[1]/coef(mod)[2])
plot(ls$temp, ls$devrate)
abline(mod)

# percent in each stage
perc <- df %>% 
  filter(temp > 15) %>% 
  filter(stage != "ovip") %>% 
  group_by(temp) %>% 
  mutate(percdays = days / sum(days))
percmean <- perc %>% 
  group_by(stage) %>% 
  summarise(mean = mean(percdays) * totdd)

# What's the variance? SD of degree-day mean is about 5-10%
# Use L1 to L3 data at best temperatures
df <- data.frame(temp = c(12.5, 15, 20, 25, 27.5),
                 n = c(125, 38, 60, 36, 72),
                 days = c(28.05, 25.44, 15.76, 10.9, 8.23),
                 se = c(0.28, 0.68, .17, .24, .12))
df <- data.frame(temp = c(15, 20, 25, 27.5),
                 n = c(38, 60, 36, 72),
                 days = c(43.29, 25.54, 17.15, 15.09),
                 se = c(.65, .22, .23, .12))


ptdf <- df %>%
  filter(temp > 15) %>% 
  group_by(temp) %>%
  do(simdata = rnorm(n = .$n, mean = .$days, sd = .$se * sqrt(.$n))) %>%
  unnest()

ptdf %>% group_by(temp) %>% 
  summarise(sd = sd(x = simdata * (temp - 12.2)), 
            dd = mean(simdata * (temp - 12.2)))

ggplot(ptdf, aes(x = simdata)) +
  geom_density() +
  facet_wrap(~ temp)

# Emergence distribution ----


dat <- read.csv("data/GCA_pheno_2013.csv", header = TRUE)
dat$Date <- mdy(dat$Date)
dat <- dat %>% 
  group_by(Site, Stage) %>% 
  mutate(RelProp = MeanCount / max(MeanCount, na.rm = TRUE),
         CDF = cumsum(MeanCount) / sum(MeanCount, na.rm = TRUE),
         Count = ifelse(Stage == "adult",   ### this only works for Fritzi's data
                        round(MeanCount * 15), round(MeanCount * 150)))

# update gdd

library(GSODR)
bs_stations <- nearest_stations(LAT = 44.98, 
                                  LON = -123.27,
                                  distance = 25)
su_stations <- nearest_stations(LAT = 43.387721, 
                                LON = -123.315854,
                                distance = 20)

tbar <- get_GSOD(years = 2013, station = bs_stations) %>% 
  dplyr::select(YDAY, MAX, MIN) %>% 
  mutate(DD10 = TriDD(tmax = MAX, tmin = MIN, LDT = 10, UDT = 30)) %>% 
  mutate(DD12 = TriDD(tmax = MAX, tmin = MIN, LDT = 12.2, UDT = 30)) %>% 
  arrange(YDAY) %>% 
  mutate(AccumDD10 = cumsum(DD10),
         AccumDD12 = cumsum(DD12),
         Site = "MorganLake")
tbar2 <- get_GSOD(years = 2013, station = su_stations) %>% 
  dplyr::select(YDAY, MAX, MIN) %>% 
  mutate(DD10 = TriDD(tmax = MAX, tmin = MIN, LDT = 10, UDT = 30)) %>% 
  mutate(DD12 = TriDD(tmax = MAX, tmin = MIN, LDT = 12.2, UDT = 30)) %>% 
  arrange(YDAY) %>% 
  mutate(AccumDD10 = cumsum(DD10),
         AccumDD12 = cumsum(DD12),
         Site = "Sutherlin")
tbar <- bind_rows(tbar,tbar2)

datcut <- dat %>% 
  filter(Site %in% c("MorganLake", "Sutherlin")) %>% 
  mutate(ydays = yday(Date)) %>% 
  left_join(tbar, by = c("ydays" = "YDAY", "Site")) %>% 
  droplevels()

plt <- ggplot(datcut, aes(x = AccumDD12, y = MeanCount, group = Site, color = Site)) + 
  geom_point() +
  geom_line() +
  scale_x_continuous(limits = c(0, NA)) +
  facet_grid(~Stage) #+
# theme_bw() #+ 
# ggtitle("Cumulative distribution function of observations")
plt
 


library(fGarch)
hist(rsnorm(10000, mean = 130, sd = 25, xi = 3), breaks = 100)


