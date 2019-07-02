pkgs <- c("lubridate", "dplyr", "tidyr",
          "stringr", "purrr", "ggplot2", "viridis")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them

theme_set(theme_bw(base_size = 18)) 


dat <- read.csv("data/galerucella_chamber_emergence_2018.csv", header = TRUE)

# reformat data from wide to long
df <- dat %>% 
  gather(date, count, X8.6.2018:X9.6.2018) %>% 
  mutate(date = mdy(gsub(pattern = "X", replacement = "", x = date)),
         ovip_date = dmy(paste0(ovip_date, "-2018"))) %>% 
  rowwise() %>% 
  mutate(degdays = (temperature - 10) * as.numeric(date - ovip_date)) %>% 
  group_by(chamber, temperature, photoperiod, population, degdays) %>% 
  summarise(emerged = sum(count, na.rm = TRUE)) %>% 
  filter(degdays <= 550 & degdays >= 300) %>% 
  group_by(chamber, population) %>% 
  mutate(tot_emerged = sum(emerged),
         prop_emerged = emerged / tot_emerged)


levels(df$population) <- c("Bellingham, WA", "Baskett-Slough, OR", "Ft. Drum, NY", "Montesano, WA",
                           "McArthur, CA", "Palermo, CA", "Sutherlin, OR", "Coeburn, VA", "West Point, NY",
                           "Yakima, WA")

# what's the average degree-day for emergence excluding cool chamber 2?
avgemerg <- df %>% 
  filter(chamber %in% c(1, 4, 5, 6)) %>% 
  group_by(population) %>% 
  summarise(meandd = weighted.mean(degdays, emerged)) %>% 
  arrange(meandd)
avgemerg

# plot of each chamber's counts by population
plt <- ggplot(df, aes(x = degdays, y = emerged, group = population, color = population)) + 
  geom_point() +
  # geom_smooth(method = "gam") +
  facet_wrap(~chamber)
plt

# better plot, emergence by chamber for each population in a separate facet
plt <- ggplot(df, aes(x = degdays, y = prop_emerged, group = interaction(population, chamber), color = as.factor(chamber))) + 
  geom_point() +
  geom_smooth(span = .15, se = FALSE) +
  facet_wrap(~population) + 
  ggtitle("Adult emergence in degree-days by growth chamber\n#2 is slightly cool, #3 is 26C treatment") +
  xlab("Cumulative degree-days base 10C") +
  ylab("Proportion of total emerged by chamber x population") +
  theme(legend.position = c(0.9, 0.175))
plt


df2 <- df %>% 
  filter(population %in% unique(df$population)[c(1, 2, 4, 7, 10)], chamber >= 4)
plt <- ggplot(df2, aes(x = degdays, y = prop_emerged, group = interaction(population, chamber), color = as.factor(chamber))) + 
  geom_point() +
  scale_color_viridis(name = "Chamber", discrete = TRUE, begin = .25, end = .75) +
  geom_smooth(span = .15, se = FALSE) +
  facet_wrap(~population) + 
  ggtitle("Adult emergence in degree-days by growth chamber") +
  xlab("Cumulative degree-days base 10C") +
  ylab("Proportion of total emerged by chamber") +
  theme(legend.position = c(0.85, 0.2), plot.title = element_text(hjust = 0.5))
plt



df3 <- df %>% 
  filter(chamber >= 4, degdays < 500) %>% 
  group_by(population, degdays) %>% 
  summarise(emerged = sum(emerged)) %>% 
  group_by(population) %>% 
  arrange(degdays) %>% 
  mutate(tot_emerged = sum(emerged),
         cum_emerged = cumsum(emerged),
         cdf = cum_emerged / tot_emerged)
  
  
plt <- ggplot(df3, aes(x = degdays, y = cdf, group = population, color = population)) +
  geom_line() +
  xlab("Accumulated degree-days (10C)") +
  ylab("Cumulative proportion emerged") +
  ggtitle("Development time to adult emergence in growth chamber") +
  theme(legend.position = c(0.8, 0.25), plot.title = element_text(hjust = 0.5))
plt

  
# 
# 
# # any difference in oviposition date?
# df <- dat %>% 
#   gather(date, count, X8.6.2018:X9.6.2018) %>% 
#   mutate(date = mdy(gsub(pattern = "X", replacement = "", x = date)),
#          ovip_date = dmy(paste0(ovip_date, "-2018"))) %>% 
#   rowwise() %>% 
#   mutate(degdays = (temperature - 10) * as.numeric(date - ovip_date)) %>% 
#   filter(chamber %in% c(1, 4, 5, 6)) %>% 
#   # group_by(temperature, photoperiod, population, degdays, ovip_date) %>% 
#   group_by(temperature, photoperiod, population, degdays) %>% 
#   summarise(emerged = sum(count, na.rm = TRUE)) %>% 
#   filter(degdays <= 550 & degdays >= 300) %>% 
#   group_by(population) %>% 
#   # group_by(population, ovip_date) %>% 
#   mutate(tot_emerged = sum(emerged),
#          prop_emerged = emerged / tot_emerged) %>% 
#   group_by(population) %>% 
#   mutate(ovip_day = as.numeric(ovip_date - max(ovip_date)))
# 
# levels(df$population) <- c("Bellingham, WA", "Baskett-Slough, OR", "Ft. Drum, NY", "Montesano, WA",
#                            "McArthur, CA", "Palermo, CA", "Sutherlin, OR", "Coeburn, VA", "West Point, NY",
#                            "Yakima, WA")
# 
# df$ovip_day <- factor(as.character(df$ovip_day), c("-2", "-1", "0"))
# levels(df$ovip_day) <- c("1st", "2nd", "3rd")
# 
# 
# plt <- ggplot(df, aes(x = degdays, y = prop_emerged, group = interaction(population, ovip_day), color = ovip_day)) + 
#   geom_point() +
#   geom_smooth(span = .15, se = FALSE) +
#   facet_wrap(~population) + 
#   ggtitle("Adult emergence in degree-days for chambers 1, 4, 5, & 6\ngrouped by day of oviposition") +
#   xlab("Cumulative degree-days base 10C") +
#   ylab("Proportion of total emerged in chamber") +
#   theme(legend.position = c(0.9, 0.15))
# plt
# 
# df %>% 
#   group_by(population, ovip_day) %>% 
#   summarise(meandd = weighted.mean(degdays, emerged)) %>% 
#   arrange(meandd)
# 
# 
# plt <- ggplot(df, aes(x = degdays, y = prop_emerged)) + 
#   geom_point() +
#   geom_smooth(span = .15, se = FALSE) +
#   facet_wrap(~population) + 
#   ggtitle("Adult emergence in degree-days by population") +
#   xlab("Cumulative degree-days base 10C") +
#   ylab("Proportion of total emerged") +
#   theme(legend.position = c(0.9, 0.15))
# plt
