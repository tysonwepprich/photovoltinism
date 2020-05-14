pkgs <- c("lubridate", "dplyr", "tidyr",
          "stringr", "purrr", "ggplot2", "viridis")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them

'%!in%' <- function(x,y)!('%in%'(x,y))

theme_set(theme_bw(base_size = 18) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())) 

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

cdfemerg <- df %>% 
  filter(chamber %in% c(1, 4, 5, 6)) %>% 
  group_by(population, photoperiod) %>% 
  arrange(population, degdays) %>% 
  mutate(cumemerg = cumsum(emerged),
         cumprop = cumemerg / tot_emerged)


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



# Pop x photo interaction ----

df <- dat %>% 
  gather(date, count, X8.6.2018:X9.6.2018) %>% 
  mutate(date = mdy(gsub(pattern = "X", replacement = "", x = date)),
         ovip_date = dmy(paste0(ovip_date, "-2018"))) %>% 
  rowwise() %>% 
  mutate(days = as.numeric(date - ovip_date),
         degdays = (temperature - 10) * as.numeric(date - ovip_date),
         id = paste(chamber, population, petri_num, sep = "_")) %>% 
  filter(chamber %!in% c(2, 3)) %>% 
  filter(!is.na(count)) %>% 
  filter(population %in% c("S", "Y", "MC", "BS", "BL", "P")) %>% 
  droplevels.data.frame() %>% 
  group_by(chamber, population, days) %>% 
  mutate(emerged = sum(count, na.rm = TRUE)) %>%
  group_by(chamber, population) %>% 
  mutate(tot_emerged = sum(emerged, na.rm = TRUE),
         prop_emerged = emerged / tot_emerged)

levels(df$population) <- c("Bellingham, WA", "Rickreall, OR", 
                            "McArthur, CA", "Palermo, CA", "Sutherlin, OR",
                            "Yakima TC, WA")

df$population <- factor(df$population, levels = c("Bellingham, WA", "Yakima TC, WA", "Rickreall, OR",
                                                          "Sutherlin, OR", "McArthur, CA", "Palermo, CA"))

df2 <- df %>% filter(days > 27 & days < 45) %>% 
  group_by(population, days) %>% 
  summarise(emerged = sum(count, na.rm = TRUE)) %>% 
  group_by(population) %>% 
  mutate(tot_emerged = sum(emerged, na.rm = TRUE),
         prop_emerged = emerged / tot_emerged)

plt <- ggplot(df2, aes(x = days, y = prop_emerged, group = population, shape = population)) + 
  # geom_smooth(span = .15, se = FALSE) +
  geom_point(size = 2) +
  ggalt::geom_xspline(alpha = .5) +
  scale_shape_discrete(name = "Population") +
  # scale_color_brewer(name = "Source population", palette = "Dark2") +
  scale_y_continuous(limits = c(-0.005, .42), expand = c(0, 0)) +
  xlab("Days from egg to adult") +
  ylab("Daily proportion eclosed") +
  theme(legend.position = c(0.75, 0.75), plot.title = element_text(hjust = 0.5))
plt


df3 <- df %>% filter(days > 27 & days < 45) %>% 
  group_by(chamber, days) %>% 
  summarise(emerged = sum(count, na.rm = TRUE)) %>% 
  group_by(chamber) %>% 
  mutate(tot_emerged = sum(emerged, na.rm = TRUE),
         prop_emerged = emerged / tot_emerged) %>% 
  ungroup() %>% 
  mutate(chamber = as.factor(chamber))
df3$chamber <- plyr::revalue(df3$chamber, 
                             c("1"="14.75/9.25",
                               "4"="15.75/8.25",
                               "5"="16.25/7.75",
                               "6"="16.75/7.25"))


plt2 <- ggplot(df3, aes(x = days, y = prop_emerged, group = chamber, shape = chamber)) + 
  # geom_smooth(span = .15, se = FALSE) +
  geom_point(size = 2) +
  ggalt::geom_xspline(alpha = 0.5) +
  scale_shape_discrete(name = "Photoperiod (L/D)") +
  # scale_color_brewer(name = "Photoperiod (L/D)", palette = "Set1") +
  scale_y_continuous(limits = c(-0.005, .42), expand = c(0, 0)) +
  xlab("Days from egg to adult") +
  ylab("Daily proportion eclosed") +
  theme(legend.position = c(0.75, 0.8), plot.title = element_text(hjust = 0.5))
plt2

library(ggpubr)
a <- ggarrange(plt, plt2, 
               labels = c("A", "B"),
               ncol = 2, nrow = 1, align = "h", common.legend = FALSE)
a

ggsave(filename = "fig4.tif", plot = a, device = "tiff", width = 12, height = 6, units = "in", dpi = 1200)


df.expanded <- df[rep(row.names(df), df$count), ]

df.expanded$zday <- as.vector(scale(df.expanded$photoperiod))
df.expanded$chamber <- as.factor(df.expanded$chamber)

library(lme4)
moda <- glm(days ~ zday + population, family = poisson, data = df.expanded)
modb <- glm(days ~ zday + population, family = quasipoisson, data = df.expanded)


mod2 <- lmer(log(days) ~ zday + population + (1 | id), data = df.expanded)

# singular fit
mod <- glmer(days ~ zday + (1 + zday|population) + (1|id), data = df.expanded,
             family = Gamma(link = "log"))

mod2 <- lmer(log(days) ~ chamber + (1 | population/id), data = df.expanded)
mod <- glmer(days ~ chamber + population + (1|id), data = df.expanded,
             family = Gamma(link = "log"), glmerControl(optimizer = "bobyqa"))
             # family = poisson(link = "log"), glmerControl(optimizer = "bobyqa"))
mod <- glmer(days ~ zday + population -1+ (1|id), data = df.expanded,
             family = Gamma(link = "log"), glmerControl(optimizer = "bobyqa"))
# family = poisson(link = "log"), glmerControl(optimizer = "bobyqa"))


drop1(mod, .~., test = "Chisq")
car::Anova(mod, type = 2)
sjPlot::tab_model(mod)

library(ggridges)
dplt <- ggplot(df, aes(x = days, y = prop_emerged)) +
  geom_point() +
  facet_wrap(~population, ncol = 2)
dplt

# Test how pairwise populations are different from each other
library(emmeans)
cells <- emmeans(mod, pairwise ~ population, adjust = "tukey")
pp <- pwpp(cells$emmeans, type = "response", method = "pairwise")
pp + ylab("Population")
emmeans::CLD(cells)

# Interaction plots of results
library(ggeffects)
gg <- ggpredict(mod, c("population", "zday")) %>% plot()
plt <- gg +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle(NULL) +
  xlab("Population") +
  ylab("Predicted days to develop at 23C") +
  scale_color_grey(name = "Day length\nin hours", start = 0.1, end = 0.7)

plt
ggsave(filename = "GCA_GC_2018.png", plot = plt, device = "png", width = 9, height = 5, units = "in")


library(fitdistrplus)

descdist(df.expanded$days, discrete = TRUE)
descdist(residuals(modb), discrete = FALSE)

fit.beta <- fitdist(df.expanded$days, "lognormal", method = "mme")


# just do boxplots for distributions
plt <- ggplot(df.expanded, aes(x = population, y = days)) +
  geom_histogram()
plt

