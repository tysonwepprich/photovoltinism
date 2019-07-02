pkgs <- c("lubridate", "dplyr", "tidyr", "lme4",
          "stringr", "purrr", "ggplot2", "viridis")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them

theme_set(theme_bw(base_size = 16)) 


dat <- read.csv("data/galerucella_chamber_diapause_2018.csv", header = TRUE)
dat2 <- read.csv("data/galerucella_chamber_diapause_2019.csv", header = TRUE)

dat <- bind_rows(dat, dat2)

dat <- dat %>% 
  rowwise() %>% 
  mutate(total = repro_female + feeder + diapause,
         perc_repro = 1 - diapause/total,
         rowid = paste(population, chamber, emerg_group, sep = "_"))  %>%
  filter(starting_lifestage == "egg") %>% 
  droplevels.data.frame()
  # arrange(population, chamber, starting_lifestage, emerg_group) %>% 
  # filter(population == "S") %>% 
  # data.frame()
dat$population <- as.factor(dat$population)

levels(dat$population) <- c("Bellingham, WA", "Baskett-Slough, OR", "Ft. Drum, NY", "Montesano, WA",
                           "McArthur, CA", "Palermo, CA", "Sutherlin, OR", "Coeburn, VA", "West Point, NY",
                           "Yakima, WA")

# scaling photoperiod for model fitting
dat$zday <- scale(dat$photoperiod)

# glm with random effects
# response is the percent reproductive, with the total count for the dish as the model weights (could do this other ways)
# trying out different model random effects
# works best with random intercepts and slopes by population (mod1)
mod <- glmer(perc_repro ~ population + zday + temperature + (1|rowid), weights = total, data = dat, family = binomial)
mod1 <- glmer(perc_repro ~ zday + temperature + (1 + zday|population), weights = total, data = dat, family = binomial)
mod2 <- glmer(perc_repro ~ zday + temperature + (1|population), weights = total, data = dat, family = binomial)
summary(mod2)
AIC(mod, mod1, mod2)

# just modeling 21C chambers here
dat_cut <- dat %>% filter(temperature == 23)
mod1 <- glmer(perc_repro ~ zday + (1 + zday|population), weights = total, data = dat_cut, family = binomial)


# model predictions for plotting
newdat <- expand.grid(list(zday = seq(-1.71, 1.71, length.out = 100),
                           population = unique(dat_cut$population),
                           temperature = 23))
newdat$pred <- predict(mod1, newdat, type = "response")
# back transform zday to photoperiod
newdat$photoperiod <- newdat$zday * attr(dat$zday, 'scaled:scale') + attr(dat$zday, 'scaled:center')
newdat <- newdat %>% filter(temperature == 23)


# plot shows raw data for all chambers and model prediction for 21C chambers by population
plt <- ggplot(dat_cut, aes(x = photoperiod, y = perc_repro)) +
  geom_point(aes(color = emerg_group, shape = as.factor(temperature)), size = 4, alpha = .5) +
  scale_color_discrete(name = "Petri dish") +
  scale_shape_discrete(name = "Temperature") +
  geom_line(data = newdat, aes(x = photoperiod, y = pred)) +
  xlab("Daylength in chamber") +
  ylab("Percent reproductive") +
  facet_wrap(~population, ncol = 3)
plt
ggsave(filename = "GCA_GC2019.png", plot = plt, device = "png", width = 12, height = 9, units = "in")

# group early/late petris together
dat <- read.csv("C:/Users/wepprict/Desktop/galerucella_chamber_diapause_2018.csv", header = TRUE)

dat <- dat %>% 
  filter(starting_lifestage == "egg") %>% 
  group_by(chamber, temperature, photoperiod, population) %>% 
  summarise(total = sum(repro_female, feeder, diapause),
         perc_repro = 1 - sum(diapause)/total,
         rowid = paste(population[1], chamber[1], sep = "_"))  %>%
  droplevels.data.frame()
# arrange(population, chamber, starting_lifestage, emerg_group) %>% 
# filter(population == "S") %>% 
# data.frame()

levels(dat$population) <- c("Bellingham, WA", "Baskett-Slough, OR", "Ft. Drum, NY", "Montesano, WA",
                            "McArthur, CA", "Palermo, CA", "Sutherlin, OR", "Coeburn, VA", "West Point, NY",
                            "Yakima, WA")

plt <- ggplot(dat, aes(x = photoperiod, y = perc_repro)) +
  geom_point(aes(shape = as.factor(temperature)), size = 3, alpha = .5) +
  scale_color_discrete(name = "Petri dish") +
  scale_shape_discrete(name = "Temperature") +
  geom_line(data = newdat, aes(x = photoperiod, y = pred)) +
  xlab("Daylength in chamber") +
  ylab("Percent reproductive") +
  facet_wrap(~population, ncol = 3)
plt


# what about combining with Fritzi previous chamber experiments?
dat <- read.csv("data/galerucella_chamber_diapause_2018.csv", header = TRUE) %>% 
  mutate(Year = 2018)
dat2 <- read.csv("data/galerucella_chamber_diapause_2019.csv", header = TRUE) %>% 
  mutate(Year = 2019)

dat <- bind_rows(dat, dat2)

dat <- dat %>% 
  rowwise() %>% 
  mutate(total = repro_female + feeder + diapause,
         perc_repro = 1 - diapause/total,
         rowid = paste(population, chamber, emerg_group, sep = "_"))  %>%
  filter(starting_lifestage == "egg") %>% 
  droplevels.data.frame()
# arrange(population, chamber, starting_lifestage, emerg_group) %>% 
# filter(population == "S") %>% 
# data.frame()
dat$population <- as.factor(dat$population)

levels(dat$population) <- c("Bellingham, WA", "Baskett-Slough, OR", "Ft. Drum, NY", "Montesano, WA",
                            "McArthur, CA", "Palermo, CA", "Sutherlin, OR", "Coeburn, VA", "West Point, NY",
                            "Yakima, WA")


dat2 <- read.csv('data/combined 13_15.csv')
dat2$population <- dat2$Site
levels(dat2$population) <- c("Bellingham, WA", "Yakima, WA", "Palermo, CA", "Sutherlin, OR")

alldat <- dat2 %>% 
  group_by(population, Treatment, Year) %>% 
  summarise(repro = sum(diapause),
            total = n(),
            perc_repro = repro/total,
            photoperiod = Treatment[1],
            temperature = 23) %>% 
  bind_rows(dat)


# scaling photoperiod for model fitting
alldat$zday <- scale(alldat$photoperiod)
alldat$ztemp <- scale(alldat$temperature)

mod1 <- glmer(perc_repro ~ zday + ztemp + (1 + zday|population), weights = total, data = alldat, family = binomial)
mod2 <- glmer(perc_repro ~ zday * ztemp + (1 + zday|population), weights = total, data = alldat, family = binomial)
mod3 <- glmer(perc_repro ~ zday * ztemp + (1|population), weights = total, data = alldat, family = binomial)
mod4 <- glmer(perc_repro ~ zday * ztemp + (1 + zday|population) + (1|Year), weights = total, data = alldat, family = binomial)
mod5 <- glmer(perc_repro ~ zday * ztemp + (1 + zday + ztemp|population) + (1|Year), weights = total, data = alldat, family = binomial)
mod6 <- glmer(perc_repro ~ zday * ztemp + (1 + zday + ztemp|population) + (1 + zday|Year), weights = total, data = alldat, family = binomial)

summary(mod4)
AIC(mod1, mod2, mod3, mod4, mod5, mod6)

# model predictions for plotting
newdat <- expand.grid(list(zday = seq(-1.75, 1.75, length.out = 100),
                           ztemp = 0.431,
                           population = unique(alldat$population),
                           Year = NA))
newdat$pred <- predict(mod4, newdat, re.form = ~(1 + zday|population), type = "response")
# back transform zday to photoperiod
newdat$photoperiod <- newdat$zday * attr(alldat$zday, 'scaled:scale') + attr(alldat$zday, 'scaled:center')
newdat$temperature <- newdat$ztemp * attr(alldat$ztemp, 'scaled:scale') + attr(alldat$ztemp, 'scaled:center')

# model predictions for plotting each Year/temperature separate lines
linedat <- expand.grid(list(zday = seq(-1.75, 1.75, length.out = 100),
                           ztemp = unique(alldat$ztemp),
                           population = unique(alldat$population),
                           Year = c(2013, 2014, 2018, 2019)))
linedat$pred <- predict(mod4, linedat, re.form = ~(1 + zday|population) + (1|Year), type = "response")
# back transform zday to photoperiod
linedat$photoperiod <- linedat$zday * attr(alldat$zday, 'scaled:scale') + attr(alldat$zday, 'scaled:center')
linedat$temperature <- linedat$ztemp * attr(alldat$ztemp, 'scaled:scale') + attr(alldat$ztemp, 'scaled:center')

linedat$temperature <- round(linedat$temperature)
linedat <- linedat %>% filter(temperature == 23)


# plot shows raw data for all chambers and model prediction for 21C chambers by population
plt <- ggplot(alldat, aes(x = photoperiod, y = perc_repro)) +
  geom_point(aes(color = as.factor(Year), shape = as.factor(temperature)), size = 4, alpha = .5) +
  # scale_color_viridis(discrete = TRUE) +
  geom_line(data = newdat, aes(x = photoperiod, y = pred), size = 1.5) +
  geom_hline(yintercept = .5, linetype = "dashed") +
  # geom_line(data = linedat, aes(x = photoperiod, y = pred,
  #                               group = interaction(Year, temperature),
  #                               color = as.factor(Year)))+
  scale_color_discrete(name = "Year") +
  scale_shape_discrete(name = "Temperature") +
  # geom_line(data = newdat, aes(x = photoperiod, y = pred)) +
  xlab("Daylength in chamber") +
  ylab("Percent reproductive") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~population, ncol = 3)
plt


findInt <- function(model, value) {
  function(x) {
    predict(model, data.frame(Treat=x), type="response") - value
  }
}

uniroot(findInt(modS, .5), range(datS$Treat))$root
uniroot(findInt(modS, .95), range(datS$Treat))$root
uniroot(findInt(modS, .05), range(datS$Treat))$root


cdl <- newdat %>% 
  group_by(population) %>% 
  summarise(cdl = as.numeric(photoperiod[which.min(abs(pred - .5))]),
            cdl_sd = (cdl - as.numeric(photoperiod[which.min(abs(pred - .05))]))/2) %>% 
  mutate(optcdl = c(16, NA, 15, 14.5, 15.25, NA, 14.5, 14.75, NA, NA, NA)) %>% 
  filter(complete.cases(.)) %>% 
  tidyr::gather(var, cp, cdl,optcdl) %>% 
  mutate(type = c(rep("growth chamber", 6), rep("optimal", 6)))

cdl <- bind_rows(cdl, data.frame(population = "Germany", cdl_sd = 0, var = "cdlgerm",
                                 cp = 17, type = "source optimal"))
cdl$type <- factor(cdl$type, levels = c("source optimal", "growth chamber", "optimal"))

library(ggrepel)
ggplot(cdl, aes(x = type, y = cp, color = population, group = population)) +
  geom_point() + 
  geom_line() +
  geom_text_repel(aes(label = population)) + 
  theme(legend.position = "none") +
  xlab(NULL) + 
  ylab("Critical photoperiod")



newdat <- newdat %>% 
  filter(population %in% unique(population)[c(1, 2, 3, 4, 6, 7)])
newdat <- droplevels.data.frame(newdat)
newdat$population = factor(newdat$population, levels(newdat$population)[c(6,1,5,2,3,4)])
cdl <- newdat %>% 
  group_by(population) %>% 
  summarise(photoperiod = as.numeric(photoperiod[which.min(abs(pred - .5))]))


plt <- ggplot(newdat, aes(x = photoperiod, y = pred, group = population)) +
  geom_line(data = newdat, aes(x = photoperiod, y = pred, color = population), 
            size = 1.5, alpha = .5) +
  # scale_color_discrete(guide = FALSE) +
  geom_point(data = cdl, aes(x = photoperiod, y = 0.5), inherit.aes = FALSE) +
  geom_label_repel(data = cdl, aes(x = photoperiod, y = 0.5,label = round(photoperiod, 2)), inherit.aes = FALSE) +
  xlab("Hours of light in chamber (23C)") +
  ylab("Percent reproductive (modeled)") +
  ggtitle("Critical photoperiod by NW population") +
  theme(legend.position = c(.2, .8))
plt



# Try new comparison of CP curve model predictions
# PNW populations with site/years analyzed separately
# Analyzed by year
# Analyzed all together


dat <- read.csv("data/galerucella_chamber_diapause_2018.csv", header = TRUE) %>% 
  mutate(Year = 2018)
dat2 <- read.csv("data/galerucella_chamber_diapause_2019.csv", header = TRUE) %>% 
  mutate(Year = 2019)

dat <- bind_rows(dat, dat2)

dat <- dat %>% 
  rowwise() %>% 
  mutate(total = repro_female + feeder + diapause,
         perc_repro = 1 - diapause/total,
         rowid = paste(population, chamber, emerg_group, sep = "_"))  %>%
  filter(starting_lifestage == "egg") %>% 
  droplevels.data.frame()
# arrange(population, chamber, starting_lifestage, emerg_group) %>% 
# filter(population == "S") %>% 
# data.frame()
dat$population <- as.factor(dat$population)

levels(dat$population) <- c("Bellingham, WA", "Baskett-Slough, OR", "Ft. Drum, NY", "Montesano, WA",
                            "McArthur, CA", "Palermo, CA", "Sutherlin, OR", "Coeburn, VA", "West Point, NY",
                            "Yakima, WA")


dat2 <- read.csv('data/combined 13_15.csv')
dat2$population <- dat2$Site
levels(dat2$population) <- c("Bellingham, WA", "Yakima, WA", "Palermo, CA", "Sutherlin, OR")

alldat <- dat2 %>% 
  group_by(population, Treatment, Year) %>% 
  summarise(repro = sum(diapause),
            total = n(),
            perc_repro = repro/total,
            photoperiod = Treatment[1],
            temperature = 23) %>% 
  bind_rows(dat)

dat23 <- alldat %>% filter(temperature == 23) %>% droplevels.data.frame()
# plot each site x year separately
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

plt <- ggplot(dat23, aes(x = photoperiod, y = perc_repro, color = log(total), group = interaction(population, Year))) +
  geom_point() +
  binomial_smooth(aes(weight = total), formula = y ~ x, se = TRUE) +
  facet_grid(population ~ Year) +
  geom_hline(yintercept = .5, linetype = "dashed") +
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("All 23C chamber experiments, each site x year analyzed separately")
plt

plt <- ggplot(dat23, aes(x = photoperiod, y = perc_repro, group = population)) +
  binomial_smooth(aes(weight = total), formula = y ~ x, se = TRUE) +
  geom_point(aes(color = as.factor(Year))) +
  facet_wrap(~population) +
  geom_hline(yintercept = .5, linetype = "dashed") +
  theme(text = element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("All 23C chamber experiments, each site analyzed separately")
plt


# scaling photoperiod for model fitting
alldat$zday <- scale(alldat$photoperiod)
alldat$ztemp <- scale(alldat$temperature)

mod1 <- glmer(perc_repro ~ zday + ztemp + (1 + zday|population), weights = total, data = alldat, family = binomial)
mod2 <- glmer(perc_repro ~ zday * ztemp + (1 + zday|population), weights = total, data = alldat, family = binomial)
mod3 <- glmer(perc_repro ~ zday * ztemp + (1|population), weights = total, data = alldat, family = binomial)
mod4 <- glmer(perc_repro ~ zday * ztemp + (1 + zday|population) + (1|Year), weights = total, data = alldat, family = binomial)
mod5 <- glmer(perc_repro ~ zday * ztemp + (1 + zday + ztemp|population) + (1|Year), weights = total, data = alldat, family = binomial)
mod6 <- glmer(perc_repro ~ zday * ztemp + (1 + zday + ztemp|population) + (1 + zday|Year), weights = total, data = alldat, family = binomial)

summary(mod4)
AIC(mod1, mod2, mod3, mod4, mod5, mod6)

# model predictions for plotting
newdat <- expand.grid(list(zday = seq(-1.75, 1.75, length.out = 100),
                           ztemp = 0.431,
                           population = unique(alldat$population),
                           Year = NA))
newdat$pred <- predict(mod4, newdat, re.form = ~(1 + zday|population), type = "response")
# back transform zday to photoperiod
newdat$photoperiod <- newdat$zday * attr(alldat$zday, 'scaled:scale') + attr(alldat$zday, 'scaled:center')
newdat$temperature <- newdat$ztemp * attr(alldat$ztemp, 'scaled:scale') + attr(alldat$ztemp, 'scaled:center')

# model predictions for plotting each Year/temperature separate lines
linedat <- expand.grid(list(zday = seq(-1.75, 1.75, length.out = 100),
                            ztemp = unique(alldat$ztemp),
                            population = unique(alldat$population),
                            Year = c(2013, 2014, 2018, 2019)))
linedat$pred <- predict(mod4, linedat, re.form = ~(1 + zday|population) + (1|Year), type = "response")
# back transform zday to photoperiod
linedat$photoperiod <- linedat$zday * attr(alldat$zday, 'scaled:scale') + attr(alldat$zday, 'scaled:center')
linedat$temperature <- linedat$ztemp * attr(alldat$ztemp, 'scaled:scale') + attr(alldat$ztemp, 'scaled:center')

linedat$temperature <- round(linedat$temperature)
linedat <- linedat %>% filter(temperature == 23)


# plot shows raw data for all chambers and model prediction for 21C chambers by population
plt <- ggplot(alldat, aes(x = photoperiod, y = perc_repro)) +
  geom_point(aes(color = as.factor(Year), shape = as.factor(temperature)), size = 4, alpha = .5) +
  # scale_color_viridis(discrete = TRUE) +
  geom_line(data = newdat, aes(x = photoperiod, y = pred), size = 1.5) +
  geom_line(data = linedat, aes(x = photoperiod, y = pred,
                                group = interaction(Year, temperature),
                                color = as.factor(Year)))+
  scale_color_discrete(name = "Year") +
  scale_shape_discrete(name = "Temperature") +
  # geom_line(data = newdat, aes(x = photoperiod, y = pred)) +
  xlab("Daylength in chamber") +
  ylab("Percent reproductive") +
  facet_wrap(~population, ncol = 3)
plt




