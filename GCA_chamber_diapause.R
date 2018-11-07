pkgs <- c("lubridate", "dplyr", "tidyr", "lme4",
          "stringr", "purrr", "ggplot2", "viridis")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them

theme_set(theme_bw(base_size = 20)) 


dat <- read.csv("data/galerucella_chamber_diapause_2018.csv", header = TRUE)

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
summary(mod)
AIC(mod)

# just modeling 21C chambers here
dat21 <- dat %>% filter(temperature == 21)
mod1 <- glmer(perc_repro ~ zday + (1 + zday|population), weights = total, data = dat21, family = binomial)


# model predictions for plotting
newdat <- expand.grid(list(zday = seq(-1.6, 1.6, length.out = 100),
                           population = unique(dat$population),
                           temperature = c(21, 26)))
newdat$pred <- predict(mod1, newdat, type = "response")
# back transform zday to photoperiod
newdat$photoperiod <- newdat$zday * attr(dat$zday, 'scaled:scale') + attr(dat$zday, 'scaled:center')
newdat <- newdat %>% filter(temperature == 21)


# plot shows raw data for all chambers and model prediction for 21C chambers by population
plt <- ggplot(dat, aes(x = photoperiod, y = perc_repro)) +
  geom_point(aes(color = emerg_group, shape = as.factor(temperature)), size = 4, alpha = .5) +
  scale_color_discrete(name = "Petri dish") +
  scale_shape_discrete(name = "Temperature") +
  geom_line(data = newdat, aes(x = photoperiod, y = pred)) +
  xlab("Daylength in chamber") +
  ylab("Percent reproductive") +
  facet_wrap(~population, ncol = 3)
plt





# what about combining with Fritzi previous chamber experiments?
dat <- read.csv("C:/Users/wepprict/Desktop/galerucella_chamber_diapause_2018.csv", header = TRUE)

dat <- dat %>% 
  rowwise() %>% 
  mutate(total = repro_female + feeder + diapause,
         perc_repro = 1 - diapause/total,
         rowid = paste(population, chamber, emerg_group, sep = "_"),
         Year = 2018)  %>%
  filter(starting_lifestage == "egg") %>% 
  droplevels.data.frame()
# arrange(population, chamber, starting_lifestage, emerg_group) %>% 
# filter(population == "S") %>% 
# data.frame()

levels(dat$population) <- c("Bellingham, WA", "Baskett-Slough, OR", "Ft. Drum, NY", "Montesano, WA",
                            "McArthur, CA", "Palermo, CA", "Sutherlin, OR", "Coeburn, VA", "West Point, NY",
                            "Yakima, WA")


dat2 <- read.csv('data/combined 13_15.csv')
dat2$population <- dat2$Site
levels(dat2$population) <- c("Bellingham, WA", "Ephrata, WA", "Palermo, CA", "Sutherlin, OR")

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

summary(mod)
AIC(mod1, mod2, mod3, mod4)

# model predictions for plotting
newdat <- expand.grid(list(zday = seq(-1.6, 1.6, length.out = 100),
                           ztemp = 0.554,
                           population = unique(alldat$population),
                           Year = NA))
newdat$pred <- predict(mod4, newdat, re.form = ~(1 + zday|population), type = "response")
# back transform zday to photoperiod
newdat$photoperiod <- newdat$zday * attr(alldat$zday, 'scaled:scale') + attr(alldat$zday, 'scaled:center')
newdat$temperature <- newdat$ztemp * attr(alldat$ztemp, 'scaled:scale') + attr(alldat$ztemp, 'scaled:center')

# model predictions for plotting each Year/temperature separate lines
linedat <- expand.grid(list(zday = seq(-1.6, 1.6, length.out = 100),
                           ztemp = c(0.554, -0.578, 2.253),
                           population = unique(alldat$population),
                           Year = c(2013, 2014, 2018)))
linedat$pred <- predict(mod4, linedat, re.form = ~(1 + zday|population) + (1|Year), type = "response")
# back transform zday to photoperiod
linedat$photoperiod <- linedat$zday * attr(alldat$zday, 'scaled:scale') + attr(alldat$zday, 'scaled:center')
linedat$temperature <- linedat$ztemp * attr(alldat$ztemp, 'scaled:scale') + attr(alldat$ztemp, 'scaled:center')

linedat$temperature <- round(linedat$temperature)
linedat <- linedat %>% filter(temperature == 23, Year == 2018)


alldat <- alldat %>% 
  filter(population %in% unique(alldat$population)[c(1, 4, 5, 8, 9)])


# plot shows raw data for all chambers and model prediction for 21C chambers by population
plt <- ggplot(alldat, aes(x = photoperiod, y = perc_repro)) +
  geom_point(aes(color = as.factor(Year), shape = as.factor(temperature)), size = 4, alpha = .5) +
  # scale_color_viridis(discrete = TRUE) +
  # geom_line(data = newdat, aes(x = photoperiod, y = pred), size = 1.5) +
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



alldat <- alldat %>% 
  filter(population %in% unique(alldat$population)[c(1, 4, 5, 8, 9)])



# model predictions for plotting
newdat <- expand.grid(list(zday = seq(-2, 2, length.out = 100),
                           ztemp = 0.554,
                           population = unique(alldat$population)[c(1, 4, 5, 8, 9)],
                           Year = NA))
# newdat$pred <- predict(mod4, newdat, , type = "response")
# back transform zday to photoperiod
newdat$photoperiod <- newdat$zday * attr(alldat$zday, 'scaled:scale') + attr(alldat$zday, 'scaled:center')
newdat$temperature <- newdat$ztemp * attr(alldat$ztemp, 'scaled:scale') + attr(alldat$ztemp, 'scaled:center')

preds <- predict(mod4, newdat, type = "response", re.form = ~(1 + zday|population))
newdat <- cbind(newdat, preds)



cdl <- newdat %>% 
  group_by(population) %>% 
  summarise(cdl = as.numeric(photoperiod[which.min(abs(preds - .5))]))

cdl$preds <- 0.5
cdl$photoperiod <- cdl$cdl

newdat$population <- factor(newdat$population, levels = levels(newdat$population)[c(2, 3, 5, 4, 1)])

library(ggrepel)
# model predictions for different latitudes and 2014 vs 2017
plt <- ggplot(newdat, aes(x = photoperiod, y = preds, group = population)) +
  geom_line(aes(color = population), 
            size = 1.5, alpha = .5) +
  geom_label(data = cdl, aes(label = round(photoperiod, 1))) +
  xlab("Hours of light in chamber (23C)") +
  ylab("Percent reproductive (modeled)") +
  ggtitle("Critical daylengths of Galerucella in PNW") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0.15, 0.85))
plt






# where are the sites?
sites <- dat %>% 
  dplyr::select(site, lat, lon) %>% 
  distinct() %>% 
  data.frame()
sites

library(mapdata)
library(raster)
library(ggrepel)
region_param <- "WEST"

# Derived parameters -----
REGION <- assign_extent(region_param = region_param)

states <- map_data("state", xlim = c(REGION@xmin, REGION@xmax),
                   ylim = c(REGION@ymin, REGION@ymax), lforce = "e")
names(states)[1:2] <- c("x", "y")
mp <- ggplot(data = sites, aes(x = lon, y = lat)) +
  geom_point() +
  geom_label_repel(aes(label = site)) +
  geom_polygon(data = states, aes(x = x, y = y, group = group), fill = NA, color = "black", size = .1) +
  coord_fixed(1.3) 
mp

ggsave("DCA_sitemap.png",
       plot = mp, device = "png", width = 9, height = 9, units = "in")

