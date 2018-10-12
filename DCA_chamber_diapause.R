pkgs <- c("lubridate", "dplyr", "tidyr", "lme4",
          "stringr", "purrr", "ggplot2", "viridis")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them
source("CDL_funcs.R")
theme_set(theme_bw(base_size = 18)) 


dat <- read.csv("data/DCA_CDL_3515_all.csv", header = TRUE)

dat <- dat %>% 
  rowwise() %>% 
  mutate(total = ndiap + nrep,
         perc_repro = nrep/total,
         rowid = paste(trimws(site, which = "both"), year, group, sep = "_"))  %>%
  droplevels.data.frame() %>% 
  group_by(site) %>% 
  mutate(lat = mean(lat), 
         lon = mean(lon))

dat$sitelat <- paste(dat$site, round(dat$lat, 1), sep = ": ")

# convert daylength time to numeric for modeling
dat$photoperiod <- sapply(strsplit(as.character(dat$day),":"),
                  function(x) {
                    x <- as.numeric(x)
                    x[1]+x[2]/60
                  }
)

# scaling photoperiod for model fitting
dat$zday <- scale(dat$photoperiod)
dat$zyear <- scale(dat$year)
dat$zlat <- scale(dat$lat)

# glm with random effects
# response is the percent reproductive, with the total count for the dish as the model weights (could do this other ways)
# trying out different model random effects
# works best with random intercepts and slopes by population (mod1)
mod <- glmer(perc_repro ~ zday + year + (1|rowid), weights = total, data = dat, family = binomial)
mod1 <- glmer(perc_repro ~ (zday + zlat + zyear)^2 + (1 + zday|site), weights = total, data = dat, family = binomial)
mod2 <- glmer(perc_repro ~ zday + zlat + zday:zlat + (1 + zday|site), weights = total, data = dat, family = binomial)
summary(mod)
AIC(mod1, mod2)


# with just sites moving south through arizona
datriver <- dat %>%
  filter(lat < 37.5 & lon < -110)
# no random effects of sites, assumed to be similar/explained by latitude
mod1 <- glm(perc_repro ~ (zday + zlat + zyear)^2, weights = total, data = datriver, family = binomial)



# model predictions for plotting
newdat <- expand.grid(list(zday = seq(-2.5, 2.5, length.out = 100),
                           zyear = unique(dat$zyear),
                           site = unique(dat$site)))
uniqdat <- dat %>% 
  dplyr::select(site, zlat) %>% 
  distinct()

newdat <- left_join(newdat, uniqdat)

newdat$pred <- predict(mod2, newdat, type = "response")
# back transform zday to photoperiod
newdat$photoperiod <- newdat$zday * attr(dat$zday, 'scaled:scale') + attr(dat$zday, 'scaled:center')
newdat$year <- newdat$zyear * attr(dat$zyear, 'scaled:scale') + attr(dat$zyear, 'scaled:center')
newdat$lat <- newdat$zlat * attr(dat$zlat, 'scaled:scale') + attr(dat$zlat, 'scaled:center')



# plot shows raw data for all chambers and model prediction for 21C chambers by population
plt <- ggplot(dat, aes(x = photoperiod, y = perc_repro)) +
  geom_point(aes(shape = as.factor(year)), alpha = .5) +
  geom_line(data = newdat, aes(x = photoperiod, y = pred, color = zlat), size = 1.5) +
  scale_color_viridis(discrete = FALSE) +
  facet_wrap(~site, ncol = 3)
plt


# plot shows raw data for all chambers and model prediction for 21C chambers by population
plt <- ggplot(newdat, aes(x = photoperiod, y = pred, group = site, color = zlat)) +
  geom_line(size = 1.5) +
  scale_color_viridis(discrete = FALSE)
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
  geom_text_repel(aes(label = site)) +
  geom_polygon(data = states, aes(x = x, y = y, group = group), fill = NA, color = "black", size = .1) +
  coord_fixed(1.3) 
mp

