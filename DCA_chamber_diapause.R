pkgs <- c("lubridate", "dplyr", "tidyr", "lme4",
          "stringr", "purrr", "ggplot2", "viridis")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them
source("CDL_funcs.R")
theme_set(theme_bw(base_size = 18)) 

# for table exports
library(stargazer)
library(sjPlot)


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

dat$sitelat <- as.factor(paste(dat$site, round(dat$lat, 1), sep = ": "))

# convert daylength time to numeric for modeling
dat$photoperiod <- sapply(strsplit(as.character(dat$day),":"),
                  function(x) {
                    x <- as.numeric(x)
                    x[1]+x[2]/60
                  }
)

# scaling photoperiod for model fitting
zday <- scale(dat$photoperiod)
zyear <- scale(dat$year)
zlat <- scale(dat$lat)
dat$zday <- as.numeric(zday[, 1])
dat$zyear <- as.numeric(zyear[, 1])
dat$zlat <- as.numeric(zlat[, 1])


# glm with random effects
# response is the percent reproductive, with the total count for the dish as the model weights (could do this other ways)
# trying out different model random effects
# works best with random intercepts and slopes by population (mod1)
mod <- glmer(perc_repro ~ zday + year + (1|rowid), weights = total, data = dat, family = binomial)
mod1 <- glmer(perc_repro ~ (zday + zlat + zyear)^2 + (1 + zday|site), weights = total, data = dat, family = binomial)
mod2 <- glmer(perc_repro ~ zday + zlat + zday:zlat + (1 + zday|site), weights = total, data = dat, family = binomial)
summary(mod)
AIC(mod1, mod2)

# stargazer(mod2, type = "text")
tab_model(mod2, file = "allsite_glmer.html", show.r2 = FALSE, show.icc = FALSE)

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
newdat$photoperiod <- newdat$zday * attr(zday, 'scaled:scale') + attr(zday, 'scaled:center')
newdat$year <- newdat$zyear * attr(zyear, 'scaled:scale') + attr(zyear, 'scaled:center')
newdat$lat <- newdat$zlat * attr(zlat, 'scaled:scale') + attr(zlat, 'scaled:center')
# newdat$lat <- as.numeric(newdat$lat)
# newdat$zlat <- as.numeric(newdat$zlat)
newdat$sitelat <- as.factor(paste(newdat$site, round(newdat$lat, 1), sep = ": "))

cdl <- newdat %>% 
  group_by(site, sitelat) %>% 
  summarise(cdl = as.numeric(photoperiod[which.min(abs(pred - .5))]))

cdl$sitelat <- factor(cdl$sitelat, levels = cdl$sitelat[order(cdl$cdl, decreasing = TRUE)])
newdat$sitelat <- factor(newdat$sitelat, levels = cdl$sitelat[order(cdl$cdl, decreasing = TRUE)])
dat$sitelat <- factor(dat$sitelat, levels = cdl$sitelat[order(cdl$cdl, decreasing = TRUE)])

# plot shows raw data for all chambers and model prediction for 21C chambers by population
plt <- ggplot(dat, aes(x = photoperiod, y = perc_repro)) +
  geom_point(aes(shape = as.factor(year)), alpha = .5, size = 1.6) +
  scale_shape_discrete(name = "Year") +
  geom_line(data = newdat, aes(x = photoperiod, y = pred, color = lat), size = 1.5, alpha = .5) +
  scale_color_viridis(discrete = FALSE, name = "Latitude") +
  facet_wrap(~sitelat, ncol = 3) +
  geom_text(data = cdl, aes(x = cdl + .3, y = .1, label = round(cdl,1))) +
  geom_vline(data = cdl, aes(xintercept = cdl)) +
  xlab("Hours of light in chamber") +
  ylab("Percent reproductive") +
  ggtitle("Critical daylengths of Diorhabda populations") +
  theme(legend.position = c(.65, .075), legend.direction = "horizontal", legend.margin = margin(1, 1, 1, 1, unit = "pt"))
plt

ggsave("DCA_diapause_allyears.png",
       plot = plt, device = "png", width = 9, height = 9, units = "in")


# with just sites moving south through arizona
datriver <- dat %>%
  filter(lat < 37.5 & lon < -110) %>% 
  mutate(year = factor(year, levels = c("2014", "2017")))
zday <- scale(datriver$photoperiod)
zlat <- scale(datriver$lat)
datriver$zday <- as.numeric(zday[, 1])
datriver$zlat <- as.numeric(zlat[, 1])
# no random effects of sites, assumed to be similar/explained by latitude
mod <- glm(perc_repro ~ zday * zlat * year, weights = total, data = datriver, family = binomial)

tab_model(mod, file = "riversite_glm.html", show.r2 = FALSE, show.icc = FALSE)

plot_model(mod)
plot_model(mod, type = "int")
stargazer(mod, type = "text")


# model predictions for plotting
newdat <- expand.grid(list(zday = seq(-4, 2.5, length.out = 500),
                           year = unique(datriver$year),
                           zlat = c(-1.865, -.7285, 1.143)))
# uniqdat <- datriver %>% 
#   dplyr::select(site, zlat) %>% 
#   distinct()
# 
# newdat <- left_join(newdat, uniqdat)
ilink <- family(mod)$linkinv
preds <- predict(mod, newdat, type = "link", se.fit = TRUE)
preds <- transform(preds, Fitted = ilink(fit), Upper = ilink(fit + (2 * se.fit)),
                   Lower = ilink(fit - (2 * se.fit)))
newdat <- cbind(newdat, preds)
# back transform zday to photoperiod
newdat$photoperiod <- newdat$zday * attr(zday, 'scaled:scale') + attr(zday, 'scaled:center')
newdat$lat <- newdat$zlat * attr(zlat, 'scaled:scale') + attr(zlat, 'scaled:center')
newdat <- newdat[-which(newdat$year == 2014 & newdat$lat < 34.2), ]

cdl <- newdat %>% 
  group_by(lat, year) %>% 
  summarise(cdl = as.numeric(photoperiod[which.min(abs(Fitted - .5))]))

cdl$lat <- as.factor(as.character(round(cdl$lat, 1)))
cdl$Fitted <- 0.5
cdl$photoperiod <- cdl$cdl
newdat$lat <- as.factor(as.character(round(newdat$lat, 1)))

# model predictions for different latitudes and 2014 vs 2017
plt <- ggplot(newdat, aes(x = photoperiod, y = Fitted, group = interaction(lat, year))) +
  geom_ribbon(data = newdat, aes(ymin = Lower, ymax = Upper, x = photoperiod, 
                                 fill = lat), alpha = .2, inherit.aes = FALSE) +
  geom_line(data = newdat, aes(x = photoperiod, y = Fitted, color = lat), 
            size = 1.5, alpha = .5) +
  scale_fill_discrete(name = "Latitude") +
  scale_color_discrete(guide = FALSE) +
  facet_wrap(~year, ncol = 1) +
  geom_point(data = cdl, aes(x = photoperiod, y = Fitted)) +
  geom_label_repel(data = cdl, aes(label = round(photoperiod, 1))) +
  xlab("Hours of light in chamber (35/15C)") +
  ylab("Percent reproductive (modeled)") +
  ggtitle("Critical daylengths of Diorhabda moving south")
plt

ggsave("DCA_diapause_south.png",
       plot = plt, device = "png", width = 9, height = 9, units = "in")


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

