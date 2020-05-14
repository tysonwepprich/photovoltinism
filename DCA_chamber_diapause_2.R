# Updating with 2019 data
# using some methods from Galerucella chamber analysis


pkgs <- c("lubridate", "dplyr", "tidyr", "lme4",
          "stringr", "purrr", "ggplot2", "viridis")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them
source("CDL_funcs.R")
'%!in%' <- function(x,y)!('%in%'(x,y))

theme_set(theme_bw(base_size = 20) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())) 

dat <- read.csv("data/DCA_CDL_3515_all.csv", header = TRUE)

# Remove sites with only 1 photo trt
dat <- dat[-which(dat$site %in% c("Parker")),]
dat <- dat[-which(dat$site == "St. George" & dat$year == 2017),]
dat$site <- plyr::revalue(dat$site, c("Virgin River"="Gold Butte"))
                          
dat <- dat %>% droplevels()
levels(dat$site) <- c("Artesia, NM","Big Bend, NV","Blythe, CA", "Cibola, CA",
                      "Delta, UT", "Gold Butte, NV","Imperial, AZ","Lake Mead, NV", 
                      "Lovell, WY", "Lovelock, NV","Princess, NV","Pueblo, CO", 
                      "St. George, UT",  "Topock Marsh, AZ")

# Reorder factors
dat$site <- factor(dat$site, levels = c("Lovell, WY", "Lovelock, NV", "Pueblo, CO", "Artesia, NM",
                      "Delta, UT", "St. George, UT", "Gold Butte, NV", "Lake Mead, NV", "Princess, NV",
                      "Big Bend, NV", "Topock Marsh, AZ", "Blythe, CA", "Cibola, CA",
                      "Imperial, AZ"))

# Other onetime sites removed
dat <- dat[-which(dat$site %in% c("Artesia, NM", "Pueblo, CO")),]
dat <- dat %>% droplevels()

sites <- dat %>% 
  # filter(year == 2019) %>%  
  dplyr::select(site, lat, lon) %>%
  distinct() %>% mutate(Species = "Diorhabda carinulata") %>% 
  group_by(site, Species) %>% 
  summarise(lat = mean(lat), lon = mean(lon))


# some had 3 groups in 2017, group all together for now
dat <- dat %>% 
  group_by(day, site, lat, lon, year) %>% 
  summarise(ndiap = sum(ndiap),
            nrep = sum(nrep))

dat <- dat %>% 
  rowwise() %>% 
  mutate(total = ndiap + nrep,
         perc_repro = nrep/total,
         perc_diap = ndiap/total,
         rowid = paste(trimws(site, which = "both"), year, sep = "_"))  %>%
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
dat$yearfact <- as.factor(as.character(dat$year))
dat$zday <- as.numeric(zday[, 1])
dat$zyear <- as.numeric(zyear[, 1])
dat$zlat <- as.numeric(zlat[, 1])

dat <- dat %>% filter(year == 2019)
# glm with random effects
# response is the percent reproductive, with the total count for the dish as the model weights (could do this other ways)
# trying out different model random effects
# works best with random intercepts and slopes by population (mod1)
mod <- glm(perc_repro ~ zday + site, weights = total, data = dat, family = binomial)
mod1 <- glm(perc_repro ~ zday * site, weights = total, data = dat, family = binomial)
mod2 <- glm(perc_repro ~ photoperiod * site, weights = total, data = dat, family = binomial)

mod3 <- glm(perc_repro ~ rowid * photoperiod, weights = total, data = dat, family = binomial(link = "probit"))
mod3 <- glmer(perc_repro ~ zday + (1 + zday|rowid), weights = total, data = dat, family = binomial(link = "probit"))


summary(mod3)
AIC(mod, mod1, mod2)
car::Anova(mod3, type = 2)

# stargazer(mod2, type = "text")
# tab_model(mod2, file = "allsite_glmer.html", show.r2 = FALSE, show.icc = FALSE)

# Main problem with GLMER is figuring out the CI for the CP. Dose.p.glmm worked for GCA
# but not with varying slopes model here.


# # GLMER model predictions for plotting
# newdat <- expand.grid(list(zday = seq(-2.3, 1.5, length.out = 500),
#                            rowid = unique(dat$rowid)))
# uniqdat <- dat %>% 
#   dplyr::select(site, lat, zlat, year, rowid) %>% 
#   distinct()
# 
# newdat <- left_join(newdat, uniqdat)
# # # GLM model predictions for plotting
# # newdat <- expand.grid(list(photoperiod = seq(10, 16, length.out = 500),
# #                            rowid = unique(dat$rowid)))
# # uniqdat <- dat %>% 
# #   dplyr::select(site, lat, year, rowid) %>% 
# #   distinct()
# # 
# # newdat <- left_join(newdat, uniqdat)
# 
# newdat$pred <- predict(mod3, newdat, type = "response")
# # back transform zday to photoperiod
# newdat$photoperiod <- newdat$zday * attr(zday, 'scaled:scale') + attr(zday, 'scaled:center')
# 
# 
# cdl <- newdat %>% 
#   group_by(site, lat, year, rowid) %>% 
#   summarise(cdl = as.numeric(photoperiod[which.min(abs(pred - .5))]))
# 
# # CP Confidence Interval----
# dose.p.glmm <-  function(obj, cf = c(1, 2), p = .5) {
#   f <- family(obj)
#   eta <- f$linkfun(p)
#   b <- fixef(obj)[cf]
#   x.p <- (eta - b[1L])/b[2L]
#   # names(x.p) <- paste("p = ", format(p), ":", sep = "")
#   pd <- -cbind(1, x.p)/b[2L]
#   SE <- sqrt(((pd %*% vcov(obj)[cf, cf]) * pd) %*% c(1, 1))
#   res <- data.frame(pred = x.p, SE = matrix(SE), p = p)
#   # class(res) <- "glm.dose"
#   res
# }
# 
# test <- list()
# for (i in 1:6){
#   j <- c(7:2)[i]
#   tmp <- dose.p.glmm(mod3, cf = c(j, 1), p = .5)
#   invdose <- c(tmp$pred - 1.96 * tmp$SE, tmp$pred, tmp$pred + 1.96 * tmp$SE)
#   df <- data.frame(t(invdose * attr(alldat$zday, 'scaled:scale') + attr(alldat$zday, 'scaled:center')))
#   names(df) <- c("lwr", "cp", "upr")
#   df$population <- levels(alldat$population)[j-1]
#   test[[length(test)+1]] <- df
#   }
# 
# library(merTools)
# exampPreds <- predictInterval(mod3, newdata = newdat, which = "full",
#                               type = "probability", n.sims = 100, returnSims = TRUE,
#                               include.resid.var = FALSE, level = 0.95)
# newdat <- bind_cols(newdat, exampPreds)
# newdat$row <- 1:nrow(newdat)
# 
# cdl <- newdat %>%
#   group_by(rowid) %>%
#   summarise(cdl = as.numeric(photoperiod[which.min(abs(fit - .5))]),
#             row = as.numeric(row[which.min(abs(fit - .5))]))
# 
# cdl_ci <- apply(attr(exampPreds, "sim.results"), MARGIN = 1, FUN = quantile, probs = c(0.05, .5, .95))

library(broom)
# Run all models at once for site x year
mods <- dat %>% 
  group_by(rowid) %>% 
  do(fit = glm(perc_diap ~ photoperiod, weights = total, data = ., family = binomial(link = "probit")))

modcoefs <- bind_rows(lapply(mods$fit, FUN = function(x) data.frame(b0 = coef(x)[1], b1 = coef(x)[2]))) %>% 
  mutate(mods$rowid) %>% 
  mutate(cp_mean = -b0/b1,
         cp_sd = 1/abs(b1))
modcps <- bind_rows(lapply(mods$fit, FUN = function(x) {
  tmp <- MASS::dose.p(x, p = c(.05, .5, .95))
  data.frame(pred = as.vector(tmp),
             se = as.vector(attr(tmp, "SE")),
             p = as.vector(attr(tmp, "p")))
                      })) %>% 
  mutate(rowid = rep(mods$rowid, each = 3)) %>% 
  pivot_wider(names_from = p, values_from = c(pred, se))

tabout <- modcps %>% 
  mutate(Site = str_split_fixed(rowid, "_", 2)[,1],
         Year = str_split_fixed(rowid, "_", 2)[,2],
         CP = pred_0.5,
         SE = se_0.5) %>% 
  left_join(distinct(dat[,c("rowid", "lat")])) %>% 
  dplyr::select(Site, lat, Year, CP, SE) %>% 
    arrange(-lat)
write.csv(tabout, file = "DCA_CP.csv", row.names = FALSE)

# GLM model predictions for plotting
newdat <- expand.grid(list(photoperiod = seq(10, 16, length.out = 500),
                           rowid = unique(dat$rowid)))
uniqdat <- dat %>%
  dplyr::select(site, lat, year, rowid) %>%
  distinct()

newdat <- left_join(newdat, uniqdat)
newdat$pred <- predict(mod3, newdat, type = "response")

cdl <- newdat %>%
  group_by(rowid) %>%
  summarise(cdl = as.numeric(photoperiod[which.min(abs(pred - .5))])) %>% 
  left_join(uniqdat)

# Gradient of CPs
sites <- dat %>% dplyr::select(site, lat, lon, year, rowid) %>% distinct()

dat <- left_join(modcps, sites)
dat$year <- as.factor(dat$year)
dat$lat <- ifelse(dat$site %in% c("Delta, UT"),
                  ifelse(dat$year == "2014", dat$lat + 0.1,
                         ifelse(dat$year == "2019", dat$lat - 0.1, dat$lat)), dat$lat)
photo <- data.frame(lat = seq(33, 45, by = 0.01),
                    dl = photoperiod(seq(33, 45, by = 0.01), 172, p = 1.5))

# make y-axis transform to match mercator?
MercLats <- function(lats){
  latRad <-  lats*pi/180
  mercN <- log(tan((pi/4)+(latRad/2)))
  y     <- mercN/(2*pi)
}
photo$lat2 <- MercLats(photo$lat)
dat$lat2 <- MercLats(dat$lat)


dat$upr <- 1.96*dat$se_0.5 + dat$pred_0.5
dat$lwr <- dat$pred_0.5 - 1.96*dat$se_0.5 

dat <- dat %>% filter(lat < 40) %>%
  filter(rowid != "Gold Butte, NV_2017") %>% 
  droplevels()

plt <- ggplot(dat, aes(x = lat, y = cp_mean, shape = year)) +
  geom_point(size = 5) +
  scale_shape_discrete(name = "Critical\nphotoperiod") +
  geom_linerange(aes(ymin = lwr, ymax = upr), size = 1.5, alpha = .6) +
  # annotate("text", label = "solstice", x = 38.5, y = 14.95, size = 5, angle = 18, fontface = "italic") +
  xlab("Latitude") +
  ylab("Day length (hours)") +
  theme(legend.position = c(.87, .1))
plt

# flipped
library(ggstance)
plt1 <- ggplot(dat, aes(x = pred_0.5, y = lat2, shape = year)) +
  geom_point(size = 5) +
  scale_shape_discrete(name = "Critical\nphotoperiod") +
  geom_linerangeh(aes(xmin = lwr, xmax = upr), size = 1.5, alpha = .6) +
  ylab("Latitude") +
  xlab("Day length (hours)") +
  scale_y_continuous(limits = MercLats(c(32, 40)), expand = c(0, 0)) +
  theme(legend.position = c(.225, .775))
plt1


library(ggpubr)
p <- ggarrange(plt, plt1 + rremove("ylab") + rremove("y.ticks") + rremove("y.text"),
               ncol = 2, nrow = 1, align = "h")
p

library(ggrepel)
# plot shows raw data for all chambers and model prediction for 21C chambers by population
plt <- ggplot(dat, aes(x = photoperiod, y = perc_repro)) +
  geom_point(aes(shape = as.factor(year)), alpha = .5, size = 1.6) +
  scale_shape_discrete(name = "Year") +
  geom_line(data = newdat, aes(x = photoperiod, y = pred, color = as.factor(year), group = as.factor(year)), size = 1, alpha = .8) +
  scale_color_brewer(palette = "Dark2", name = "Year") +
  facet_wrap(~site, ncol = 3) +
  geom_text_repel(data = cdl, aes(x = cdl, y = .5, label = round(cdl,1), color = as.factor(year))) +
  geom_hline(yintercept = .5, linetype = "dotted") +
  xlab("Hours of light in chamber") +
  ylab("Percent reproductive") +
  ggtitle("Critical daylengths of Diorhabda populations")
  # theme(legend.position = c(.65, .075), legend.direction = "horizontal", legend.margin = margin(1, 1, 1, 1, unit = "pt"))
plt

ggsave("DCA_diapause_allyears.png",
       plot = plt, device = "png", width = 9, height = 9, units = "in")


# One panel 2019
plt <- ggplot(dat, aes(x = photoperiod, y = perc_repro, group = site)) +
  geom_point(aes(color = site), alpha = .5, size = 1.6) +
  # scale_shape_discrete(name = "Site") +
  geom_line(data = newdat, aes(x = photoperiod, y = pred, color = site, group = site), size = 1, alpha = .8) +
  scale_color_brewer(palette = "Dark2", name = "Population") +
  geom_label_repel(data = cdl, aes(x = cdl, y = .5, label = round(cdl,1), color = site), label.padding = .15) +
  # geom_hline(yintercept = .5, linetype = "dotted") +
  xlab("Hours of light in chamber") +
  ylab("Percent reproductive") +
  ggtitle("Critical daylengths of Diorhabda populations")
# theme(legend.position = c(.65, .075), legend.direction = "horizontal", legend.margin = margin(1, 1, 1, 1, unit = "pt"))
plt


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
  coord_fixed(1.3) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
mp

ggsave("DCA_sitemap2.png",
       plot = mp, device = "png", width = 9, height = 9, units = "in")






