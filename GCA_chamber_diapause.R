pkgs <- c("lubridate", "dplyr", "tidyr", "lme4", "lmerTest",
          "stringr", "purrr", "ggplot2", "viridis")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load them

'%!in%' <- function(x,y)!('%in%'(x,y))

theme_set(theme_bw(base_size = 20) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())) 

# Prepare data ####
# 2014/2019 data for photoperiod response

# Previously had included 2018, but temp different at 21/26, no chamber temp calibration
# dat2018 <- read.csv("data/galerucella_chamber_diapause_2018.csv", header = TRUE) %>% 
#   mutate(Year = 2018)
dat2019 <- read.csv("data/galerucella_chamber_diapause_2019.csv", header = TRUE) %>% 
  mutate(Year = 2019)

dat <- dat2019 %>% 
  rowwise() %>% 
  mutate(total = repro_female + feeder + diapause,
         perc_repro = 1 - diapause/total,
         rowid = paste(population, chamber, emerg_group, sep = "_")) %>% 
  droplevels.data.frame()
levels(dat$population) <- c("Bellingham, WA", "Baskett-Slough, OR", 
                            "McArthur, CA", "Palermo, CA", "Sutherlin, OR",
                            "Yakima, WA")
dat$population <- plyr::revalue(dat$population, c("Yakima, WA"="Yakima TC, WA", 
                                                  "Baskett-Slough, OR"="Rickreall, OR"))

# From Fritzi's experiments, keep 2014 but not 2013 (no gen in common env first as control)
dat2 <- read.csv('data/combined 13_15.csv')
dat2$population <- dat2$Site
levels(dat2$population) <- c("Bellingham, WA", "Ephrata, WA", "Palermo, CA", "Sutherlin, OR")

alldat <- dat2 %>% 
  filter(Year == 2014) %>% 
  group_by(population, Treatment, Year) %>% 
  summarise(repro = sum(diapause),
            total = n(),
            perc_repro = repro/total,
            photoperiod = Treatment[1],
            temperature = 23,
            rowid = paste(population[1], Treatment[1], sep = "_")) %>% 
  bind_rows(dat)

# Fix factor orders

alldat$population <- factor(alldat$population, levels = c("Bellingham, WA", "Ephrata, WA", "Yakima TC, WA", "Rickreall, OR",
                                                    "Sutherlin, OR", "McArthur, CA", "Palermo, CA"))


# scaling photoperiod for model fitting
alldat$zday <- scale(alldat$photoperiod)



# Fit GLM ----

# GLM vs GLMER????
# options(contrasts = c("contr.treatment", "contr.poly"))
alldat <- alldat %>% filter(Year == 2019) %>% droplevels()
# Using varying slopes model had singular fit, 
mod0 <- glmer(perc_repro ~ zday +  (1 + zday| population) + (1|population:rowid), weights = total, 
              data = alldat, family = binomial(link = "logit"))
# assuming population random effects, rowid nested within population
mod1 <- glmer(perc_repro ~ zday +  (1 | population) + (1|population:rowid), weights = total, 
              data = alldat, family = binomial(link = "logit"))
# Population as fixed effect, petri dish random to control for pseudoreplication
mod2 <- glmer(perc_repro ~ zday + (population-1) + emerg_group +  (1|rowid), weights = total, 
              data = alldat, family = binomial(link = "logit"), glmerControl(optimizer = "bobyqa"))
mod3 <- glmer(perc_repro ~ zday + (population-1) + (1 | rowid), weights = total, 
              data = alldat, family = binomial(link = "logit"))
mod4 <- glmer(perc_repro ~ zday * (population) + (1 | rowid), weights = total, 
              data = alldat, family = binomial(link = "logit"), glmerControl(optimizer = "bobyqa"))

# 2014 AIC, Chisq says to drop interaction term
# 2019 same, also emerg_group not important
AIC(mod0, mod1, mod2, mod3, mod4)
drop1(mod3, .~., test = "Chisq")
anova(mod4, mod3)
summary(mod4)
MuMIn::r.squaredGLMM(mod3)
car::Anova(mod4, type = 2)

sjPlot::tab_model(mod3a, mod3, transform = NULL, show.ci = FALSE, show.se = TRUE,
                  show.est = TRUE, show.stat = TRUE, show.r2 = TRUE)


# GLM overdispersion, justifies adding rowid to 2014 data
mod <- glm(perc_repro ~ zday + (population-1) , weights = total, 
           data = alldat, family = binomial(link = "logit"))
mod1 <- glm(perc_repro ~ zday + (population-1) , weights = total, 
            data = alldat, family = quasibinomial(link = "logit"))

pchisq(summary(mod1)$dispersion * mod$df.residual, mod$df.residual, lower = F) # significance for overdispersion

library(merTools)




## grab the inverse link function
ilink <- family(mod3)$linkinv

# model predictions for plotting
newdat <- expand.grid(list(zday = seq(-1.8, 2, length.out = 500),
                           population = unique(alldat$population),
                           Year = unique(alldat$Year),
                           rowid = unique(alldat$rowid)[1]))

# seems to only work on model with intercept
mod3 <- glmer(perc_repro ~ zday + population + (1 | rowid), weights = total, 
              data = alldat, family = binomial(link = "logit"))
exampPreds <- predictInterval(mod3, newdata = newdat, which = "fixed",
                              type = "probability",
                              include.resid.var = FALSE, level = 0.95)
newdat <- bind_cols(newdat, exampPreds)
# GLM confidence intervals
# ## add fit and se.fit on the **link** scale
# newdat <- bind_cols(newdat, setNames(as_tibble(predict(mod3, newdat, re.form = NA, se.fit = TRUE)[1:2]),
#                                    c('fit_link','se_link')))
# ## create the interval and backtransform
# newdat <- mutate(newdat,
#                 fit_resp  = ilink(fit_link),
#                 right_upr = ilink(fit_link + (1.96 * se_link)),
#                 right_lwr = ilink(fit_link - (1.96 * se_link)))# back transform zday to photoperiod

newdat$photoperiod <- newdat$zday * attr(alldat$zday, 'scaled:scale') + attr(alldat$zday, 'scaled:center')

saveRDS(newdat, "GCpredictions.rds")

# cdl <- newdat %>% 
  # group_by(population) %>% 
  # summarise(cdl = as.numeric(photoperiod[which.min(abs(fit_resp - .5))]))

# CP Confidence Interval----
# dose.p.glmm <-  function(obj, cf = c(7, 1), p = .5) {
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
# cdl2014 <- bind_rows(test) %>% 
#   mutate(year = 2014)
# cdl2019 <- bind_rows(test) %>% 
#   mutate(year = 2019)
# cdls <- round_df(bind_rows(cdl2014, cdl2019), 2)
# write.csv(cdls, file = "cps.csv", row.names = FALSE)

cdls <- read.csv(file = "cps.csv")


# Test how pairwise populations are different from each other
library(emmeans)
cells <- emmeans(mod3, pairwise ~ population, adjust = "tukey")
pp <- pwpp(cells$emmeans, type = "response", method = "pairwise", values = TRUE)
pp + ylab("Population")
CLD(cells)

emm1 <- emmeans(mod3, specs = ~ population, adjust = "tukey")

North2014 <- c(1, 1, 0, 0)/2
South2014 <- c(0, 0, 1, 1)/2
North2019 <- c(1, 1, 1, 0, 0, 0)/3
South2019 <- c(0, 0, 0, 1, 1, 1)/3
contrast(emm1, method = list(South2014 - North2014), type = "response", ratios = TRUE)
contrast(emm1, method = list(South2019 - North2019), type = "response", ratios = FALSE)


# Photoperiod figures -----
# plot shows raw data for all chambers and model prediction for 21C chambers by population
plt <- ggplot(alldat, aes(x = photoperiod, y = perc_repro)) +
  geom_point(aes(size = total), alpha = 0.5) +
  scale_size_continuous(name = "Sample size", range = c(1, 4), breaks = c(10, 20, 30, 40)) +
  scale_x_continuous(limits = c(14.25, 17.25), expand = c(0, 0), breaks = seq(14.5, 17, 0.5)) +
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  geom_line(data = newdat, aes(x = photoperiod, y = fit)) +
  geom_ribbon(data = newdat, aes(x = photoperiod, ymin = lwr, ymax = upr, group = as.factor(Year)), alpha = 0.1, inherit.aes = FALSE) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  # geom_point(data = cdl, aes(x = cp, y = 0.5), size = 3, shape = 17) +
  # geom_errorbarh(data = cdl, aes(xmin = lwr, y = 0.5, xmax = upr), size = .8, height = .025, inherit.aes = FALSE) +
  ylab("Proportion reproductive") +
  xlab("Treatment: day length (hours)") +
  facet_wrap(~population, ncol = 2) 
#   theme(legend.position = c(0.1, .9), 
#         legend.text = element_text(size = rel(.5)),
#         legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
# )
plt
ggsave(filename = "GCA_GC_2019.png", plot = plt, device = "png", width = 8, height = 8, units = "in")


# Facets, single point, cp highlighted
alldat <- alldat %>% 
  group_by(population, photoperiod) %>% 
  summarise(perc_repro = weighted.mean(perc_repro, total, na.rm = TRUE),
            n = sum(total, na.rm = TRUE))

ns <- alldat %>% group_by(population) %>% summarise(n = paste0("n=", sum(n, na.rm = TRUE)))
cdls <- cdls %>% filter(year == 2014) %>% droplevels()

plt <- ggplot(alldat, aes(x = photoperiod, y = perc_repro)) +
  geom_point(alpha = 0.5) +
  scale_x_continuous(limits = c(14.25, 17.25), expand = c(0, 0), breaks = seq(14.5, 17, 0.5)) +
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  geom_line(data = newdat, aes(x = photoperiod, y = fit)) +
  geom_ribbon(data = newdat, aes(x = photoperiod, ymin = lwr, ymax = upr, group = as.factor(Year)), alpha = 0.1, inherit.aes = FALSE) + 
  geom_segment(data = cdls, aes(x = 14.25, y = 0.5, xend = cp, yend = 0.5), linetype = "dotted", inherit.aes = FALSE) +
  geom_segment(data = cdls, aes(x = cp, y = 0, xend = cp, yend = 0.5), linetype = "dotted", inherit.aes = FALSE) +
  geom_label(data = cdls, aes(x = cp, y = 0.05, label = round(cp, 1)), inherit.aes = FALSE) +
  geom_text(data = ns, aes(x = 14.5, y = .95, label = n), inherit.aes = FALSE) +
    # geom_point(data = cdls, aes(x = cp, y = 0.5), size = 3, shape = 17) +
  # geom_errorbarh(data = cdls, aes(xmin = lwr, y = 0.5, xmax = upr), size = .8, height = .025, inherit.aes = FALSE) +
  ylab("Proportion reproductive") +
  xlab("Photoperiod (light hours)") +
  facet_wrap(~population, ncol = 2) 
#   theme(legend.position = c(0.1, .9), 
#         legend.text = element_text(size = rel(.5)),
#         legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
# )
plt
ggsave(filename = "fig2.tif", plot = plt, device = "tiff", width = 8, height = 8, units = "in", dpi = 600)
ggsave(filename = "fig1.tif", plot = plt, device = "tiff", width = 8, height = 6, units = "in", dpi = 600)


# No Facets
alldat <- alldat %>% 
  group_by(population, photoperiod) %>% 
  summarise(perc_repro = weighted.mean(perc_repro, total))

plt <- ggplot(alldat, aes(x = photoperiod, y = perc_repro, group = population, shape = population)) +
  # geom_point(aes(size = total), alpha = 0.5) +
  geom_point(alpha = 0.5, size = 4) +
  scale_shape_manual(values=c(0, 1, 2, 4, 5, 6))+
  scale_size_continuous(name = "Sample size", range = c(1, 4), breaks = c(10, 20, 30, 40)) +
  scale_x_continuous(limits = c(14.25, 17.25), expand = c(0, 0), breaks = seq(14.5, 17, 0.5)) +
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  geom_line(data = newdat, aes(x = photoperiod, y = fit, linetype = population), alpha = 0.5, size = 1) +
  # scale_linetype_manual(values=c(1, 5, 3, 4)) +
  # geom_ribbon(data = newdat, aes(x = photoperiod, ymin = lwr, ymax = upr, group = population), alpha = 0.1, inherit.aes = FALSE) + 
  # geom_hline(yintercept = 0.5, linetype = "dashed") +
  # geom_point(data = cdl, aes(x = cp, y = 0.5), size = 3, shape = 17) +
  # geom_errorbarh(data = cdl, aes(xmin = lwr, y = 0.5, xmax = upr), size = .8, height = .025, inherit.aes = FALSE) +
  ylab("Proportion reproductive") +
  xlab("Treatment: photoperiod (hours)")
  # theme(legend.position = c(.2, .8)) 
  
plt
ggsave(filename = "GCA_GC_2019_1panel.png", plot = plt, device = "png", width = 10, height = 8, units = "in")


# Gradient of CPs
sites <- data.frame(population = c(
  "Ephrata, WA", "Yakima TC, WA",
  "Sutherlin, OR", "Bellingham, WA",
  "McArthur, CA", "Palermo, CA", "Rickreall, OR"),
  long = c(-119.655253, -119.986530,
           -123.315854, -122.479482,
           -121.41, -121.58, -123.27),
  lat = c(47.161647, 46.756318,
          43.387721, 48.756105,
          41.10, 39.41, 44.98))

dat <- left_join(cdls, sites)
dat$year <- as.factor(dat$year)
dat$lat <- ifelse(dat$population %in% c("Palermo, CA", "Sutherlin, OR", "Bellingham, WA"),
  ifelse(dat$year == "2014", dat$lat + 0.09, dat$lat - 0.09), dat$lat)
photo <- data.frame(lat = seq(38.7, 49.1, by = 0.01),
                    dl = photoperiod(seq(38.7, 49.1, by = 0.01), 172, p = 1.5))

# make y-axis transform to match mercator?
MercLats <- function(lats){
latRad <-  lats*pi/180
mercN <- log(tan((pi/4)+(latRad/2)))
y     <- mercN/(2*pi)
}
photo$lat2 <- MercLats(photo$lat)
dat$lat2 <- MercLats(dat$lat)

plt <- ggplot(dat, aes(x = lat, y = cp, shape = year)) +
  geom_point(size = 5) +
  scale_shape_discrete(name = "Critical\nphotoperiod") +
  geom_linerange(aes(ymin = lwr, ymax = upr), size = 1.5, alpha = .6) +
  geom_line(data = photo, aes(x = lat, y = dl), inherit.aes = FALSE) +
  annotate("text", label = "solstice", x = 38.5, y = 14.95, size = 5, angle = 18, fontface = "italic") +
  xlab("Latitude") +
  ylab("Day length (hours)") +
  theme(legend.position = c(.87, .1))
plt

# flipped
library(ggstance)
plt1 <- ggplot(dat, aes(x = cp, y = lat2, shape = year)) +
  geom_point(size = 5) +
  scale_shape_discrete(name = NULL) +
  geom_linerangeh(aes(xmin = lwr, xmax = upr), size = 1.5, alpha = .6) +
  geom_line(data = photo, aes(x = dl, y = lat2), inherit.aes = FALSE, alpha = 0.5) +
  annotate("text", label = "solstice", x = 15.65, y = .14, size = 5, angle = 74, fontface = "italic", alpha = 0.5) +
  # annotate("text", label = "reproductive\nmultivoltine", x = 38, y = 14.5, size = 5, angle = 0) +
  # annotate("text", label = "diapause\nunivoltine", x = 38, y = 15.5, size = 5, angle = 0) +
  ylab("Latitude") +
  xlab("Critical photoperiod (day hours)") +
  scale_y_continuous(limits = MercLats(c(38.7, 49.1)), expand = c(0, 0)) +
  theme(legend.position = c(.175, .925))
plt1


library(ggpubr)
p <- ggarrange(plt, plt1 + rremove("ylab") + rremove("y.ticks") + rremove("y.text"),
                ncol = 2, nrow = 1, align = "h")
p

ggsave(filename = "GCA_GC_map.tif", plot = p, device = "tiff", dpi = 600, width=10.5, height = 7.5, units = "in")


# Sex ratio issue ----
# alternative smooth based on GLM
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

gam_smooth <- function(...) {
  geom_smooth(method = "gam", method.args = list(family = "binomial"), ...)
}

alldat <- alldat %>% filter(Year == 2019) %>% droplevels()
# if 50/50 sex ratio assumed
alldat <- alldat %>% 
  mutate(total_repro = repro_female + feeder,
    sex_ratio = repro_female / total_repro,
         fem_tot = ceiling(total / 2),
         fem_diap = fem_tot - repro_female) %>% 
  mutate(fem_diap = ifelse(fem_diap < 0, 0, fem_diap),
         perc_repro_fem = repro_female / (repro_female + fem_diap),
         fem_tot = repro_female + fem_diap)

# 50/50 sex ratio approached when all reproductive
# implies females more conservative, stay in diapause at low % reproductive
mod <- gam(sex_ratio ~ s(perc_repro, k = 5), weights = total, family = binomial, data = alldat)

plt <- ggplot(alldat, aes(x = perc_repro, y = sex_ratio)) + 
  geom_point(aes(size = log(total)), alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  gam_smooth(aes(weight = total), formula = y ~ s(x, k = 4), se = TRUE) +
  ylab("Female proportion of reproductive adults") +
  xlab("Total proportion of adults reproductive") +
  theme(legend.position = c(.875, .175))
plt


mod3 <- glmer(perc_repro_fem ~ zday + (population-1) + (1 | rowid), weights = fem_tot, 
              data = alldat, family = binomial(link = "logit"), glmerControl(optimizer = "bobyqa"))



## grab the inverse link function
ilink <- family(mod3)$linkinv

# model predictions for plotting
newdat <- expand.grid(list(zday = seq(-1.8, 2, length.out = 500),
                           population = unique(alldat$population),
                           Year = unique(alldat$Year),
                           rowid = unique(alldat$rowid)[1]))

# seems to only work on model with intercept
mod3 <- glmer(perc_repro_fem ~ zday + population + (1 | rowid), weights = fem_tot, 
              data = alldat, family = binomial(link = "logit"), glmerControl(optimizer = "bobyqa"))
exampPreds <- predictInterval(mod3, newdata = newdat, which = "fixed",
                              type = "probability",
                              include.resid.var = FALSE, level = 0.95)
newdat <- bind_cols(newdat, exampPreds)
# GLM confidence intervals
# ## add fit and se.fit on the **link** scale
# newdat <- bind_cols(newdat, setNames(as_tibble(predict(mod3, newdat, re.form = NA, se.fit = TRUE)[1:2]),
#                                    c('fit_link','se_link')))
# ## create the interval and backtransform
# newdat <- mutate(newdat,
#                 fit_resp  = ilink(fit_link),
#                 right_upr = ilink(fit_link + (1.96 * se_link)),
#                 right_lwr = ilink(fit_link - (1.96 * se_link)))# back transform zday to photoperiod

newdat$photoperiod <- newdat$zday * attr(alldat$zday, 'scaled:scale') + attr(alldat$zday, 'scaled:center')

# cdl <- newdat %>% 
# group_by(population) %>% 
# summarise(cdl = as.numeric(photoperiod[which.min(abs(fit_resp - .5))]))

dose.p.glmm <-  function(obj, cf = c(7, 1), p = .5) {
  f <- family(obj)
  eta <- f$linkfun(p)
  b <- fixef(obj)[cf]
  x.p <- (eta - b[1L])/b[2L]
  # names(x.p) <- paste("p = ", format(p), ":", sep = "")
  pd <- -cbind(1, x.p)/b[2L]
  SE <- sqrt(((pd %*% vcov(obj)[cf, cf]) * pd) %*% c(1, 1))
  res <- data.frame(pred = x.p, SE = matrix(SE), p = p)
  # class(res) <- "glm.dose"
  res
}

test <- list()
for (i in 1:6){
  j <- c(7:2)[i]
  tmp <- dose.p.glmm(mod3, cf = c(j, 1), p = .5)
  invdose <- c(tmp$pred - 1.96 * tmp$SE, tmp$pred, tmp$pred + 1.96 * tmp$SE)
  df <- data.frame(t(invdose * attr(alldat$zday, 'scaled:scale') + attr(alldat$zday, 'scaled:center')))
  names(df) <- c("lwr", "cp", "upr")
  df$population <- levels(alldat$population)[j-1]
  test[[length(test)+1]] <- df
  }
# cdl2014 <- bind_rows(test) %>% 
#   mutate(year = 2014)
cdl2019 <- bind_rows(test) %>%
  mutate(year = 2019)
# cdls <- round_df(bind_rows(cdl2014, cdl2019), 2)
# write.csv(cdls, file = "cps.csv", row.names = FALSE)

cdls <- read.csv(file = "cps.csv")

# Test how pairwise populations are different from each other
library(emmeans)
cells <- emmeans(mod3, pairwise ~ population, adjust = "tukey")
pp <- pwpp(cells$emmeans, type = "response", method = "pairwise", values = TRUE)
pp + ylab("Population")


# Photoperiod figures -----
# plot shows raw data for all chambers and model prediction for 21C chambers by population
plt <- ggplot(alldat, aes(x = photoperiod, y = perc_repro_fem)) +
  geom_point(aes(size = fem_tot), alpha = 0.5) +
  scale_size_continuous(name = "Sample size", range = c(1, 4), breaks = c(10, 20, 30, 40)) +
  scale_x_continuous(limits = c(14.25, 17.25), expand = c(0, 0), breaks = seq(14.5, 17, 0.5)) +
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  geom_line(data = newdat, aes(x = photoperiod, y = fit)) +
  geom_ribbon(data = newdat, aes(x = photoperiod, ymin = lwr, ymax = upr, group = as.factor(Year)), alpha = 0.1, inherit.aes = FALSE) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  # geom_point(data = cdl, aes(x = cp, y = 0.5), size = 3, shape = 17) +
  # geom_errorbarh(data = cdl, aes(xmin = lwr, y = 0.5, xmax = upr), size = .8, height = .025, inherit.aes = FALSE) +
  ylab("Proportion reproductive") +
  xlab("Treatment: day length (hours)") +
  facet_wrap(~population, ncol = 2) 
#   theme(legend.position = c(0.1, .9), 
#         legend.text = element_text(size = rel(.5)),
#         legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
# )
plt


# GLMER with all data ----
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

# just modeling selected temperature here
sel_temp <- 23
sel_year <- 2019
# sel_year <- c(2014, 2019)
alldat <- alldat %>% filter(temperature %in% sel_temp,
                             Year %in% sel_year) #%>% 
# filter(population %!in% c("McArthur, CA", "Rickreall, OR")) %>% 
  # mutate(rowid = paste(population, Treatment, sep = "_"))




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
# Site/year analyzed separately ----
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
newdat <- expand.grid(list(zday = seq(-2.5, 2.5, length.out = 200),
                           ztemp = 0.431,
                           population = unique(alldat$population),
                           Year = NA))
newdat$pred <- predict(mod4, newdat, re.form = ~(1 + zday|population), type = "response")
# back transform zday to photoperiod
newdat$photoperiod <- newdat$zday * attr(alldat$zday, 'scaled:scale') + attr(alldat$zday, 'scaled:center')
newdat$temperature <- newdat$ztemp * attr(alldat$ztemp, 'scaled:scale') + attr(alldat$ztemp, 'scaled:center')

# model predictions for plotting each Year/temperature separate lines
linedat <- expand.grid(list(zday = seq(-2.5, 2.5, length.out = 200),
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

# ENTSOC figures ----

lats <- data.frame(population = c("Bellingham, WA", "Vantage, WA", "Rickreall, OR",
                                  "Sutherlin, OR", "McArthur, CA", "Palermo, CA"), 
                   sitelat = c(48.8, 46.7, 45.0, 43.4, 41.1, 39.4),
                   maxhour = photoperiod(c(48.8, 46.7, 45.0, 43.4, 41.1, 39.4), doy = 173),
                   corvhour = photoperiod(44.56, 173))

# scaling photoperiod for model fitting
dat23 <- alldat %>% filter(temperature == 23) %>% droplevels.data.frame()

dat23$zday <- scale(dat23$photoperiod)
dat23$population <- plyr::revalue(dat23$population, c("Yakima, WA"="Vantage, WA", 
                            "Baskett-Slough, OR"="Rickreall, OR"))

dat23$population <- factor(dat23$population, levels = c("Bellingham, WA", "Vantage, WA", "Rickreall, OR",
                                                        "Sutherlin, OR", "McArthur, CA", "Palermo, CA"))


mod2 <- glmer(perc_repro ~ zday + (1|population) + (1|Year), weights = total, data = dat23, family = binomial)

# model predictions for plotting
newdat <- expand.grid(list(zday = seq(-2.5, 2.5, length.out = 100),
                           population = unique(dat23$population),
                           Year = NA))
newdat$pred <- predict(mod2, newdat, re.form = ~(1 |population), type = "response")
# back transform zday to photoperiod
newdat$photoperiod <- newdat$zday * attr(dat23$zday, 'scaled:scale') + attr(dat23$zday, 'scaled:center')

# plot shows raw data for all chambers and model prediction for 21C chambers by population
plt <- ggplot(dat23, aes(x = photoperiod, y = perc_repro)) +
  geom_point(aes(color = as.factor(Year), shape = as.factor(temperature)), size = 4, alpha = .5) +
  # scale_color_viridis(discrete = TRUE) +
  geom_line(data = newdat, aes(x = photoperiod, y = pred), size = 1.5) +
  # geom_line(data = linedat, aes(x = photoperiod, y = pred,
  #                               group = interaction(Year, temperature),
  #                               color = as.factor(Year))) +
  scale_color_discrete(name = "Year") +
  scale_shape_discrete(name = "Temperature") +
  # geom_line(data = newdat, aes(x = photoperiod, y = pred)) +
  xlab("Daylength in chamber") +
  ylab("Percent reproductive") +
  facet_wrap(~population, ncol = 3)
plt

cdl <- newdat %>% 
  group_by(population) %>% 
  summarise(cdl = as.numeric(photoperiod[which.min(abs(pred - .5))]))

# cdl$sitelat <- factor(cdl$sitelat, levels = cdl$sitelat[order(cdl$cdl, decreasing = TRUE)])
# newdat$sitelat <- factor(newdat$sitelat, levels = cdl$sitelat[order(cdl$cdl, decreasing = TRUE)])
# dat$sitelat <- factor(dat$sitelat, levels = cdl$sitelat[order(cdl$cdl, decreasing = TRUE)])


lats$zday <- (lats$corvhour - attr(dat23$zday, "scaled:center")) / attr(dat23$zday, "scaled:scale")
lats$pred <- predict(mod2, lats, re.form = ~(1 |population), type = "response")


# plot shows raw data for all chambers and model prediction for 21C chambers by population
plt <- ggplot(dat23, aes(x = photoperiod, y = perc_repro)) +
  geom_point(alpha = .3, size = 1.6) +
  scale_x_continuous(limits = c(min(newdat$photoperiod), max(newdat$photoperiod))) +
  # scale_shape_discrete(name = "Year") +
  # geom_line(data = newdat, aes(x = photoperiod, y = pred), size = 1.5, alpha = .5) +
  facet_wrap(~population, ncol = 2) +
  # geom_text(data = cdl, aes(x = cdl - .3, y = .95, label = round(cdl,1)), size = 5) +
  # geom_vline(data = cdl, aes(xintercept = cdl)) +
  # geom_segment(data = lats, aes(x = maxhour, xend = corvhour, y = .05, yend = .05), size = 2,
               # arrow = arrow(length = unit(0.2,"cm"))) +
  # geom_vline(data = lats, aes(xintercept = corvhour), linetype = "dashed", size = 1) +
  # geom_hline(data = lats, aes(yintercept = pred), linetype = "dashed", size = 1) +
  # geom_text(data = lats, aes(x = 13.75, y = pred + 0.1, label = round(pred,1)), size = 5) +
  xlab("Hours of light in chamber (23C)") +
  ylab("Percent reproductive") +
  # ggtitle("Critical photoperiod of Galerucella") +
  theme(legend.margin = margin(1, 1, 1, 1, unit = "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plt



# Each site separate GLM ----
library(broom)
# Run all models at once for site x year
mods <- dat %>% 
  group_by(population) %>% 
  do(fit = glm(perc_repro ~ photoperiod, weights = total, data = ., family = binomial(link = "probit")))

modcoefs <- bind_rows(lapply(mods$fit, FUN = function(x) data.frame(b0 = coef(x)[1], b1 = coef(x)[2]))) %>% 
  mutate(mods$population) %>% 
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
write.csv(modcoefs, file = "GCA_CP_2019.csv", row.names = FALSE)

