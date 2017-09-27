# GCA phenology from observations in 2013 in OR
library(lubridate)
library(ggplot2)
library(dplyr)

dat <- read.csv("GCA_pheno_2013.csv", header = TRUE)
dat$Date <- mdy(dat$Date)
dat <- dat %>% 
  group_by(Site, Stage) %>% 
  mutate(RelProp = MeanCount / sum(MeanCount, na.rm = TRUE),
         CDF = cumsum(MeanCount) / sum(MeanCount, na.rm = TRUE),
         Count = round(MeanCount * 15))
aggdat <- dat %>% 
  group_by(Stage, Date) %>% 
  summarise(MeanRelProp = mean(RelProp))

plt <- ggplot(dat, aes(x = GDD, y = CDF, group = Stage, color = Stage)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~Site, ncol = 1) +
  theme_bw()
plt

# this seems like good smoothed data to use
plt <- ggplot(dat, aes(x = GDD, y = RelProp, group = Stage, color = Site)) + 
  geom_point() +
  geom_smooth() +
  facet_wrap(~Stage, ncol = 1) +
  theme_bw()
plt

# use loess like above plot
# from stackoverflow to mimic ggplot smooth
model <- loess(wt ~ hp, data=mtcars)
xrange <- range(mtcars$hp)
xseq <- seq(from=xrange[1], to=xrange[2], length=80)
pred <- predict(model, newdata = data.frame(hp = xseq), se=TRUE)
y = pred$fit
ci <- pred$se.fit * qt(0.95 / 2 + .5, pred$df)
ymin = y - ci
ymax = y + ci
loess.DF <- data.frame(x = xseq, y, ymin, ymax, se = pred$se.fit)




adults <- dat %>% filter(Stage == "adult")

# probably not necessary to model this
library(mgcv)
mod <- gam(Count ~ s(GDD, by = Stage, k = 5),
           family = poisson, data = dat)

  
# what's the best pre-OP time?
dat <- dat %>% 
  


