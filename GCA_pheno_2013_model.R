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
         Count = ifelse(Stage == "adult", round(MeanCount * 15), round(MeanCount * 150)))
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


stagedat <- dat %>% filter(Stage == "adult")

# use loess like above plot
# from stackoverflow to mimic ggplot smooth
model <- loess(RelProp ~ GDD, data = stagedat)
xrange <- range(stagedat$GDD)
xseq <- seq(from=xrange[1], to=xrange[2], length=1000)
pred <- predict(model, newdata = data.frame(GDD = xseq), se=TRUE)
y = pred$fit
ci <- pred$se.fit * qt(0.95 / 2 + .5, pred$df)
ymin = y - ci
ymax = y + ci
loess.DF <- data.frame(x = xseq, y, ymin, ymax, se = pred$se.fit)

loess.DF$y[which(loess.DF$y < 0)] <- 0
loess.DF$CDF <- cumsum(loess.DF$y) / sum(loess.DF$y, na.rm = TRUE)


ggplot(loess.DF, aes(x = xseq, y = CDF)) + geom_line() #+ geom_hline(yintercept = .95)


ggplot(stagedat, aes(x = GDD, y = RelProp)) + 
  geom_point() +
  geom_smooth(aes_auto(loess.DF), data=loess.DF, stat="identity")
  # stat_smooth(method="glm", method.args = list(family="binomial"))

# find peak of adults
loess.DF[which(loess.DF$y == max(loess.DF$y)),]
# find first adults
# loess.DF <- loess.DF[which(loess.DF$x < 200), ]
# loess.DF[which(abs(loess.DF$y) == min(abs(loess.DF$y))),]

# using an empirical distribution, split into substage means and weights
# distribution needs x and cdf
gamdist <- preddat %>% 
  mutate(CDF = cumsum(y) / sum(y, na.rm = TRUE),
         x = GDD)

SubstageDistrib <- function(dist, numstage, perc){
  # approximation from loess predictions
  ReturnClosestValue <- function(dist, xval){
    out <- dist$CDF[which(dist$x > xval)[1]]
  }
  
  low <- (1 - perc)/2
  high <- 1 - (1 - perc) / 2
  low <- dist$x[which(dist$CDF > low)[1]]
  high <- dist$x[which(dist$CDF > high)[1]]
  
  bounds <- seq(low, high, length.out = numstage + 1)
  means <- (bounds[1:numstage] + bounds[2:(numstage + 1)]) / 2
  weights <- diff(sapply(X = bounds, FUN = ReturnClosestValue, dist = dist), lag = 1)
  return(data.frame(means, weights))
}





# probably not necessary to model this with a GAM
library(mgcv)
stagedat <- dat %>% filter(Stage == "egg")



mod <- gam(Count ~ s(GDD, k = 8) + s(Site, bs = "re"), method = "REML",
           family = poisson, data = stagedat)


xrange <- range(stagedat$GDD)
xseq <- seq(from=xrange[1], to=xrange[2], length=1000)
preddat <- data.frame(GDD = xseq, Site = "Sutherlin")
pred <- predict(mod, newdata = preddat, exclude = "s(Site)", type = "response")
pred <- predict(mod, newdata = preddat, type = "response")
preddat$y = pred
# preddat$y = exp(pred$fit)

preddat$RelProp <- preddat$y / max(preddat$y, na.rm = TRUE)

plot(preddat$GDD, preddat$y)  



plt <- ggplot(stagedat, aes(x = GDD, y = RelProp, color = Site)) + 
  geom_point() +
  geom_path(data = preddat) +
  # facet_wrap(~Stage, ncol = 1) +
  theme_bw()
plt
  


