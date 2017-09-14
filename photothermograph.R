# make photothermograph
library(tidyverse)

gdd <- readRDS("../multivoltine/data/gdd_covs.rds")

sitegdd <- gdd %>%
  group_by(SiteID, year) %>%
  summarise(maxGDD = max(cumdegday)) %>%
  ungroup() %>%
  arrange(maxGDD)

# plots with daylength at different gdd
# compares norhern/southernmost sites for gdd and photoperiod
pltdat <- gdd %>%
  filter(SiteID %in% c("016", "138"))

p <- ggplot(pltdat, aes(x = cumdegday, y = photo, group = year, color = year))
p + geom_line() + facet_wrap(~ SiteID, ncol = 1) + theme_bw()


