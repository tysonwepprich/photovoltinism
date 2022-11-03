# Plots from optimal CP model
theme_set(theme_bw(base_size = 16) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())) 


outlist <- readRDS("ibm_outlist_DCA_12sites_daymet.rds")
# Photothermographs/consequences plots ----
# Plot voltinism photothermograph with ridge plots
inds <- outlist %>% 
  flatten() %>%
  purrr::map(1) %>% dplyr::bind_rows()


inds <- inds %>% 
  filter(Year %in% c(2011, 2013, 2015), 
         cdl_mu == 15,
         lambda == 1.5)
gdd_yr <- gdd_all %>% 
  filter(year %in% c(2011, 2013, 2015))


# Marginal histograms ----
mean_cdl <- mean(inds$cdl[inds$cdl > 0], na.rm = TRUE)
maxgdd <- max(gdd_yr$accumdegday, na.rm = TRUE)
maxphoto <- max(c(max(gdd_yr$daylength), max(inds$cdl, na.rm = TRUE)))
minphoto <- min(gdd_yr$daylength)

cdl_hist <- conseq %>% # conseq from below
  group_by(SiteID) %>% 
  filter(cdl > 0) %>% 
  mutate(breaks = cut(cdl, breaks=seq(min(cdl),max(cdl),length.out = 30), 
                      labels=seq(min(cdl),max(cdl), length.out = 30)[-1], 
                      include.lowest=TRUE)) %>% 
  mutate(daylength = as.numeric(as.character(breaks))) %>%
  group_by(SiteID, daylength) %>% 
  summarise(n = n()) %>% 
  mutate(accumdegday = -50 + 50 * n/max(n),
         Site = SiteID)

# oa <- inds %>% filter(OA_DOY > 0) %>% 
#   left_join(gdd_all, by = c("Year" = "year", "SiteID" = "Site", "OA_DOY" = "yday"))
# teneral <- inds %>% filter(TA_DOY > 0) %>% 
#   left_join(gdd_all, by = c("Year" = "year", "SiteID" = "Site", "TA_DOY" = "yday"))
# sens_stage <- bind_rows(oa, teneral)
sens_hist <- conseq %>% # conseq from below
  group_by(SiteID) %>% 
  mutate(breaks = cut(accumdegday, breaks = seq(min(accumdegday), max(accumdegday), length.out = 150),
                      labels = seq(min(accumdegday), max(accumdegday), length.out = 150)[-1],
                      include.lowest = TRUE)) %>% 
  mutate(accumdegday = as.numeric(as.character(breaks))) %>% 
  group_by(SiteID, accumdegday) %>% 
  summarise(n = n()) %>% 
  mutate(daylength = minphoto - .5 + .5 * n/max(n),
         Site = SiteID)
sens_lines <- sens_hist %>% filter(accumdegday > 400)

gdd_yr$Site <- forcats::fct_relevel(gdd_yr$Site, levels(gdd_yr$Site)[c(1, 7, 2, 5, 6, 3, 4)])

plt <- ggplot(gdd_yr %>% filter(Site %in% inc_sites), 
              aes(x = accumdegday, y = daylength, group = as.factor(year), color = as.factor(year))) +
  geom_hline(data = cdl_hist, aes(yintercept = daylength, alpha = n/100000000), color = "light grey", size = 1, show.legend = FALSE) +
  geom_vline(data = sens_lines, aes(xintercept = accumdegday, alpha = n/100000000), color = "light grey", size = 1, show.legend = FALSE) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(-50, maxgdd + 100), ylim = c(minphoto - .5, 17.2), expand = FALSE) +
  # scale_color_viridis(discrete = TRUE, name = "Year", option = "D") +
  scale_color_brewer(name = "Year", type = "qual") + 
  scale_y_continuous(breaks = c(9:17)) +
  guides(color = guide_legend(reverse=FALSE)) +
  geom_segment(data=cdl_hist, size=1.5, show.legend=FALSE,
               aes(x=-50, xend=accumdegday, y=daylength, yend=daylength), inherit.aes = FALSE) +
  geom_segment(data=sens_hist, size=2, show.legend=FALSE,
               aes(x=accumdegday, xend=accumdegday, y=minphoto - .5, yend=daylength), inherit.aes = FALSE) +
  # geom_rug(data = results, aes(x = accumdegday), color = "gray", alpha = .1) +
  # geom_rug(data = results, aes(y = cdl), color = "gray", alpha = .1) +
  # ggtitle("Photothermographs with critical photoperiod\nand simulated sensitive stage emergence") +
  xlab("Accumulated degree-days") +
  ylab("Daylength (hours)") +
  facet_wrap(~Site, nrow = 3) +
  theme(legend.position = c(.85, .8))
plt

ggsave(filename = "GCA_example_PT.png", plot = plt, device = "png", width = 6, height = 9, units = "in")



# Consequences ridgelines (stacked bar doesn't look as good here, 
# but don't know how to combine these plots with photothermograph)

uniqsims <- distinct(firsts[,1:5])

outcomelist <- list()
for (nr in 1:nrow(uniqsims)){
  ids <- uniqsims[nr, ]
  temp <- inds %>% filter(Year == ids$Year, 
                          SiteID == ids$SiteID, 
                          cdl_mu == ids$cdl_mu,
                          cdl_sd == ids$cdl_sd,
                          lambda == ids$lambda)
  gdd <- gdd_all %>% filter(year == ids$Year, 
                            Site == ids$SiteID) %>% 
    mutate(doy = yday)
  
  temp_conseq <- temp %>%
    mutate(outcome = case_when(active == 1 ~ "lost",
                               diapause == 1 ~ "diapause",
                               diapause != 1 ~ "fits"))
  
  
  out_lost <- temp_conseq %>% 
    filter(indiv %in% parent[outcome == "lost"])
  oa <- out_lost %>% filter(OA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "OA_DOY" = "yday"))
  teneral <- out_lost %>% filter(TA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "TA_DOY" = "yday"))
  out_lost <- bind_rows(oa, teneral) %>% 
    mutate(offspring = "lost")
  
  out_diap <- temp_conseq %>% 
    filter(indiv %in% parent[outcome == "diapause"])
  oa <- out_diap %>% filter(OA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "OA_DOY" = "yday"))
  teneral <- out_diap %>% filter(TA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "TA_DOY" = "yday"))
  out_diap <- bind_rows(oa, teneral) %>% 
    mutate(offspring = "diapause")
  
  out_fits <- temp_conseq %>% 
    filter(indiv %in% parent[outcome == "fits"])
  oa <- out_fits %>% filter(OA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "OA_DOY" = "yday"))
  teneral <- out_fits %>% filter(TA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "TA_DOY" = "yday"))
  out_fits <- bind_rows(oa, teneral) %>% 
    mutate(offspring = "fits")
  
  outcomes <- bind_rows(out_lost, out_diap, out_fits)
  
  no_off <- anti_join(temp_conseq, outcomes, by = "indiv") %>% 
    filter(outcome == "diapause") %>% 
    filter(TA_DOY > 0) %>% 
    left_join(gdd, by = c("Year" = "year", "TA_DOY" = "yday")) %>% 
    mutate(offspring = "none, diapause")
  
  outcomes <- bind_rows(outcomes, no_off)
  outcomelist[[length(outcomelist) + 1]] <- outcomes
}

outcomes <- bind_rows(outcomelist)
conseq <- outcomes %>%
  mutate(Consequence = case_when(outcome == "diapause" ~ "Diapauses",
                                 offspring == "lost" ~ "Reproduces, offspring lost",
                                 offspring != "lost" ~ "Reproduces, offspring fits"),
         date = as.Date(doy, origin=as.Date("2015-12-31"))) %>% 
  filter(cdl_mu == 15, 
         lambda == 1.5,
         Year %in% c(2011, 2013, 2015))

conseq$Consequence <- factor(conseq$Consequence, levels = c("Diapauses", "Reproduces, offspring lost", "Reproduces, offspring fits"))
conseq$Year <- factor(as.character(conseq$Year), levels = c("2015", "2013", "2011"))
conseq$SiteID <- forcats::fct_relevel(conseq$SiteID, levels(conseq$SiteID)[c(1, 7, 2, 5, 6, 3, 4)])

library(ggridges)
library(viridis)

plt2 <- ggplot(conseq, aes(x = accumdegday, y = Year, 
                           fill = Consequence, height = ..count..)) +
  # geom_histogram() +
  geom_density_ridges(scale = .95, rel_min_height=0.001, stat = "binline", bins = 100) +
  scale_fill_viridis(discrete = TRUE, alpha = .25, drop = FALSE) +
  # scale_x_continuous(limits = c(-50, maxgdd + 100), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0.1, .1)) +
  ylab("Year") + xlab("Accumulated degree-days (base 10C)") +
  facet_wrap(~SiteID, nrow = 3) +
  theme(legend.position = c(.7, .8))
plt2


ggsave(filename = "GCA_example_conseq.png", plot = plt2, device = "png", width = 6, height = 9, units = "in")





# TODO
# 1. Adult "counts" for each generation
# 2. Egg "counts" for each generation as alternative monitoring stages
# 3. Proportion diapause in the field from each generation ("typical" penultimate most important)
# 4. Adult counts to voltinism based on relative size to "typical" penultimate
# 5. Lineage analysis too much, what about comparing parents in penultimate generations and their CP/emergence time


# REDO added columns since last time
census <- firsts %>% 
  group_by(SiteID, Year, cdl_mu, cdl_sd, lambda) %>% 
  pivot_longer(
    cols = one_of(names(firsts)[c(9:14, 16:26)]),
    names_to = c("choice", "numgen"),
    names_sep = "_",
    values_to = "count",
    values_drop_na = FALSE
  )

# voltinism by proportion diapause
# assumption: count the first generation where >.5 diapause
# problem: voltinism decimal never above x.5
census$count[is.na(census$count)] <- 0

diapvolt <- census %>% 
  # filter(cdl_mu %in% c(14, 15, 16)) %>% 
  group_by(SiteID, Year, cdl_mu, cdl_sd, lambda, numgen) %>% 
  summarise(propdiap = count[choice == "diap"] / count[choice == "att"],
            n = count[choice == "att"],
            voltdiap = as.numeric(numgen[1]) + (1 - propdiap),
            lost = lost[1]) %>% 
  ungroup() %>% 
  filter(complete.cases(.)) %>% 
  group_by(SiteID, Year, cdl_mu, cdl_sd, lambda) %>%  
  mutate(maxgen = max(numgen),
         keeprow = case_when(numgen == 1 ~ "yes",
                             propdiap > .95 ~ "no",
                             n < 200 ~ "no",
                             TRUE ~ "yes")) %>% 
  filter(keeprow != "no") %>% 
  arrange(-n) %>% 
  slice(1)



counts <- seconds %>% 
  group_by(SiteID, Year, cdl_mu, cdl_sd, lambda, numgen, lifestage) %>% 
  summarise(peakn = max(n), 
            peakdoy = doy[n == max(n)][1])

# how to remove weird years? what is typical voltinism
relvolt <- counts %>% 
  # filter(cdl_mu %in% c(14, 15, 16)) %>% 
  filter(lifestage == 5) %>%
  filter(lambda > 1) %>% 
  # filter(SiteID == "Vantage, WA", Year == 1994, cdl_mu <= 14, lifestage == 5) %>% 
  group_by(SiteID, Year, cdl_mu, cdl_sd, lambda) %>% 
  arrange(SiteID, Year, cdl_mu, cdl_sd, lambda, numgen) %>% 
  mutate(dpeak = c(0, diff(peakn))) %>% 
  summarise(maxgen = max(numgen),
            lastn = case_when(maxgen == 1 ~ peakn[1],
                              maxgen == 2 ~ peakn[2],
                              min(dpeak) < 0 ~ peakn[which(dpeak < 0)[1]],
                              TRUE ~ peakn[maxgen]),
            penultn = case_when(maxgen == 1 ~ peakn[1]-peakn[1],
                                maxgen == 2 ~ peakn[1],
                                maxgen > 2 & min(dpeak) < 0 ~ peakn[which(dpeak < 0)[1] - 1L],
                                # maxgen > 2 & min(dpeak) >= 0 ~ peakn[maxgen - 1],
                                TRUE ~ NA_integer_),
            penultgen = case_when(maxgen == 1 ~ as.integer(0L),
                                  maxgen == 2 ~ as.integer(1L),
                                  min(dpeak) < 0 ~ as.integer(numgen[which(dpeak < 0)[1]-1]),
                                  TRUE ~ as.integer(maxgen - 1L)),
            lastprop = lastn / (lastn + penultn),
            voltcount = penultgen + lastprop)


dat <- firsts %>% 
  mutate(voltmodel = attvoltinism) %>% 
  dplyr::select(SiteID, Year, cdl_mu, lambda, ann_lam, voltmodel) %>% 
  right_join(relvolt) %>% 
  dplyr::select(-maxgen) %>% 
  right_join(diapvolt) %>% 
  mutate(lost = log(lost + 1),
         ann_lam = log(ann_lam + 0.1))

pltdat <- dat %>% 
  group_by(SiteID, Year, cdl_mu, lambda) %>% 
  pivot_longer(
    cols = one_of(c("ann_lam", "voltmodel", "voltcount", "voltdiap", "lost")),
    names_to = "var",
    values_to = "value",
    values_drop_na = FALSE
  ) #%>% 
# filter(lambda == 1.5)

pltdat$var <- as.factor(pltdat$var)
pltdat$var <- plyr::revalue(pltdat$var, c("ann_lam"="Log(Annual lambda)", 
                                          "voltmodel"="Modeled voltinism",
                                          "voltcount"="Counted voltinism",
                                          "voltdiap"="Diapause voltinism",
                                          "lost"="Log(Lost)"))
pltdat$var <- factor(pltdat$var, levels = c("Modeled voltinism", "Diapause voltinism", 
                                            "Counted voltinism", "Log(Lost)", "Log(Annual lambda)"))

ggplot(pltdat, aes(x = Year, y = value, group = cdl_mu, color = as.factor(cdl_mu))) +
  geom_line() +
  scale_color_brewer(type = "qual", palette = "Dark2", name = "Critical\nphotoperiod") +
  facet_grid(var~SiteID, scales = "free") +
  ylab("Annual value") +
  theme_bw()

cor.test(pltdat$value[which(pltdat$var == "Log(Annual lambda)")],
         pltdat$value[which(pltdat$var == "Modeled voltinism")], method = "kendall")
cor.test(pltdat$value[which(pltdat$var == "Log(Annual lambda)")],
         pltdat$value[which(pltdat$var == "Counted voltinism")], method = "kendall")
cor.test(pltdat$value[which(pltdat$var == "Log(Annual lambda)")],
         pltdat$value[which(pltdat$var == "Diapause voltinism")], method = "kendall")

cordat <- pltdat %>% 
  pivot_wider(names_from = var, values_from = value)

ggplot(cordat, aes(x = `Diapause voltinism`, y = `Log(Annual lambda)`, color = `Log(Lost)`)) +
  geom_point() +
  theme_bw()

ggplot(relvolt, aes(x = peakn, group = numgen, color = numgen)) +
  geom_density() +
  facet_grid(SiteID ~ cdl_mu, scales = "free")

ggplot(firsts %>% filter(lambda == 2, cdl_mu %in% c(14, 15, 16)), aes(x = Year, y = attvoltinism, group = cdl_mu, color = as.factor(cdl_mu))) +
  geom_line() +
  facet_wrap(~SiteID) +
  theme_bw()
ggplot(firsts %>% filter(lambda == 2, cdl_mu %in% c(14, 15, 16)), aes(x = Year, y = log(ann_lam), group = cdl_mu, color = as.factor(cdl_mu))) +
  geom_line() +
  facet_wrap(~SiteID) +
  theme_bw()
ggplot(firsts %>% filter(cdl_mu == 15, lambda %in% c(1, 2, 3)), aes(x = Year, y = attvoltinism, group = lambda, color = as.factor(lambda))) +
  geom_line() +
  facet_wrap(~SiteID) +
  theme_bw()
ggplot(firsts %>% filter(cdl_mu == 15, lambda %in% c(1, 2, 3)), aes(x = Year, y = log(ann_lam), group = lambda, color = as.factor(lambda))) +
  geom_line() +
  facet_wrap(~SiteID) +
  theme_bw()

library(ggpubr) # for combining plots

pltdat2 <- seconds %>% 
  # filter(cdl_mu == 13.5, lambda == 1.5, SiteID == "Sutherlin, OR", Year == 2015, lifestage == 5)
  filter(cdl_mu == 14, lambda == 1.5, SiteID == "Palermo, CA", Year == 2015, lifestage == 5)
# filter(cdl_mu == 16, lambda == 1, SiteID == "Bellingham, WA", Year == 2014, lifestage == 5)
# filter(cdl_mu == 16, lambda == 2, SiteID == "Bellingham, WA", Year == 2014, lifestage == 5)
# filter(cdl_mu == 15.5, lambda == 1, SiteID == "Bellingham, WA", Year == 2014, lifestage == 5)
# filter(cdl_mu == 16, lambda == 1, SiteID == "Bellingham, WA", Year == 2015, lifestage == 5)

pltdat2$diapause <- as.factor(pltdat2$diapause)
pltdat2$diapause <- plyr::revalue(pltdat2$diapause, c("-1"="Emergence", 
                                                      "0"="Reproduction",
                                                      "1"="Diapause"))
pltdat2$diapause <- factor(pltdat2$diapause, levels = c("Diapause", "Reproduction", "Emergence"))

b <- ggplot(pltdat2,
            aes(x = doy, y = n, fill = diapause)) +
  scale_fill_brewer(type = "qual", palette = "Dark2", name = "Adult decision") +
  # scale_x_continuous(limits = c(170, 350)) +
  # scale_y_continuous(limits = c(0, 1000)) +
  geom_histogram(stat = "identity", width = 1) +
  theme_bw(base_size = 16) +
  theme(legend.position = c(.8, .8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Count") +
  xlab("Day of year")
b

e <- ggarrange(a,b,c,d, 
               labels = c("A", "B", "C", "D"),
               ncol = 2, nrow = 2, align = "v", common.legend = TRUE)
e

e <- ggarrange(a,b, 
               labels = c("A", "B"),
               ncol = 2, nrow = 1, align = "h", common.legend = TRUE)
e

pltdat %>% 
  # filter(cdl_mu == 14, lambda == 1.5, SiteID == "Palermo, CA", Year == 2015) %>% 
  filter(cdl_mu == 13.5, lambda == 1.5, SiteID == "Sutherlin, OR", Year == 2015) %>%
  # filter(cdl_mu == 16, lambda == 1, SiteID == "Bellingham, WA", Year == 2014) %>% 
  # filter(cdl_mu == 16, lambda == 2, SiteID == "Bellingham, WA", Year == 2014) %>%
  # filter(cdl_mu == 15.5, lambda == 1, SiteID == "Bellingham, WA", Year == 2014) %>%
  # filter(cdl_mu == 16, lambda == 1, SiteID == "Bellingham, WA", Year == 2015) %>%
  data.frame()


# Optimal IBM plots ----
res <- bind_rows(outlist)
# res <- bind_rows(flatten(outlist)) %>% 
# group_by(Year, SiteID, cdl_mu, cdl_sd, lambda) %>% 
# summarise_at(vars(voltinism:diap_2), mean, na.rm = TRUE)
saveRDS(res, "ibm_results_GCA_Yakima_daymet.rds")

# Quanitfying voltinism
library(viridis)

# CP vs voltinism with log(lambda) over years for a single site
plt <- ggplot(firsts %>% dplyr::filter(lambda == 1.5, cdl_sd == .5, SiteID == "Delta, UT"), aes(x = cdl_mu, y = attvoltinism, color = log(ann_lam_adj_5))) +
  geom_point() +
  scale_color_viridis(name = "Log(lambda)") +
  facet_wrap(~Year, ncol = 5) +
  theme_bw(base_size = 14) +
  xlab("Simulated critical photoperiod") +
  ylab("Mean number attempted") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

plt

res <- readRDS("ibm_results_dca.rds")


# Optimal CDL plots

###############
res <- firsts %>% 
  filter(cdl_sd == .25, lambda == 1.5, SiteID == "Delta, UT", Year > 2009)
# filter(cdl_sd == .33, lambda == 1.5, SiteID %in% unique(SiteID)[c(1, 2, 6, 7, 9, 10)])

plt <- ggplot(res, aes(x = cdl_mu, y = attvoltinism)) +
  # geom_point(aes(color = ann_lam), size = 3.5) +
  geom_point(aes(color = log(ann_lam + 0.01))) +
  scale_x_reverse() +
  scale_color_viridis(begin = 0, end = 1) +
  facet_wrap(~Year, ncol = 3) +
  # theme_bw() +
  ylab("Voltinism: weighted average of generations attempted") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated voltinism and projected growth rate by critical photoperiod")
plt

ggsave(filename = paste0("GCA_voltinism_", site, ".png"), plot = plt, device = "png", width = 9, height = 6, units = "in")


# heatmap
plt <- ggplot(res, aes(x = cdl_mu, y = lambda, fill = log(ann_lam))) +
  geom_raster() +
  scale_x_reverse() +
  scale_fill_viridis(begin = 1, end = 0) +
  facet_wrap(~Year, ncol = 3) +
  theme_bw() +
  ylab("Between generation lambda") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated log(Annual Population Growth Rate)", subtitle = "Accounting for lost generations")
plt


# lines
plt <- ggplot(res, aes(x = cdl_mu, y = log(ann_lam), group = lambda, color = lambda)) +
  geom_line() +
  # scale_x_reverse() +
  scale_color_viridis(begin = 1, end = 0) +
  facet_wrap(~Year, ncol = 3) +
  theme_bw() +
  ylab("Log(annual lambda)") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated log(Annual Population Growth Rate)", subtitle = "Accounting for lost generations")
plt

# across years
res <- readRDS("ibm_results_west.rds") # %>%  filter(SiteID == 123)
res$SiteID <- factor(res$SiteID, levels = c("Lovell, WY",
                                        "Delta, UT", "St. George, UT", "Princess, NV",
                                        "Big Bend, NV", "Blythe, CA", "Cibola, CA",
                                        "Imperial, AZ"))
res$SiteID <- plyr::revalue(res$SiteID, c("Yakima Training Center, WA"="Yakima TC, WA", 
                                          "Baskett Slough, OR"="Rickreall, OR"))
res$SiteID <- factor(res$SiteID, levels = c("Bellingham, WA", "Yakima TC, WA", "Rickreall, OR",
                                            "Sutherlin, OR", "McArthur, CA", "Palermo, CA"))

cps <- read.csv("DCA_CP.csv", header = TRUE) %>%
  filter(Site %in% res$SiteID, Year == 2019) %>%
  rename(SiteID = Site)
cps <- read.csv("GCA_CP.csv", header = TRUE) %>% 
  filter(year == 2019) %>% 
  rename(SiteID = population, CP = cp)

res <- res %>% filter(SiteID %in% cps$SiteID) %>% droplevels()

firsts <- firsts %>% filter(SiteID %in% cps$SiteID) %>% droplevels()

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geom_lam <- firsts %>% 
  # filter(Year > 2007) %>% 
  group_by(SiteID, cdl_mu, cdl_sd, lambda) %>% 
  # tidyr::complete(Year, cdl_mu, cdl_sd, lambda, fill = list(ann_lam = 0)) %>% 
  summarise(mean_annual_lambda_0 = gm_mean(ann_lam + 0.001),
            mean_annual_lambda_1 = gm_mean(ann_lam_adj_5 + 0.001))
# mean_annual_lambda_2 = gm_mean(ann_lam_adj_2500 + 0.001))

geom_lam$SiteID <- factor(geom_lam$SiteID, levels = c("Lovell, WY",
                                            "Delta, UT", "St. George, UT", "Princess, NV",
                                            "Big Bend, NV", "Blythe, CA", "Cibola, CA",
                                            "Imperial, AZ"))

plt <- ggplot(geom_lam, aes(x = cdl_mu, y = lambda, fill = log(mean_annual_lambda_1))) +
  geom_raster() +
  # scale_x_reverse() +
  scale_fill_viridis(begin = 0, end = 1) +
  theme_bw() +
  facet_grid(~SiteID) +
  ylab("Between generation lambda") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated critical photoperiod and mean growth rates", subtitle = "1994-2018 geometric mean log(lambda)")
plt

# lines
plt <- ggplot(geom_lam %>% filter(lambda == 2, cdl_sd == 1), 
              aes(x = cdl_mu, y = log(mean_annual_lambda_1))) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7) +
  # scale_color_continuous(name = "CP std. dev.", breaks = unique(geom_lam$cdl_sd), labels = unique(geom_lam$cdl_sd)) +
  geom_vline(data = cps, aes(xintercept = CP), linetype = "dotted", alpha = 0.7) +
  scale_x_continuous(breaks = c(8, 10, 12, 14, 16, 18)) +
  # scale_color_viridis(begin = 1, end = 0) +
  xlab("Modeled critical daylength") +
  ylab("Modeled log mean annual growth rate") +
  facet_wrap(~SiteID)
plt

# Geometric mean vs annual
# lines
res$SiteID <- factor(res$SiteID, levels = levels(res$SiteID)[c(1,6,4,5,2,3)])


plt <- ggplot(res %>% filter(lambda == 2, cdl_sd == .25), aes(x = cdl_mu, y = log(ann_lam), group = Year, color = Year)) +
  geom_line(alpha = .5) +
  # scale_color_viridis(begin = .9, end = 0) +
  geom_line(data = geom_lam %>% filter(lambda == 2, cdl_sd == .25), aes(x = cdl_mu, y = log(mean_arith_lambda)), size = 1, inherit.aes = FALSE) +
  facet_wrap(~SiteID, ncol = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Log(annual lambda)") +
  xlab("Mean critical photoperiod")
# theme(legend.position = c(.85, .15)) +
ggtitle("Simulated mean log(annual growth rate) 1993-2018 at constant critical photoperiod distributions", subtitle = "Vertical lines are critical photoperiod (mean) measured in 2019 growth chamber experiments")
plt



# Optimal CP paper ----
outlist <- readRDS("ibm_outlist_GCA_6sites_daymet.rds")
sites <- data.frame(ID = c("Ephrata, WA", "Yakima TC, WA",
  "Sutherlin, OR", "Bellingham, WA",
  "McArthur, CA", "Palermo, CA", "Rickreall, OR"),
  x = c(-119.655253, -119.986530,
          -123.315854, -122.479482,
          -121.41, -121.58, -123.27),
  y = c(47.161647, 46.756318,
          43.387721, 48.756105,
          41.10, 39.41, 44.98))[-1,]


firsts <- outlist %>% 
  flatten() %>%
  purrr::map(1) %>% dplyr::bind_rows()
seconds <- outlist %>% 
  flatten() %>%
  purrr::map(2) %>% dplyr::bind_rows()



firsts <- firsts %>% 
  left_join(sites, by = c("SiteID" = "ID")) %>% 
  mutate(SiteID = forcats::fct_reorder(SiteID, y, .desc = TRUE))

# best cdl combo
summ_lam <- firsts %>% 
  group_by(SiteID, cdl_sd, lambda, Year) %>% 
  mutate(maxattgen = max(attvoltinism, na.rm = TRUE),
         maxcomgen = max(voltinism, na.rm = TRUE),
         maxlambda = max(ann_lam),
         maxlambda5 = max(ann_lam_adj_5)) %>% 
  group_by(SiteID, cdl_mu, cdl_sd, lambda) %>% 
  summarise(mean_annual_lambda = gm_mean(ann_lam + 0.001),
            mean_annual_lambda_5 = gm_mean(ann_lam_adj_5 + 0.001),
            mean_annual_lambda_10 = gm_mean(ann_lam_adj_10 + 0.001),
            mean_lost = mean(lost),
            mean_wvolt = mean(ifelse(is.na(voltinism), 0, voltinism)),
            mean_attvolt = mean(attvoltinism),
            mean_mismatch = mean(attvoltinism - maxcomgen),
            mean_lambda5_mm = mean(ann_lam_adj_5 - maxlambda5),
            sd_annual_lambda = sd(ann_lam + 0.001),
            sd_lost = sd(lost),
            sd_wvolt = sd(ifelse(is.na(voltinism), 0, voltinism)),
            sd_attvolt = sd(attvoltinism),
            sd_mismatch = sd(attvoltinism - maxcomgen))


year_best_cp <- firsts %>% 
  group_by(SiteID, cdl_sd, lambda, Year) %>% 
  mutate(maxattgen = max(attvoltinism, na.rm = TRUE),
         maxcomgen = max(voltinism, na.rm = TRUE),
         maxlambda = max(ann_lam),
         maxlambda5 = max(ann_lam_adj_5),
         ann_lam_mm = ann_lam_adj_5 - maxlambda5) %>% 
  summarise(best_ann_cp = mean(cdl_mu[which(ann_lam_mm == 0)], na.rm = TRUE))
    
best_cp <- summ_lam %>% 
  group_by(SiteID, cdl_sd, lambda) %>% 
  mutate(maxlam = max(mean_annual_lambda_5),
         difflam = mean_annual_lambda_5 - maxlam) %>% 
  summarise(best_cp = mean(cdl_mu[which(difflam == 0)], na.rm = TRUE),
            mean_lost = mean(mean_lost[which(difflam == 0)], na.rm = TRUE))

# inc_sites <- c("Lovell, WY", "Delta, UT", "Imperial, AZ")
# inc_sites <- sites$ID
inc_sites <- c("Bellingham, WA", "Sutherlin, OR", "Palermo, CA")
dat <- firsts %>% filter(lambda == 2.5, cdl_sd == .25, SiteID %in% inc_sites)
bestline <- summ_lam %>% filter(lambda == 2.5, cdl_sd == .25, SiteID %in% inc_sites)
yearcp <- year_best_cp %>% filter(lambda == 2.5, cdl_sd == .25, SiteID %in% inc_sites)
bestcp <- best_cp %>% filter(lambda == 2.5, cdl_sd == .25, SiteID %in% inc_sites)

plt <- ggplot(dat, aes(x = cdl_mu, y = log(ann_lam_adj_5 + 0.001), group = Year)) +
  geom_line(alpha = .2, color = "black") +
  # scale_color_viridis(begin = .9, end = 0) +
  scale_x_continuous(breaks = c(10:20)) +
  geom_line(data = bestline, aes(x = cdl_mu, y = log(mean_annual_lambda_5)), size = 1, inherit.aes = FALSE) +
  facet_wrap(~SiteID, ncol = 1) +
  geom_rug(data = yearcp, aes(x = best_ann_cp), inherit.aes = FALSE) +
  geom_vline(data = bestcp, aes(xintercept = best_cp), linetype = "dashed") +
  coord_cartesian(ylim = c(-4, 4), xlim = c(10, 19)) +  
  ylab("Log(population growth rate)") +
  xlab("Critical photoperiod (day hours)")
# theme(legend.position = c(.85, .15))
plt

ggsave(filename = "GCA_AnnualLambda_3site.png", plot = plt, device = "png", width = 6, height = 12, units = "in")



# Does best CP vary by SD and lambda?
inc_sd <- c(0, .5, 1, 1.5)
# dat <- firsts %>% filter(lambda == 1.5, cdl_sd == .25, SiteID %in% inc_sites)
bestline <- summ_lam %>% filter(SiteID %in% inc_sites, lambda == 1.5, cdl_sd %in% inc_sd)
yearcp <- year_best_cp %>% filter(SiteID %in% inc_sites, lambda == 1.5, cdl_sd %in% inc_sd)
bestcp <- best_cp %>% filter(SiteID %in% inc_sites, lambda == 1.5, cdl_sd %in% inc_sd)



plt <- ggplot(bestline, aes(x = cdl_mu, y = log(mean_annual_lambda_5), group = cdl_sd, color = cdl_sd)) +
  geom_line(size = 1) +
  scale_color_continuous(name = "Standard deviation of\nphotoperiod response") +
  scale_x_continuous(breaks = c(10:19)) +
  facet_wrap(~SiteID, ncol = 1) +
  geom_vline(data = bestcp, aes(xintercept = best_cp, color = cdl_sd), linetype = "dashed") +
  coord_cartesian(ylim = c(-4, 2), xlim = c(10, 19)) +  
  ylab("Log(population growth rate)") +
  xlab("Critical photoperiod (day hours)") +
  theme(legend.position = "bottom")
plt

ggsave(filename = "GCA_AnnualLambda_3site_SD.png", plot = plt, device = "png", width = 6, height = 12, units = "in")


inc_sd <- .5
# dat <- firsts %>% filter(lambda == 1.5, cdl_sd == .25, SiteID %in% inc_sites)
bestline <- summ_lam %>% filter(SiteID %in% inc_sites, cdl_sd %in% inc_sd)
yearcp <- year_best_cp %>% filter(SiteID %in% inc_sites, cdl_sd %in% inc_sd)
bestcp <- best_cp %>% filter(SiteID %in% inc_sites, cdl_sd %in% inc_sd)



plt <- ggplot(bestline, aes(x = cdl_mu, y = log(mean_annual_lambda_5), group = lambda, color = lambda)) +
  geom_line(size = 1) +
  scale_color_continuous(name = "Between\ngeneration\ngrowth rate") +
  scale_x_continuous(breaks = c(10:19)) +
  facet_wrap(~SiteID, ncol = 1) +
  geom_vline(data = bestcp, aes(xintercept = best_cp, color = lambda), linetype = "dashed") +
  coord_cartesian(ylim = c(-5, 3), xlim = c(10, 19)) +  
  ylab("Log(population growth rate)") +
  xlab("Critical photoperiod (day hours)")
# theme(legend.position = c(.85, .15))
plt

ggsave(filename = "GCA_AnnualLambda_3site_lambda.png", plot = plt, device = "png", width = 6, height = 12, units = "in")


# Parameter effect on best CP ----

summ_lam <- firsts %>% 
  pivot_longer(cols = starts_with("ann_lam"), names_to = "penalty", values_to = "ann_lambda") %>% 
  group_by(SiteID, cdl_mu, cdl_sd, lambda, penalty) %>% 
  summarise(mean_annual_lambda = gm_mean(ann_lambda + 0.001))
          
best_cp <- summ_lam %>% 
  group_by(SiteID, cdl_sd, lambda, penalty) %>% 
  mutate(maxlam = max(mean_annual_lambda),
         difflam = mean_annual_lambda - maxlam) %>% 
  summarise(best_cp = mean(cdl_mu[which(difflam == 0)], na.rm = TRUE))

best_cp <- best_cp %>% left_join(gdd_seas[,c("Site", "dd_fn1")], by = c("SiteID" = "Site"))

mod <- lm(best_cp ~ (dd_fn1 + cdl_sd + lambda + as.factor(penalty))^2, data = best_cp)
mod <- lmer(best_cp ~ (dd_fn1 + cdl_sd + lambda + as.factor(penalty))^2 + (1|SiteID), data = best_cp)

anova(mod)
car::Anova(mod, type = 2)

bestmod <- step(mod, direction = "backward")


# Plot GC vs Best CP ----
optcp <- read.csv("data/DCA_CP.csv", header = TRUE) 
optcp <- optcp[-9, ] %>% 
  group_by(SiteID) %>% 
  arrange(-Year) %>% 
  slice(1)

optcp <- read.csv("data/GCA_CP_2019.csv", header = TRUE)

cpdistrib <- lapply(1:nrow(optcp), FUN = function(x){
  tmp <- data.frame(SiteID = optcp$mods.population[x],
                    cp = seq(10,20,length.out = 1000),
                    dens = dnorm(seq(10,20,length.out = 1000), mean = optcp$cp_mean[x], sd = optcp$cp_sd[x]))
  return(tmp)
})

cpdist <- bind_rows(cpdistrib) %>% 
  left_join(optcp, by = c("SiteID" = "mods.population")) %>% 
  left_join(sites, by = c("SiteID" = "ID")) %>% 
  filter(cp > cp_mean - 2*cp_sd & cp < cp_mean + 2*cp_sd) %>% 
  mutate(SiteID = forcats::fct_reorder(SiteID, y, .desc = FALSE))


library(ggridges)

plt <- ggplot(cpdist, aes(x = cp, y = SiteID, height = dens)) +
  geom_ridgeline(scale = 1, alpha = .3) + 
  scale_x_continuous(breaks = c(10, 12, 14, 16, 18, 20)) +
  geom_point(data = best_cp, #%>% filter(cdl_sd == 0, lambda == 1.5, penalty == "ann_lam_adj_5"),
             aes(x = best_cp, y = SiteID), size = 3, alpha = .1, inherit.aes = FALSE) +
  scale_y_discrete(expand = expand_scale(mult = c(0.01, .1))) +
  xlab("Photoperiod threshold (day hours)") +
  ylab(NULL)
plt

ggsave(filename = "GCA_CPcompare.png", plot = plt, device = "png", width = 6, height = 8, units = "in")



# better than lat or gdd?

dat <- optcp %>% 
  left_join(gdd_seas[,c("Site", "dd_fn1")], by = c("SiteID" = "Site")) %>% 
  left_join(best_cp %>% filter(cdl_sd == 0.5, lambda == 1.5, penalty == "ann_lam_adj_5")) %>% 
  left_join(sites, by = c("SiteID" = "ID"))


summary(lm(cp_mean ~ y, data = dat))
summary(lm(cp_mean ~ y + dd_fn1, data = dat))
summary(lm(cp_mean ~ y * dd_fn1, data = dat))
summary(lm(cp_mean ~ best_cp, data = dat))





year_lam <- res %>% 
  left_join(maxvolt, by = c("Year" = "year", "SiteID" = "Site")) %>% 
  group_by(Year, SiteID, cdl_mu, cdl_sd, lambda) %>% 
  mutate(mean_mismatch = mean(attvoltinism - maxgen))
best_yr <- year_lam %>% 
  filter(lambda == 2, cdl_sd == .25) %>%
  group_by(SiteID, Year) %>%
  mutate(diffbest = ann_lam - max(ann_lam)) %>% 
  filter(diffbest == 0) %>% 
  arrange(SiteID, -cdl_mu) %>% 
  slice(1) %>% 
  mutate(type = "yr")
opt_site <- bestcdltab %>% 
  ungroup() %>% 
  dplyr::select(SiteID, cdl_mu) %>% 
  rename(bestcdl = cdl_mu)
opt_yr <- res %>% 
  left_join(opt_site) %>% 
  filter(cdl_mu == bestcdl) %>% 
  filter(lambda == 2, cdl_sd == .25) %>% 
  mutate(type = "avg",
         opt_volt = voltinism, 
         opt_ann_lam = ann_lam) %>% 
  dplyr::select(SiteID, Year, bestcdl, opt_volt, opt_ann_lam)

# mismatch <- bind_rows(best_yr, opt_yr) %>% 
#   group_by(SiteID, Year) %>% 
#   arrange(SiteID, Year, type) %>% 
#   summarise(cdl_diff = diff(cdl_mu),
#             lam_diff = diff(ann_lam),
#             volt_diff = diff(voltinism))

meanmm <- mismatch %>% 
  group_by(SiteID) %>% 
  summarise_at(vars(cdl_diff:volt_diff), mean)

mismatch <- left_join(best_yr, opt_yr) %>% 
  mutate(voltmm = opt_volt - voltinism,
         cdlmm = bestcdl - cdl_mu )
mismatch$SiteID <- factor(mismatch$SiteID, levels = levels(res$SiteID))
best_yr$SiteID <- factor(best_yr$SiteID, levels = levels(res$SiteID))

# comparing yearly best choice versus geom mean best choice
plt <- ggplot(mismatch, aes(x = log(ann_lam), y = log(opt_ann_lam), color = cdlmm)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~SiteID, ncol = 2)
plt


plt <- ggplot(best_yr, aes(x = cdl_mu)) +
  geom_histogram(bins = 40, fill = "grey") +
  facet_wrap(~SiteID, ncol = 2) +
  geom_vline(data = opt_site, aes(xintercept = bestcdl), linetype = "dashed") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Mean critical photoperiod")

plt




# maximize geometric mean of annual lambda (non including overwintering)
best_cdl <- geom_lam %>% 
  # filter(lambda == 1.5) %>%
  group_by(SiteID, lambda, cdl_sd) %>%
  # mutate(diffbest = mean_annual_lambda - max(mean_annual_lambda)) %>% 
  mutate(diffbest = mean_arith_lambda - max(mean_arith_lambda)) %>% 
  filter(diffbest == 0) %>% 
  arrange(SiteID, cdl_mu) %>% 
  slice(1) #%>%

best_cdl$SiteID <- factor(best_cdl$SiteID, levels = levels(best_cdl$SiteID)[c(1,6,4,5,2,3)])


plt <- ggplot(best_cdl, aes(x = lambda, y = cdl_mu, group = as.factor(cdl_sd), color = as.factor(cdl_sd))) +
  geom_point(  position = position_dodge(width = .2)) +
  scale_color_viridis(name = "Std. dev. CP", discrete = TRUE, begin = .2, end = .8) +
  facet_wrap(~SiteID, ncol = 2) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Optimal mean critical photoperiod") +
  xlab("Between-generation population growth rate")
plt

bestcdltab <- best_cdl %>% 
  filter(cdl_sd == .25, lambda == 2) %>% 
  arrange(SiteID)
write.csv(round_df(bestcdltab), "bestcdl.csv", row.names = FALSE)


best_wvolt <- geom_lam %>% 
  filter(lambda == 1.5) %>% 
  group_by(SiteID) %>% 
  mutate(diffbest = mean_wvolt - max(mean_wvolt)) %>% 
  filter(diffbest == 0) %>% 
  arrange(SiteID, cdl_mu) %>% 
  slice(1)

best_avolt <- geom_lam %>% 
  filter(lambda == 1.5) %>% 
  group_by(SiteID) %>% 
  mutate(diffbest = mean_attvolt - max(mean_attvolt)) %>% 
  filter(diffbest == 0) %>% 
  arrange(SiteID, cdl_mu) %>% 
  slice(1)

best_lost <- geom_lam %>% 
  filter(lambda == 1.5) %>% 
  group_by(SiteID) %>% 
  mutate(diffbest = mean_lost - min(mean_lost)) %>% 
  filter(diffbest == 0) %>% 
  arrange(SiteID, cdl_mu) %>% 
  slice(1)

best_mm <- geom_lam %>% 
  filter(lambda == 1.5) %>% 
  group_by(SiteID) %>% 
  mutate(diffbest = abs(mean_mismatch) - min(abs(mean_mismatch))) %>% 
  filter(diffbest == 0) %>% 
  arrange(SiteID, cdl_mu) %>% 
  slice(1)

plot(jitter(best_cdl$cdl_mu), jitter(best_lost$cdl_mu), ylab = "CP minimizing lost gen", xlab = "CP maximizing lambda")
abline(0,1)
plot(jitter(best_cdl$cdl_mu), jitter(best_mm$cdl_mu), ylab = "CP minimizing mismatch", xlab = "CP maximizing lambda")
abline(0,1)
# 
# # extract voltinism/lost/attempted info for the row with best cdl
# for (site in unique(res$SiteID)){
#   best <- best_cdl %>% filter(SiteID == site)
#   tmp <- res %>% 
#     filter(SiteID == site, cdl_mu == best$cdl_mu, cdl_sd == best$cdl_sd, lambda == best$lambda)
#   
#   attempts <- tmp %>% 
#     tidyr::gather(gen, num, starts_with("att"))
#   
#   
# }


# site variables from gdd data
site_gdd <- gdd_all %>% 
  group_by(Site) %>% 
  summarise(mean_accumdd = mean(accumdegday[yday == 365]),
            early_accumdd = mean(accumdegday[yday == 170]),
            late_accumdd = mean(accumdegday[yday == 365] - accumdegday[yday == 170]),
            percgdd_early = early_accumdd / late_accumdd,
            lat = Latitude[1],
            lon = Longitude[1])

# map
dat <- left_join(best_cdl, site_gdd, by = c("SiteID" = "Site"))

# # gdd better predictor of cdl than latitude
summary(lm(cdl_mu ~ scale(mean_accumdd), data = dat))
summary(lm(cdl_mu ~ scale(lat), data = dat))

region <- rgdal::readOGR("./src/ref/ne_50m_admin_1_states_provinces_lakes", 'ne_50m_admin_1_states_provinces_lakes', encoding='UTF-8')
# region <-spTransform(region, CRS(proj4string(template)))
reg.points = fortify(region, region="name_en")
reg.df = left_join(reg.points, region@data, by = c("id" = "name_en"))
theme_set(theme_bw(base_size = 20))

xs <- range(dat$lon) + c(-.5, .5)
ys <- range(dat$lat) + c(-.5, .5)

tmpplt <- ggplot() +
  geom_polygon(data = reg.df, aes(x = long, y = lat, group = group), fill = NA, color = "light gray", inherit.aes = FALSE, size = 1, alpha = .3) +
  # geom_point(data = dat, aes(x = lon, y = lat, color = as.factor(round(mean_attvolt))), alpha = 1, size = 3, inherit.aes = FALSE) +
  # scale_color_viridis(discrete = TRUE, option = "C") +
  geom_point(data = dat, aes(x = lon, y = lat, color = mean_attvolt), alpha = 1, size = 4, inherit.aes = FALSE) +
  scale_color_viridis(discrete = FALSE, name = "Attempted") +
  coord_fixed(1.3, xlim = xs, ylim = ys, expand = FALSE, clip = "on") +
  # coord_fixed(1.3, xlim = c(-126.5, -119.5), ylim = c(36, 52), expand = FALSE, clip = "on") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
tmpplt



# adjust for lambda
lamlist <- list()
lamrange <- seq(.9, 3, length.out = 25)
for (i in 1:nrow(res)){
  tmp <- res[i, ]
  
  comp <- tmp %>% 
    dplyr::select(contains("comp")) %>% 
    tidyr::gather(gen, num) %>% 
    mutate(gen = as.numeric(str_split_fixed(gen, pattern = coll("_"), 2)[,2]))
  
  for (j in 1:length(lamrange)){
    lam <- comp %>% 
      filter(complete.cases(.)) %>% 
      mutate(lambda = (num * (lamrange[j]) ^ gen))
    
    lamlist[[(length(lamlist) + 1)]] <- tmp %>% 
      dplyr::select(Year, cdl_mu, cdl_sd, voltinism, lost) %>% 
      mutate(gen_lambda = lamrange[j],
             annual_lambda = sum(lam$lambda) / 100)
  }
}


lams <- bind_rows(lamlist)

plt <- ggplot(lams, aes(x = CDL, y = gen_lambda, fill = log(annual_lambda))) +
  geom_raster() +
  scale_x_reverse() +
  scale_fill_viridis(begin = 1, end = 0) +
  facet_wrap(~Year, ncol = 5) +
  theme_bw() +
  ylab("Between generation lambda") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated log(Annual Population Growth Rate)", subtitle = "Accounting for lost generations")
plt
ggsave(filename = paste0("GCA_AnnualLambda_", site, ".png"), plot = plt, device = "png", width = 9, height = 6, units = "in")



gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geom_lam <- lams %>% 
  # filter(Year > 2007) %>% 
  group_by(CDL, gen_lambda) %>% 
  tidyr::complete(Year, CDL, gen_lambda, fill = list(annual_lambda = 0)) %>% 
  summarise(mean_annual_lambda = gm_mean(annual_lambda))

plt <- ggplot(geom_lam, aes(x = CDL, y = gen_lambda, fill = log(mean_annual_lambda))) +
  geom_raster() +
  scale_x_reverse() +
  scale_fill_viridis(begin = 1, end = 0) +
  theme_bw() +
  ylab("Between generation lambda") +
  xlab("Critical photoperiod") +
  ggtitle("Simulated mean log(Annual Population Growth Rate)", subtitle = "1993-2017 geometric mean log(lambda)")
plt

ggsave(filename = paste0("GCA_meanlambda_", site, ".png"), plot = plt, device = "png", width = 9, height = 6, units = "in")






# GRAVEYARD----
# Compare to cohort model results

# UNKNOWN: Does voltinism from IBM == voltinism from cohort model?
# Cohort has weighted completed x weighted diapause, is this the same as
# weighted mean of [completed(1) x diapause(1), completed(2) x diapause(2), ...]?
# Go to cohort model and test by putting voltinism calc (1st for each cohort)
# in matrix of weighted results 

coh <- readRDS("site_array_gca_17coh.rds")

# extract gdd, voltinism, lost, 
cohlist <- list()
for (w in 1:53){
  
  coh_att <- data.frame(coh[, w, 7,])
  names(coh_att) <- paste0("X", 2016:2020)
  coh_att <- coh_att %>% 
    mutate(ID = sites$ID) %>% 
    tidyr::gather(year, attempted, X2016:X2020) %>% 
    mutate(year = gsub(pattern = "X", replacement = "", x = year, fixed = TRUE),
           week = w)
  
  
  coh_volt <- data.frame(coh[, w, 8,])
  names(coh_volt) <- paste0("X", 2016:2020)
  coh_volt <- coh_volt %>% 
    mutate(ID = sites$ID) %>% 
    tidyr::gather(year, completed, X2016:X2020) %>% 
    mutate(year = gsub(pattern = "X", replacement = "", x = year, fixed = TRUE),
           week = w)
  
  coh_wvolt <- data.frame(coh[, w, 9,])
  names(coh_wvolt) <- paste0("X", 2016:2020)
  coh_wvolt <- coh_wvolt %>% 
    mutate(ID = sites$ID) %>% 
    tidyr::gather(year, wvolt, X2016:X2020) %>% 
    mutate(year = gsub(pattern = "X", replacement = "", x = year, fixed = TRUE),
           week = w)
  
  coh_diap <- data.frame(coh[, w, 10,])
  names(coh_diap) <- paste0("X", 2016:2020)
  coh_diap <- coh_diap %>% 
    mutate(ID = sites$ID) %>% 
    tidyr::gather(year, diap, X2016:X2020) %>% 
    mutate(year = gsub(pattern = "X", replacement = "", x = year, fixed = TRUE),
           week = w)
  
  coh_all <-  coh_att %>% 
    left_join(coh_volt) %>% 
    left_join(coh_wvolt) %>% 
    left_join(coh_diap)
  cohlist[[w]] <- coh_all
}
cohplot <- bind_rows(cohlist) %>% 
  tidyr::gather(var, val, attempted, completed, wvolt, diap)

# plot helps show that diapause * completed doesn't make sense, if 

ggplot(cohplot, aes(x = week, y = val, group = ID, color = ID)) +
  geom_line() +
  facet_grid(var ~ year, scales = "free") +
  theme_bw(base_size = 14)

# best diapause cutoff to make IBM??
# cohplot <- bind_rows(cohlist)
# ggplot(cohplot, aes(x = week, y = ifelse(diap > .75, NA, 1) * completed, group = ID, color = ID)) +
#   geom_line(aes(size = 1-diap)) +
#   facet_grid(. ~ year) +
#   theme_bw(base_size = 14)
# 
# coh_volt <- cohplot %>% 
#   mutate(d25 = ifelse(diap > .25, NA, 1) * completed,
#          d50 = ifelse(diap > .5, NA, 1) * completed,
#          d75 = ifelse(diap > .75, NA, 1) * completed,
#          d90 = ifelse(diap > .9, NA, 1) * completed,
#          d95 = ifelse(diap > .95, NA, 1) * completed,
#          d98 = ifelse(diap > .98, NA, 1) * completed) %>% 
#   group_by(ID, year) %>% 
#   summarise_all(max, na.rm = TRUE) %>% 
#   arrange(year, ID)
# 
# # good correlation between the cohort and IBM voltinism estimates now
# # but what does it really mean? 
# # is this just weighted voltinism? Why does it work for the cohort model?
# coh_volt <- coh_volt %>% ungroup() %>%  mutate(res$voltinism)
# plot(coh_volt$`res$voltinism`, coh_volt$d95)
# abline(0, 1)

# this matches IBM well if attempted used instead of completed, why?
cohdf <- bind_rows(cohlist) %>% 
  group_by(ID, year) %>% 
  arrange(week) %>% 
  # mutate(diapnewgen = cumsum(c(0, diff(attempted)) * (1-diap))) %>% 
  filter(week == 53) %>% 
  ungroup() %>% 
  rename(SiteID = ID, Year = year) %>% 
  mutate(Year = as.numeric(Year),
         mismatch = wvolt - completed)

ind_diap <- res %>% 
  rowwise() %>% 
  mutate(ind_comp = sum(c(comp_1, comp_2, comp_3, comp_4), na.rm = TRUE),
         ind_diap = ind_comp / (ind_comp + lost)) %>% 
  left_join(cohdf)

plot(ind_diap$voltinism, ind_diap$wvolt)
abline(0, 1)
plot(ind_diap$ind_diap, ind_diap$diap)
abline(0, 1)

