# 1. Setup -------
# packages, options, functions loaded
# TODO: find ways to take out dplyr, purrr, mixsmsn functions?
pkgs <- c("lubridate", "daymetr",
          "stringr", "purrr", 
          "ggplot2", "dplyr", "tidyr", "ggpubr")
# install.packages(pkgs) # install if needed
inst = lapply(pkgs, library, character.only = TRUE) # load packages

# load collection of functions for this model
source('CDL_funcs.R')
source('species_params_20220114.R')
theme_set(theme_bw(base_size = 14) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())) 


# 2. User input -----
spparam <- read.csv("phototherm_params_20220114.csv", header = TRUE)

for (i in 1:nrow(spparam)){
# Species Specific Phenology Model Parameters:

# model scope
# yr           <- 2018
start_doy    <- 1
end_doy      <- 365
species      <- spparam$species[i] # GCA/APH/DCA
biotype      <- spparam$biotype[i] 
nsim <- 1 # can be more for individual variation ~100
sss <- spparam$site[i]
startyear <- 2010
endyear <- 2020

# Load sites
# sites <- read.csv("data/DCA_sitecoords.csv", header = TRUE)[,c(2, 4, 5)]
# sites <- read.csv("data/GCA_modeling_sites.csv", header = TRUE)[1:3,]
# 

# photoperiod decision inclusion
# 1 for single value CDL, 0 for none
model_CDL  <- 1


# Derived parameters -----
params <- species_params(species, biotype, nsim, model_CDL)
# params <- species_params(mod_type = "ibm", species, biotype, nsim, model_CDL, dd_sd = 0)

ldt <- params$stage_ldt[1]
udt <- params$stage_udt[1]

# Download weather data ----

# for dataframe of sites
# siteslist <- list()
# for(i in 1:nrow(sites)){
  temp <- download_daymet(site = sss, lat = spparam$y[i], lon = spparam$x[i],
                          start = startyear, end = endyear, internal = TRUE,
                          silent = TRUE, force = FALSE)
  outdf <- temp$data %>%
    mutate(elev = temp$altitude,
           Site = sss,
           Latitude = spparam$y[i],
           Longitude = spparam$x[i])
  # siteslist[[i]] <- outdf
# }
gdd <- outdf
# gdd <- bind_rows(siteslist)


# Calculate degree days and photoperiod ----

# Degree-day functions: 
# TriDD: horizontal cutoff
# TriDDvert: vertical cutoff
# degree_days: from TrenchR package, no option for vertical cutoff

gdd_all <- gdd %>%
  dplyr::filter(yday <= end_doy & yday >= start_doy) %>%
  rowwise() %>%
  mutate(degdayOLD = TriDD(tmax..deg.c., tmin..deg.c., ldt, udt),
         degdayTRI = degree_days(tmin..deg.c., tmax..deg.c., ldt, udt, method = "single.triangulation"),
         degday = TriDDvert(tmax..deg.c., tmin..deg.c., ldt, udt)) %>%
  ungroup() %>%
  group_by(Site, year) %>%
  arrange(yday) %>%
  mutate(accumdegdayOLD = cumsum(degdayOLD),
         accumdegdayTRI = cumsum(degdayTRI),
         accumdegday = cumsum(degday),
         daylength = photoperiod(Latitude, yday),
         daylength_daymet = dayl..s. / 60 / 60,
         lightfrost = tmin..deg.c. <= 0,
         hardfrost = tmin..deg.c. <= -2)


# naive max voltinism by gdd alone
ow_dd <- params$stage_dd[,1]
gen_dd <- sum(params$stage_dd[-1])

# calculated in different ways, average away individual variation or not?
maxvolt <- gdd_all %>% 
  dplyr::filter_all(all_vars(!is.na(.)))  %>% 
  group_by(Site, year) %>% 
  arrange(yday) %>%
  mutate(nday = n(),
         accumdegday = cumsum(degday)) %>% 
  dplyr::filter(yday == max(yday),
                nday > 350) %>% 
  summarise(maxgen = floor((accumdegday - mean(ow_dd)) / mean(gen_dd)))

# # Non-integer output of maxgen accounting for ind. variation
# maxvolt <- gdd_all %>% 
#   dplyr::filter(yday == 365) %>% 
#   summarise(maxgen = mean(floor((accumdegday - ow_dd) / gen_dd)),
#             accumdegday = accumdegday)

# Plot: Voltinism by degree-days (constant) ----
# example years 2011,13,15 for cool to warm
# G&C2015 uses 100dd OW ovip, 497 1st gen, 523 subsequent adult emergence

expand_df <- function(df){
  newdf <- data.frame(gen = 0:df$maxgen) %>% 
    mutate(ddreq = median(ow_dd) + gen * median(gen_dd)) %>% 
    rowwise() %>% 
    mutate(accumdegday = gdd_all %>% 
             dplyr::filter(Site == df$Site[1] & year == df$year[1]) %>% 
             dplyr::filter(accumdegday - ddreq[1] >= 0) %>% 
             arrange(yday) %>% 
             slice(1) %>% pull(accumdegday),
           yday = gdd_all %>% 
             dplyr::filter(Site == df$Site[1] & year == df$year[1]) %>% 
             dplyr::filter(accumdegday - ddreq[1] >= 0) %>% 
             arrange(yday) %>% 
             slice(1) %>% pull(yday))
  return(newdf)
}

milestones <- maxvolt %>% 
  # dplyr::filter(year %in% c(2011, 2013, 2015)) %>% 
  group_by(Site, year) %>% 
  do(expand_df(.))
milestones$year <- as.character(milestones$year)
milestones$date <- as.Date(milestones$yday, origin=as.Date("2015-12-31"))




dds <- gdd_all %>% 
  # dplyr::filter(year %in% c(2011, 2013, 2015)) %>% 
  mutate(date = as.Date(yday, origin=as.Date("2015-12-31")))
dds$year <- as.character(dds$year)

glab <- milestones %>% 
  group_by(Site, gen) %>% 
  summarise(accumdegday = mean(accumdegday),
            date = min(dds$date) + 10,
            lab = paste("Gen", gen, sep = " ")[1]) 

plt <- ggplot(dds, aes(x = date, y = accumdegday, group = year, color = year)) +
  geom_line(size = 1.5) +
  scale_color_brewer(name = "Year", type = "qual", palette = "Paired") + 
  # scale_fill_manual(name = "Generations",
  #                   values = Turbo(out.colors = 13))
  geom_point(data = milestones, aes(x = date, y = accumdegday), color = "black") +
  geom_hline(data = glab, aes(yintercept = accumdegday), linetype = "dotted") +
  geom_label(data = glab, aes(x = date, y = accumdegday, label = lab), fill = "white", inherit.aes = FALSE) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", limits = c(as.Date("2016-01-01"), as.Date("2016-12-01"))) +
  facet_wrap(~Site, nrow = 3) +
  xlab("Calendar date") + ylab("Accumulated degree-days")
# plt
# saved at 850x700 pixels
fname <- paste("gdd10yr", species, milestones$Site[1], "png", sep = ".")
ggsave(filename = fname, plot = plt, device = "png", width = 10, height = 8, units = "in")

# Plot: Voltinism by photothermograph (constant)----
# add photosensitive time to milestones
cdl <- spparam$cp_mean[i]
expand_df <- function(df){
  cycle <- "reproduce"
  newdf <- data.frame(gen = 1, start = median(ow_dd), end = median(ow_dd) + 1 * median(gen_dd),
                      sense = median(ow_dd) + 1 * median(gen_dd) - median(sum(params$stage_dd[c(params$photo_sens:length(params$stage_dd))]))) %>% 
    mutate(daylength = df %>% dplyr::filter(accumdegday - sense >= 0) %>% 
             arrange(yday) %>% 
             slice(1) %>% pull(daylength),
           choice = ifelse(daylength < cdl, "diapause", "reproduce"))
  while(cycle == "reproduce"){
    if(newdf$choice[nrow(newdf)] == "reproduce"){
      newrow <- newdf[nrow(newdf),] %>% 
        mutate(gen = gen + 1,
               start = end,
               end = end + median(gen_dd),
               sense = sense + median(gen_dd),
               daylength = ifelse(sense < max(df$accumdegday),
                                  df %>% dplyr::filter(accumdegday - sense >= 0) %>% 
                                    arrange(yday) %>% 
                                    slice(1) %>% pull(daylength), NA),
               choice = ifelse(daylength < cdl, "diapause", "reproduce"))
      newdf <- bind_rows(newdf, newrow)
    }else{
      cycle = "diapause"
    }
    if(newdf$end[nrow(newdf)] > max(df$accumdegday)){
      cycle = "lost"
      newdf$choice[nrow(newdf)] <- "lost"
    }
  }
  return(newdf)
}

events <- dds %>% 
  group_by(Site, year) %>% 
  do(expand_df(.)) %>% 
  mutate(ymin = as.numeric(year) - 0.5,
         ymax = as.numeric(year) + 0.5,
         shading = ifelse(choice == "lost", .4, 1),
         sense = ifelse(choice == "lost", NA, sense),
         choice = as.character(ifelse(choice == "lost", NA, choice)))

sens_time <- events %>% 
  dplyr::select(Site,year,sense) %>% 
  distinct()

# for (sss in sites$ID){

dds_site <- dds %>% filter(Site == sss)
events_site <- events %>% filter(Site == sss)
sens_time <- events_site %>% 
  dplyr::select(Site,year,sense) %>% 
  distinct()

maxdd <- max(c(dds_site$accumdegday, events$end)) + 10
xlabel <- paste0("Accumulated degree-days ", "(Base ", ldt, "C)")
cplab <- paste0("CP = ", cdl)

#photothermograph
plt1 <- ggplot(dds_site, aes(x = accumdegday, y = daylength, group = year, color = year)) +
  geom_line(size = 1.5) +
  scale_x_continuous(limits = c(0, maxdd)) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_color_brewer(name = "Year", type = "qual", palette = "Paired") + 
  geom_hline(aes(yintercept = cdl), linetype = "dashed") +
  annotate("text", x = maxdd - 50, y = cdl + .15, label = cplab, color = "black", size = 4) +
  geom_vline(data = sens_time, aes(xintercept = sense), linetype = "dotted") +
  theme(legend.position = "none") +
  ylab("Daylength (hours)") + xlab(xlabel)
# plt1

dd_simple <- dds_site %>% 
  ungroup() %>% 
  filter(yday %in% c(1, 365)) %>% 
  dplyr::select(Site, year, accumdegday) %>% 
  distinct() %>% 
  mutate(year = as.numeric(year))

genrect <- events_site %>% dplyr::select(-sense, -daylength, -choice)
genchoice <- events_site %>% na.omit

# bargraph of generations and diapause
plt2 <- ggplot(dd_simple, aes(x = accumdegday, y = year, group = year, color = year)) +
  scale_x_continuous(limits = c(0, maxdd)) +
  geom_rect(data = genrect, aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, 
                               fill = year, alpha = shading), color = "black", inherit.aes = FALSE) +
  scale_alpha_continuous(guide = "none") +
  scale_fill_brewer(guide = "none", name = "Year", type = "qual", palette = "Paired") + 
  geom_point(data = genchoice, aes(x = sense, y = (ymin + ymax)/2, shape = choice), color = "black", size = 3, inherit.aes = FALSE) +
  scale_shape_manual(guide = "none", name = "Choice", values=c(4, 1), na.translate=FALSE)+
  scale_y_continuous(breaks = c(2010:2020)) +
  # theme(legend.position = c(.2, .9)) +
  # legend.direction = "horizontal") +
  ylab("Year") + xlab(xlabel)
# plt2

site_plt <- ggarrange(plt1,plt2, 
               labels = c("A", "B"),
               ncol = 1, nrow = 2, align = "v", heights = c(2, 1))
fname <- paste("phototherm10yr", species, substr(sss, 1, 4), cdl, "png", sep = ".")
ggsave(filename = fname, plot = site_plt, device = "png", width = 8, height = 10, units = "in")

}

