# helper functions for photothermograph Shiny App


#### if then else raster function [sim. to GRASS r.mapcalc if(x,a,b)]:
Cond=function(condition, trueValue, falseValue){
  return(condition * trueValue + (!condition)*falseValue)
}



#### Single triangle with upper threshold (Sevachurian et al. 1977) - also a good substitution for 
#  single sine method
TriDD=function(tmax, tmin, LDT, UDT){
  Tmp1=6*((tmax-LDT)*(tmax-LDT))/(tmax-tmin)
  Tmp2=6*((tmax-UDT)*(tmax-UDT))/(tmax-tmin)
  Cond(tmax < LDT,0,
       Cond(tmin >= UDT,UDT-LDT,
            Cond((tmax < UDT) & (tmin <= LDT), Tmp1/12,
                 Cond((tmin <= LDT) & (tmax >= UDT), (Tmp1-Tmp2)/12,
                      Cond((tmin > LDT) & (tmax >= UDT), 6*(tmax+tmin-2*LDT)/12 - (Tmp2/12),
                           Cond((tmin > LDT) & (tmax < UDT), 6*(tmax+tmin-2*LDT)/12,0))))))
} 



# photoperiod functions
# default 1.5 degrees of twilight (6 maximum for civil twilight)
photoperiod <- function(lat, doy, p = 1.5){
  theta <- 0.2163108 + 
    2 * atan(0.9671396 * tan(0.00860 * (doy - 186)))
  phi <- asin(0.39795 * cos(theta))
  D <- 24 - (24 / pi) * acos(
    (sin(p * pi / 180) + sin(lat * pi / 180) * sin(phi))/
      (cos(lat * pi / 180) * cos(phi))
  )
}

# Download Daymet data and calculate photoperiod
get_temp_data <- function(userinput = NA){
  
  startyear <- userinput$year_range[1]
  endyear <- userinput$year_range[2]
  lat <- userinput$lat
  lon <- userinput$lon
  # Download weather data
  temp <- download_daymet(site = "default", lat = lat, lon = lon,
                          start = startyear, end = endyear, internal = TRUE,
                          silent = TRUE, force = FALSE)
  # Daylength for each unique day at site
  uniq_photo <- temp$data %>% 
    filter(yday < 350) %>% 
    dplyr::select(yday) %>% 
    distinct() %>% 
    mutate(daylength = photoperiod(lat, yday))
  
  out <- temp$data %>% 
    filter(yday < 350) %>% 
    left_join(uniq_photo, by = "yday")
  return(out)
}

get_gdd_data <- function(temp_data = NA, userinput = NA){
  ldt <- userinput$dev_temp[1]
  udt <- userinput$dev_temp[2]
  
  gdd <- temp_data %>% 
    rowwise() %>% 
    mutate(degday = TriDD(tmax..deg.c., tmin..deg.c., ldt, udt)) %>% 
    ungroup() %>% 
    group_by(year) %>% 
    arrange(yday) %>% 
    mutate(accumdegday = cumsum(degday))
  return(gdd)
}



SimLifecycle <- function(pars, yr, gdd){
  gdd <- gdd %>% filter(year == yr)
  maxgdd <- max(gdd$accumdegday, na.rm = TRUE)
  firstfrost <- min(gdd$accumdegday[which(gdd$tmin..deg.c. <= 0 & gdd$yday > 200)])
  lastfrost <- max(gdd$accumdegday[which(gdd$tmin..deg.c. <= 0 & gdd$yday < 180)])
  
  dflist <- list()
  
  # Simulate overwinter emergence
  emerg_dd <- rnorm(pars$nsim, mean = pars$emerg_dd, sd = pars$emerg_dd_sd)
  gen_dd <- rnorm(pars$nsim, mean = pars$gen_dd, sd = pars$gen_dd_sd)
  povip_dd <- rnorm(pars$nsim, mean = pars$povip_dd, sd = pars$povip_dd_sd)
  ovip_dd <- rpois(pars$nsim, lambda = pars$ovip_dd)
  genlength1 <- emerg_dd + gen_dd
  genlength2 <- gen_dd + povip_dd + ovip_dd
  df <- data.frame(index = seq_len(pars$nsim), sens_dd = emerg_dd, gen = 0, reproduce = TRUE) %>% 
    mutate(rowid = row_number()) %>% 
    group_by(rowid) %>% 
    mutate(chooserow = which.max(sens_dd <= gdd$accumdegday),
           accumdegday = gdd$accumdegday[chooserow],
           yday = gdd$yday[chooserow],
           daylength = gdd$daylength[chooserow],
           cdl = NA,
           diap_dd = genlength1[rowid],
           repro_dd = genlength2[rowid],
           fit_diapause_gdd = maxgdd >= (accumdegday + genlength1[rowid]),
           fit_repro_gdd = maxgdd >= (accumdegday + genlength2[rowid]),
           diapause_before_frost = firstfrost > (accumdegday + genlength1[rowid]),
           emerge_after_frost = lastfrost < sens_dd) %>% 
    ungroup() %>%
    dplyr::select(-rowid)
  dflist[[1]] <- df
  
  # Maximal potential voltinism without photoperiod
  maxvolt <- vector(length = length(emerg_dd))
  lastsensgdd <- vector(length = length(emerg_dd))
  
  for (i in seq_along(emerg_dd)){
    maxvolt[i] <- length(seq(from = genlength1[i],
                             to = maxgdd, 
                             by = genlength2[i]))
    lastsensgdd[i] <- max(seq(from = genlength1[i],
                              to = maxgdd, 
                              by = genlength2[i]))
  }
  
  voltinism <- data.frame(index = seq_len(pars$nsim), emerg_dd = emerg_dd, gen_dd = gen_dd, povip_dd = povip_dd,
                          ovip_dd = ovip_dd, maxvolt = maxvolt, maxgdd = maxgdd, firstfrost = firstfrost, lastfrost = lastfrost,
                          lastsensgdd = lastsensgdd)
  
  
  eggs <- ceiling(pars$lambda)
  mortality <- (eggs - pars$lambda) / eggs
  
  SimExpand <- function(dfrow){
    outdf <- dfrow[rep(seq_len(nrow(dfrow)), eggs), ]
    outdf$eggs <- seq_len(eggs)
    outdf$sens_dd <- Cond(outdf$gen[1] == 0, outdf$accumdegday + outdf$diap_dd,
                          outdf$accumdegday + outdf$repro_dd)
    outdf <- outdf %>% 
      rowwise() %>% 
      mutate(chooserow = which.max(sens_dd <= gdd$accumdegday),
             accumdegday = ifelse(chooserow == 1, NA, gdd$accumdegday[chooserow]),
             yday = ifelse(chooserow == 1, NA, gdd$yday[chooserow]),
             daylength = ifelse(chooserow == 1, NA, gdd$daylength[chooserow]),
             gen = gen + 1,
             cdl = rnorm(1, mean = pars$cdl_mu, sd = pars$cdl_sd),
             reproduce = daylength >= cdl,
             fit_diapause_gdd = maxgdd >= (accumdegday + diap_dd),
             fit_repro_gdd = maxgdd >= (accumdegday + repro_dd),
             diapause_before_frost = firstfrost > (accumdegday + diap_dd),
             emerge_after_frost = NA)
    
    return(outdf)
  }
  
  while(length(which(df$reproduce == TRUE & df$fit_repro_gdd == TRUE)) > 0){
    df <- df %>% 
      ungroup() %>% 
      filter(reproduce == TRUE & fit_repro_gdd == TRUE) %>% 
      group_by(row_number()) %>% 
      do(SimExpand(.)) %>% 
      ungroup() %>% 
      sample_frac(size = (1 - mortality), replace = FALSE)
    dflist[[length(dflist)+1]] <- df
  }
  
  simresults <- bind_rows(dflist) %>% 
    left_join(voltinism, by = "index") %>% 
    dplyr::select(-`row_number()`, -eggs, -chooserow)
  
  return(simresults)  
}


# run_sims <- function(nsim=NA, cdl_mu=NA, cdl_sd=NA, lambda=NA, emerg_dd=NA, 
#                      emerg_dd_sd=NA, gen_dd=NA, gen_dd_sd=NA, povip_dd=NA, 
#                      povip_dd_sd=NA, ovip_dd=NA, gdd = NA){
run_sims <- function(pars = NA, gdd = NA){
  # repackage parameters
  # pars <- data.frame(nsim = nsim, cdl_mu = cdl_mu, cdl_sd = cdl_sd,
  #                    lambda = lambda, emerg_dd = emerg_dd, emerg_dd_sd = emerg_dd_sd, 
  #                    gen_dd = gen_dd, gen_dd_sd = gen_dd_sd, povip_dd = povip_dd, 
  #                    povip_dd_sd = povip_dd_sd, ovip_dd = ovip_dd)
  pars <- data.frame(nsim = pars$nsim, cdl_mu = pars$cdl_mu, cdl_sd = pars$cdl_sd,
                     lambda = pars$lambda, emerg_dd = pars$emerg_dd, emerg_dd_sd = pars$emerg_dd_sd, 
                     gen_dd = pars$gen_dd, gen_dd_sd = pars$gen_dd_sd, povip_dd = pars$povip_dd, 
                     povip_dd_sd = pars$povip_dd_sd, ovip_dd = pars$ovip_dd)
  
  sims <- expand.grid(year = unique(gdd$year))
  results <- sims %>% 
    group_by(year) %>% 
    do(SimLifecycle(pars, .$year, gdd))
  return(results)
}


plot_sim <- function(results, gdd){
  
  # Marginal histograms ----
  mean_cdl <- mean(results$cdl, na.rm = TRUE)
  maxgdd <- max(gdd$accumdegday, na.rm = TRUE)
  maxphoto <- max(c(max(gdd$daylength, na.rm = TRUE), max(results$cdl, na.rm = TRUE)))
  minphoto <- min(gdd$daylength, na.rm = TRUE)
  
  cdl_hist <- results %>% 
    ungroup() %>% 
    dplyr::select(cdl) %>%
    filter(complete.cases(.)) %>% 
    mutate(breaks = cut(cdl, breaks=seq(min(cdl),max(cdl),length.out = 20), 
                        labels=seq(min(cdl),max(cdl), length.out = 20)[-1], 
                        include.lowest=TRUE)) %>% 
    mutate(daylength = as.numeric(as.character(breaks))) %>%
    group_by(daylength) %>% 
    summarise(n = n()) %>% 
    mutate(accumdegday = -50 + 50 * n/max(n))
  
  sens_hist <- results %>% 
    ungroup() %>% 
    dplyr::select(accumdegday) %>%
    filter(complete.cases(.)) %>% 
    mutate(breaks = cut(accumdegday, breaks = seq(min(accumdegday), max(accumdegday), length.out = 100),
                        labels = seq(min(accumdegday), max(accumdegday), length.out = 100)[-1],
                        include.lowest = TRUE)) %>% 
    mutate(accumdegday = as.numeric(as.character(breaks))) %>% 
    group_by(accumdegday) %>% 
    summarise(n = n()) %>% 
    mutate(daylength = minphoto - .5 + .5 * n/max(n))
  
  # Photothermograph ----
  theme_set(theme_bw(base_size = 16)) 
  
  
  plt <- ggplot(gdd, aes(x = accumdegday, y = daylength, group = as.factor(year), color = as.factor(year))) +
    geom_line(size = 1) +
    coord_cartesian(xlim = c(-50, maxgdd + 100), ylim = c(minphoto - .5, maxphoto), expand = FALSE) +
    scale_color_viridis(discrete = TRUE, name = "Year", option = "D") +
    guides(color = guide_legend(reverse=FALSE)) +
    geom_segment(data=cdl_hist, size=2.5, show.legend=FALSE,
                 aes(x=-50, xend=accumdegday, y=daylength, yend=daylength), inherit.aes = FALSE) +
    geom_segment(data=sens_hist, size=2.5, show.legend=FALSE,
                 aes(x=accumdegday, xend=accumdegday, y=minphoto - .5, yend=daylength), inherit.aes = FALSE) +
    geom_hline(yintercept = mean_cdl, linetype = "dashed", size = 1, alpha = .5) +
    # geom_rug(data = results, aes(x = accumdegday), color = "gray", alpha = .1) +
    # geom_rug(data = results, aes(y = cdl), color = "gray", alpha = .1) +
    ggtitle("Photothermographs with critical photoperiod\nand simulated sensitive stage emergence") +
    xlab("Accumulated degree-days") +
    ylab("Daylength (hours)")
  plt
  
  
}