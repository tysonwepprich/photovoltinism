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

'%!in%' <- function(x,y)!('%in%'(x,y))


plot_sim <- function(results, gdd){
  # on hold: selecting plot years after simulation
  # if(is.null(years)){
  #   results <- results
  #   gdd <- gdd
  # }else{
  #   results <- results %>% 
  #     filter(year %in% years)
  #   gdd <- gdd %>% 
  #     filter(year %in% years)
  # }
  # Marginal histograms ----
  mean_cdl <- mean(results$cdl, na.rm = TRUE)
  maxgdd <- max(gdd$accumdegday, na.rm = TRUE)
  maxphoto <- max(c(max(gdd$daylength, na.rm = TRUE), max(results$cdl, na.rm = TRUE)))
  minphoto <- min(gdd$daylength, na.rm = TRUE)
  
  cdl_hist1 <- hist(results$cdl[-which(is.na(results$cdl))], breaks = 20, plot = FALSE)
  cdl_hist <- data.frame(daylength = cdl_hist1$breaks[1:(length(cdl_hist1$breaks)-1)],
                         n = cdl_hist1$counts) 
  cdl_hist <- cdl_hist %>% 
    mutate(accumdegday = -50 + 50 * n/max(n))
  cdl_interval <- cdl_hist$daylength[2] - cdl_hist$daylength[1]
  
  sens_hist1 <- hist(results$accumdegday[which(is.na(results$accumdegday) == FALSE)], breaks = 100, plot = FALSE)
  sens_hist <- data.frame(accumdegday = sens_hist1$breaks[1:(length(sens_hist1$breaks)-1)],
                          n = sens_hist1$counts) 
  sens_hist <- sens_hist %>% 
    mutate(daylength = minphoto - .5 + .5 * n/max(n))
  sens_interval <- sens_hist$accumdegday[2] - sens_hist$accumdegday[1]
  
  
  # Photothermograph ----
  theme_set(theme_bw(base_size = 16)) 
  
  
  plt <- ggplot(gdd, aes(x = accumdegday, y = daylength, group = as.factor(year), color = as.factor(year))) +
    geom_line(size = 1) +
    coord_cartesian(xlim = c(-50, maxgdd + 100), ylim = c(minphoto - .5, maxphoto), expand = FALSE) +
    scale_color_viridis(discrete = TRUE, name = "Year", option = "D") +
    guides(color = guide_legend(reverse=FALSE)) +
    geom_rect(data=cdl_hist, show.legend=FALSE,
              aes(xmin = -50, xmax = accumdegday, ymin = daylength, ymax = daylength+cdl_interval), inherit.aes = FALSE) +
    geom_rect(data=sens_hist, show.legend=FALSE,
              aes(xmin=accumdegday, xmax=accumdegday + sens_interval, ymin=minphoto - .5, ymax=daylength), inherit.aes = FALSE) +
    geom_hline(yintercept = mean_cdl, linetype = "dashed", size = 1, alpha = .5) +
    ggtitle("Photothermographs with critical photoperiod\nand simulated sensitive stage emergence") +
    xlab("Accumulated degree-days") +
    ylab("Daylength (hours)") +
    theme(legend.position = c(0.9, 0.85))
  
  plt
  
  
}


plot_conseq <- function(results, gdd){
  maxgdd <- max(gdd$accumdegday, na.rm = TRUE)
  
  # Voltinism cost ----
  conseq <- results %>%
    mutate(consequence = ifelse(reproduce == TRUE & fit_diapause_gdd == TRUE & diapause_before_frost == TRUE,
                                "reproduce, F1 fits",
                                # ifelse(reproduce == TRUE & fit_diapause_gdd == TRUE & diapause_before_frost == FALSE,
                                #        "reproduce, F1 dies from frost",
                                ifelse(reproduce == TRUE & (fit_diapause_gdd == FALSE | diapause_before_frost == FALSE),
                                       "reproduce, F1 doesn't fit",
                                       "diapause")),
           date = as.Date(yday, origin=as.Date("2015-12-31")))
  
  
  
  
  plt2 <- ggplot(conseq, aes(x = accumdegday, y = consequence, 
                             fill = consequence, height = ..count..)) +
    # geom_density_ridges(scale = .9, rel_min_height=0.01, stat = "density") +
    geom_density_ridges(scale = .95, rel_min_height=0.001, stat = "binline", bins = 100) +
    scale_fill_viridis(discrete = TRUE, alpha = .5) +
    # scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
    scale_y_discrete(expand = c(0.1, .1)) +
    ggtitle("Diapause decisions and consequences") +
    xlab("Accumulated degree-days") +
    ylab("Frequency of simulations") +
    coord_cartesian(xlim = c(-50, maxgdd + 100), expand = FALSE) +
    theme(legend.position = c(0.75, 0.9),
          axis.text.y=element_blank()) +
    guides(fill = guide_legend(reverse=T))
  plt2
  
  
  # plt3 <- ggplot(conseq, aes(x = date, y = consequence, 
  #                            fill = consequence, height = ..count..)) +
  #   # geom_density_ridges(scale = .9, rel_min_height=0.01, stat = "density") +
  #   geom_density_ridges(scale = .95, rel_min_height=0.001, stat = "binline", bins = 50) +
  #   scale_fill_viridis(discrete = TRUE, alpha = .5) +
  #   scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  #   # scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
  #   scale_y_discrete(expand = c(0.1, .1)) 
  # # facet_wrap(~year, nrow = 2)
  # plt3
  
  
}
