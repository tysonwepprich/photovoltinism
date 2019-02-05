# biocontrol species parameters
# to source for model_lifecycle.R
species_params <- function(mod_type, species, biotype, nsim, model_CDL, dd_sd){
  
  # Galerucella calmariensis ----
  if (species == "GCA"){
    # life cycle parameters
    # life stage order starting with overwintering stages
    stgorder   <- c("OA","E","L","P", "TA", "A") # reassigned as 1, 2, 3, 4, 5 in model
    # Degree day thresholds
    # need to match length and order of stgorder
    # Could add biotype differences here if needed
    stage_ldt <- rep(10, 6)
    stage_udt <- rep(37.8, 6)
    stage_dd  <- c(167.6, 93.3, 136.4, 137.7, 125, 10)
    
    # diapause response to photoperiod
    photo_sens <- 5 # integer sensitive stages, can have more than one
    
    if (model_CDL == 0){
      coefs <- c("NA", "NA", "NA")
    }
    if (model_CDL == 1){
      coefs <- c(switch(biotype, 
                        "N" = 16.23,
                        "S" = 15.52), NA, NA)
    }
    # Logistic regression photoperiod
    # TODO: automatically check if regression models % diapause or reproductive
    # Preferred output is % diapause
    if (model_CDL == 2){
      coefs <- c(NA, switch(biotype,
                            "N" = c(56.9745, -3.5101),
                            "S" = c(60.3523, -3.8888)))
    }
    
    # Take distribution and calculate substages for oviposition distribution
    # TRY 1 with normal truncated at beginning of season
    # arg1 = list(mu = 97.94, sigma2 = 2241.7)
    # xdist = seq(50, 350, length.out = 1000)
    # ydist = dnorm(xdist, mean = arg1$mu, sd = sqrt(arg1$sigma2))
    
    # Try 2 with skewed t
    arg1 = list(mu = 97.94, sigma2 = 2241.7, shape = 3.92, nu = 9.57)
    xdist = seq(50, 350, length.out = 1000)
    ydist = mixsmsn:::dt.ls(xdist, loc = arg1$mu, sigma2 = arg1$sigma2, shape = arg1$shape, nu = arg1$nu)
    
    # plot(xdist, ydist)
    # TODO: empirical distribution from data
  } # end GCA 
  
  # Diorhabda carinulata ----
  if (species == "DCA"){
    # life cycle parameters
    # life stage order starting with overwintering stages
    stgorder   <- c("OA","E","L","P", "A") # reassigned as 1, 2, 3, 4, 5 in model
    # Degree day thresholds
    # need to match length and order of stgorder
    # Could add biotype differences here if needed
    stage_ldt <- rep(11.1, 5)
    stage_udt <- rep(36.7, 5)
    stage_dd  <- c(275, 95, 186, 188, 51)
    
    # diapause response to photoperiod
    photo_sens <- 5 # integer sensitive stages, can have more than one
    
    if (model_CDL == 0){
      coefs <- c("NA", "NA", "NA")
    }
    if (model_CDL == 1){
      coefs <- c(switch(biotype, 
                        "Big Bend" = 13.52,
                        "Delta" = 14.32,
                        "Gold Butte" = 14.04,
                        "Lovelock" = 14.24,
                        "Topock Marsh" = 12.75
      ), NA, NA)
    }
    # Logistic regression photoperiod
    # TODO: automatically check if regression models % diapause or reproductive
    # Preferred output is % diapause
    if (model_CDL == 2){
      coefs <- c(NA, switch(biotype,
                            "Original" = c(113.039752, -7.850763 ),
                            "Evolved" = c(31.384807, -2.461436),
                            "Big Bend" = c(53.8274, -3.9395),
                            "Delta" = c(102.6457, -7.1704),
                            "Gold Butte" = c(83.8342, -5.9676),
                            "Lovelock" = c(98.4441, -6.9109),
                            "Topock Marsh" = c(31.2262, -2.4481)
      )
      )
    }
    
    # Take distribution and calculate substages for oviposition distribution
    # TRY 1 with normal truncated at beginning of season
    arg1 = list(mu = 275, sigma2 = 1000)
    xdist = seq(150, 450, length.out = 1000)
    ydist = dnorm(xdist, mean = arg1$mu, sd = sqrt(arg1$sigma2))
    
    # plot(xdist, ydist)
    # TODO: empirical distribution from data
  } # end DCA
  
  # Aphalara itadori ----
  if (species == "APHA"){
    # life cycle parameters
    # life stage order starting with overwintering stages
    stgorder   <- c("OA","E","L", "TA", "A") # reassigned as 1, 2, 3, 4 in model
    # Degree day thresholds
    # need to match length and order of stgorder
    # Could add biotype differences here if needed
    stage_ldt <- rep(6.9, 5)
    stage_udt <- rep(32, 5)
    stage_dd  <- c(306, 136, 409, 75, 10)
    
    # diapause response to photoperiod
    photo_sens <- 4 # integer sensitive stages, can have more than one
    
    if (model_CDL == 0){
      coefs <- c("NA", "NA", "NA")
    }
    if (model_CDL == 1){
      coefs <- c(switch(biotype, 
                        "N" = 14.73,
                        "S" = 14.12), NA, NA)
    }
    # Logistic regression photoperiod
    # TODO: automatically check if regression models % diapause or reproductive
    # Preferred output is % diapause
    if (model_CDL == 2){
      coefs <- c(NA, switch(biotype,
                            "N" = c(90.0916, -6.1151),
                            "S" = c(47.0957, -3.3456)))
    }
    
    # Take distribution and calculate substages for oviposition distribution
    # Very little known about emergence
    # TRY 1 with normal distrib
    arg1 = list(mu = 220, sigma2 = 2500)
    xdist = seq(90, 350, length.out = 1000)
    ydist = dnorm(xdist, mean = arg1$mu, sd = sqrt(arg1$sigma2))
    
    # plot(xdist, ydist)
    # TODO: empirical distribution from data
  } # end APHA 
  
  
  if (mod_type == "cohort"){
    # for any species
    # output parameter file with
    # substages and photoperiod response
    inputdist <- data.frame(x = xdist, y = ydist) %>% 
      arrange(x) %>% 
      mutate(CDF = cumsum(y/sum(y)))
    substages <- SubstageDistrib(dist = inputdist, numstage = nsim, perc = .99)
    
    # parameters of required degree-days and coefficients of photoperiod response model
    if (nsim == 1){
      ddpar <- matrix(stage_dd, nrow = 1, byrow = TRUE)
    }else{
      ddpar <- cbind(substages$means, matrix(rep(stage_dd[-1], nrow(substages)), nrow = nrow(substages), byrow = TRUE))
    }
    
    params <- list(
      stgorder = stgorder,
      relpopsize = substages$weights,
      stage_dd = ddpar,
      stage_ldt = stage_ldt,
      stage_udt = stage_udt,
      photo_sens = photo_sens,
      CDL = coefs)
    
    return(params)
  }
  
  
  if(mod_type == "ibm"){
    
    # for any species
    # output parameter file with
    emerg <- base::sample(xdist, size = nsim, prob = ydist)
    
    # parameters of required degree-days and coefficients of photoperiod response model
    if (nsim == 1){
      ddpar <- matrix(stage_dd, nrow = 1, byrow = TRUE)
    }else{
      ddpar <- cbind(emerg, matrix(rep(stage_dd[-1], nsim), nrow = nsim, byrow = TRUE))
    }
    
    if(dd_sd > 0){
      ddpar[, 2:ncol(ddpar)] <- apply(ddpar[, 2:ncol(ddpar)], MARGIN = 2, FUN = function(x){rnorm(length(x), mean = mean(x), sd = dd_sd * mean(x) / 100)})
    }
    
    params <- list(
      stgorder = stgorder,
      stage_dd = ddpar,
      stage_ldt = stage_ldt,
      stage_udt = stage_udt,
      photo_sens = photo_sens,
      CDL = coefs)
    
    return(params)
    
  }
}