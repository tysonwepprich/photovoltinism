# biocontrol species parameters
# to source for ibm_model_optimal_cdl.R
species_params <- function(mod_type, species, biotype, nsim, model_CDL, dd_sd){
  
  # Galerucella calmariensis ----
  if (species == "GCA"){
    # life cycle parameters
    # life stage order starting with overwintering stages
    stgorder   <- c("OA","E","L","P", "TA", "A") # reassigned as 1, 2, 3, 4, 5 in model
    # Degree day thresholds
    # need to match length and order of stgorder
    # Could add biotype differences here if needed
    stage_ldt <- rep(12.2, 6)
    stage_udt <- rep(30, 6)
    stage_dd  <- c(100, 87.8, 128.2, 126.0, 72.9, 50)
    
    # diapause response to photoperiod
    photo_sens <- 5 # integer sensitive stages, can have more than one
    
    if (model_CDL == 0){
      coefs <- c("NA", "NA", "NA")
    }
    if (model_CDL == 1){
      coefs <- c(switch(as.character(biotype), 
                        "Bellingham" = c(16.6, 1),
                        "Yakima" = c(15.9, 1),
                        "Sutherlin" = c(14.5, 1),
                        "Palermo" = c(15.1, 1),
                        "Montesano" = c(16.2, 1),
                        "Ft. Drum" = c(15.4, 1),
                        "West Point" = c(15.85, 1)))
    }
    

    # Try 3 with skew normal
    xdist <- seq(50, 350, length.out = 1000)
    ydist = fGarch::dsnorm(xdist, mean = 120, sd = 25, xi = 3)
    
    # plot(xdist, ydist)
    # TODO: empirical distribution from data
  } # end GCA 
  
  # Diorhabda carinulata ----
  if (species == "DCA"){
    # life cycle parameters
    # life stage order starting with overwintering stages
    stgorder   <- c("OA","E","L","P", "TA", "A") # reassigned as 1, 2, 3, 4, 5 in model
    # Degree day thresholds
    # need to match length and order of stgorder
    # Could add biotype differences here if needed
    stage_ldt <- rep(12.0, 6)
    stage_udt <- rep(40, 6)
    stage_dd  <- c(75, 91.4, 176.7, 174.0, 47.2, 50)
    
    # diapause response to photoperiod
    photo_sens <- 5 # integer sensitive stages, can have more than one
    
    if (model_CDL == 0){
      coefs <- c("NA", "NA", "NA")
    }
    if (model_CDL == 1){
      coefs <- c(switch(as.character(biotype), 
                        "Big Bend" = c(11.94, 1.00),
                        "Delta" = c(14.44, 0.48),
                        "St. George" = c(13.46, 0.40),
                        "Lovell" = c(14.72, 0.32),
                        "Imperial" = c(11.19, 0.58),
                        "Ft. Carson" = c(14.23, 1),
                        "YumaPG" = c(11.19, 1),
                        "Pinon Canyon" = c(14.23, 1)))
    }
    
    # Take distribution and calculate substages for oviposition distribution
    # TRY 1 with normal truncated at beginning of season
    # arg1 = list(mu = 200, sigma2 = 1000)
    # xdist = seq(120, 350, length.out = 1000)
    # ydist = dnorm(xdist, mean = arg1$mu, sd = sqrt(arg1$sigma2))
    
    # Try 3 with skew normal
    xdist <- seq(50, 350, length.out = 1000)
    ydist = fGarch::dsnorm(xdist, mean = 100, sd = 25, xi = 3)
    # plot(xdist, ydist)
  } # end DCA
  
  # Aphalara itadori ----
  if (species == "APH"){
    # life cycle parameters
    # life stage order starting with overwintering stages
    stgorder   <- c("OA","E","L", "TA", "A") # reassigned as 1, 2, 3, 4 in model
    # Degree day thresholds
    # need to match length and order of stgorder
    # Could add biotype differences here if needed
    stage_ldt <- rep(6.9, 5)
    stage_udt <- rep(30, 5)
    stage_dd  <- c(220, 413, 133, 70, 50) 
    # adjusted to match manuscript
    # stages: OW emergence mean, egg + Nymph 1-4, Nymph 5 (sensitive), pre-ovip, ovip delay
    
    # diapause response to photoperiod
    photo_sens <- 3 # integer sensitive stages, can have more than one
    
    if (model_CDL == 0){
      coefs <- c("NA", "NA", "NA")
    }
    if (model_CDL == 1){
      coefs <- c(switch(as.character(biotype), 
                        "N" = c(14.7, 0.284),
                        "S" = c(14.1, 0.284)))
    }
    
    # Take distribution and calculate substages for oviposition distribution
    # Very little known about emergence
    # TRY 1 with normal distrib
    arg1 = list(mu = 220, sigma2 = 2500)
    xdist = seq(90, 350, length.out = 1000)
    ydist = dnorm(xdist, mean = arg1$mu, sd = sqrt(arg1$sigma2))
    
    # plot(xdist, ydist)
  } # end APHA 
  
  
  if (mod_type == "cohort"){
    # for any species
    # output parameter file with
    # substages and photoperiod response
    inputdist <- data.frame(x = xdist, y = ydist) %>% 
      arrange(x) %>% 
      mutate(CDF = cumsum(y/sum(y)))
    substages <- SubstageDistrib(dist = inputdist, numstage = nsim, perc = .999)
    
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
    emerg <- base::sample(xdist, size = nsim, prob = ydist, replace = TRUE)
    
    # parameters of required degree-days and coefficients of photoperiod response model
    if (nsim == 1){
      ddpar <- matrix(stage_dd, nrow = 1, byrow = TRUE)
    }else{
      ddpar <- cbind(emerg, matrix(rep(stage_dd[-1], nsim), nrow = nsim, byrow = TRUE))
    }
    
    if(dd_sd > 0){
      ddpar[, 2:(ncol(ddpar)-1)] <- apply(ddpar[, 2:(ncol(ddpar)-1)], MARGIN = 2, FUN = function(x){rnorm(length(x), mean = mean(x), sd = dd_sd * mean(x) / 100)})
      ddpar[, ncol(ddpar)] <- rpois(nrow(ddpar), lambda = 15)
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