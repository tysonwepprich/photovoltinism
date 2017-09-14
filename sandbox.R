# discretization of distributions
# simpler version

meanx <- 100
sdx <- 6
plot(pnorm(x = seq(.01, .99, .01), mean = meanx, sd = sdx))

SubstageVals <- function(numstage, perc){
  low <- qnorm((1 - perc)/2)
  high <- qnorm(1 - (1 - perc) / 2)
  bounds <- seq(low, high, length.out = numstage + 1)
  means <- (bounds[1:numstage] + bounds[2:(numstage + 1)]) / 2
  weights <- diff(pnorm(bounds), lag = 1)
  return(data.frame(means, weights))
}



# make rasters go faster
library(ff)
mat <- ff(vmode="double",dim=c(ncell(template), 20), 
                               filename=paste0(getwd(),"/stack.ffdata"))

test1 <- test2 <- test3 <- test4 <- test5 <- list()
for (i in 21:50){
  # current accumlating way
  if (!exists("teststack")){
    teststack <- stack(template)
  } else {
    test1[[i]] <- system.time({
      teststack <- stack(teststack, template)
    })
  }
  
  # addLayer
  if (!exists("teststack1")){
    teststack1 <- stack(template)
  } else {
    test2[[i]] <- system.time({
      teststack1 <- addLayer(teststack1, template)
    })
  }
  
  #shortcut
  if (!exists("teststack2")){
    teststack2 <- stack(template)
  } else {
    test3[[i]] <- system.time({
      teststack2@layers[[i]] <- template     
      })
  }
  
  # test4[[i]] <- system.time({
  #   mat[,i] <- template[]
  # })
  
  if (!exists("teststack3")){
    teststack3 <- writeRaster(stack(template), filename = "teststack3",
                              overwrite = TRUE)
  } else {
    test5[[i]] <- system.time({
      teststack3 <- addLayer(teststack3, template)
    })
  }
}
