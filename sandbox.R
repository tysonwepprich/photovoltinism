


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
