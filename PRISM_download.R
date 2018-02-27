# PRISM download
# using ropensci 'prism' package to access webservice
# set download directory?


library(prism)
library(lubridate)

startdate <- "2017-10-01"
enddate <- "2017-12-31"
yr <- year(ymd(startdate))
options(prism.path = paste("prismDL", yr, sep = "/"))

get_prism_dailys("tmin", minDate = startdate, 
                 maxDate = enddate, keepZip = FALSE)
get_prism_dailys("tmax", minDate = startdate, 
                 maxDate = enddate, keepZip = FALSE)
