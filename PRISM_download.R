# PRISM download
# using ropensci 'prism' package to access webservice
# set download directory?


library(prism)
startdate <- "2015-01-01"
enddate <- "2015-12-31"
get_prism_dailys("tmin", minDate = startdate, 
                 maxDate = enddate, keepZip = FALSE)

