# PRISM download
# using ropensci 'prism' package to access webservice
# set download directory?


library(prism)
options(prism.path = "prismDL")

startdate <- "2017-01-01"
enddate <- "2017-10-21"
get_prism_dailys("tmax", minDate = startdate, 
                 maxDate = enddate, keepZip = FALSE)

