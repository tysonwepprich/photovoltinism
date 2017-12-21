library(raster)
library(RCurl)

test <- scp(host = "ento:732", path = "~REPO/photovoltinism/results/gale/20072013/meanGDD_07_13.grd",
            user = "tyson")


system('scp -P 732 tyson@ento.hort.oregonstate.edu:')

scp your_username@remotehost.edu:foobar.txt /some/local/directory
