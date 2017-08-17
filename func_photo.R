# photoperiod function from Forsythe
# matches 'geosphere' package if no twilight

lat <- 44.56
doy <- 226
# p <- 0.8333 # only option in 'geosphere' package
p <- 6 # civil twilight as in Grevstad & Coop 2015
theta <- 0.2163108 + 
  2 * atan(0.9671396 * tan(0.00860 * (doy - 186)))
phi <- asin(0.39795 * cos(theta))
D <- 24 - (24 / pi) * acos(
  (sin(p * pi / 180) + sin(lat * pi / 180) * sin(phi))/
    (cos(lat * pi / 180) * cos(phi))
)
