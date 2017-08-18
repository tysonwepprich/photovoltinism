# photoperiod function from Forsythe et al. 1995
# matches 'geosphere' package output if no twilight

lat <- 44.56
doy <- 226
# p <- 0.8333 # only option in 'geosphere' package
p <- 6 # civil twilight as in Taylor 1986 recommendation
# p <- 1.5 # 25% twilight as in Grevstad & Coop 2015

photoperiod <- function(lat, doy, p){
  theta <- 0.2163108 + 
    2 * atan(0.9671396 * tan(0.00860 * (doy - 186)))
  phi <- asin(0.39795 * cos(theta))
  D <- 24 - (24 / pi) * acos(
    (sin(p * pi / 180) + sin(lat * pi / 180) * sin(phi))/
      (cos(lat * pi / 180) * cos(phi))
  )
}


# testing logistic regression

df <- data.frame(daylength = seq(13, 16, .2),
                 diapause = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1))
mod <- glm(diapause ~ daylength, data = df, family = binomial(link = "logit"))

newdat <- data.frame(daylength = seq(min(df$daylength), max(df$daylength), len=100))
newdat$diapause = predict(mod, newdata=newdat, type="response")
plot(diapause~daylength, data=df, col="red4")
lines(diapause~daylength, newdat, col="green4", lwd=2)

# inverse prediction at 50% diapause
# Dan Bean did this in JMP
doses <- MASS::dose.p(mod, cf = c(1,2), p = seq(0.1, .9, .01))

predict(mod, newdata = data.frame(daylength = 14), type = "link", se.fit = TRUE)
predict(mod, newdata = data.frame(daylength = 14), type = "response", se.fit = TRUE)

X = rnorm(n = 1000, mean = 14.5, sd = .3) 
P = ecdf(X)   
plot(P) 


# from Broatch

LogEq <- function(x){
  a <- 0
  b <- .05
  c <- 1
  m <- 400
  y <- a + c / (1 + exp(-b * (x - m)))
}

x <- seq(200, 800, 5)
plot(x, LogEq(x))


sims <- rnorm(1000, 500, 50)
isims <- boot::inv.logit(sims)




