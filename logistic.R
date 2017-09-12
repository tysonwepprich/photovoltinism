logit <- function(p){log(p/(1-p))}
inverse.logit <- function(or){1/(1 + exp(-or))}

set.seed(601)
# Let's generate some fake data to work with
n <- 100
x <- sort(rnorm(n))
# Calculate the 'true' value for the link
link <- .1 + 2*x
# Calculate the 'true' probability for y
# Using inverse logit transformation
p <- inverse.logit(link)
# Generate y
y <- rbinom(n, 1, p)
# Stick x and y in data frame
dat <- data.frame(x = x, y = y)

# Fit model
o <- glm(y ~ x, data = dat, family = binomial)


# Find true x for p = .7
true.x <- (logit(.7) - .1)/2

# estimated x for p = .7
estimated.x <- unname((logit(.7) - coef(o)[1])/coef(o)[2])

library(msm) # needed for deltamethod function

p <- .7
standerr <- deltamethod(~ (log(p/(1-p)) - x1)/x2, coef(o), vcov(o))
standerr
ci <- estimated.x + c(-1, 1) * 1.96 * standerr
c(true = true.x, estimated = estimated.x, ci = ci)
