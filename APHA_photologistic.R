# Aphalara photoperiod logistic regression

library(ggplot2)
library(dplyr)
dat <- read.csv("data/Aphalara_diapause_R.csv")

mod0 <- glm(cbind(Diap, Repro) ~ Treat, data = dat, family = binomial(link = "logit"))
mod1 <- glm(cbind(Diap, Repro) ~ Treat + Pop, data = dat, family = binomial(link = "logit"))
mod2 <- glm(cbind(Diap, Repro) ~ Treat * Pop, data = dat, family = binomial(link = "logit"))
AIC(mod0, mod1, mod2)
# fakedat <- data.frame(Pop = c("North", "North", "South", "South"),
#                       Treat = c(12, 17, 12, 17),
#                       Repro = c(0, 1000, 0, 1000),
#                       Diap = c(1000, 0, 1000, 0),
#                       fraction = c(0, 1, 0, 1))
# dat <- rbind(dat, fakedat)

datN <- dat %>% filter(Pop == "North")
modN <- glm(cbind(Diap, Repro) ~ Treat, data = datN, family = binomial(link = "logit"))

predN <- data.frame(Treat = rep(seq(12, 17, length.out = 100), 1), Pop = "North")
predN$pred <- predict(modN, predN, type = "response")             


datS <- dat %>% filter(Pop == "South")
modS <- glm(cbind(Diap, Repro) ~ Treat, data = datS, family = binomial(link = "logit"))
# coefs <- coef(modS)
predS <- data.frame(Treat = rep(seq(12, 17, length.out = 100), 1), Pop = "South")
predS$pred <- predict(modS, predS, type = "response")       
# predS$pred2 <- exp(coefs[1]+coefs[2]*predS$Treat)/(1+exp(coefs[1]+coefs[2]*predS$Treat))

preds <- rbind(predN, predS)


# best model has additive Pop by AIC, not interaction or omitted
mod <- glm(cbind(Diap, Repro) ~ Treat + Pop, data = dat, family = binomial(link = "logit"))
# coefs <- coef(modS)
pred <- data.frame(Treat = rep(seq(12, 17, length.out = 100), 2), Pop = c(rep("South", 100), rep("North", 100)))
pred$pred <- predict(mod, pred, type = "response")      



plt <- ggplot(dat, aes(x = Treat, y = 1 -fraction, group = Pop, color = Pop)) +
  geom_point() +
  geom_line(data = pred, aes(x = Treat, y = pred, group = Pop, color = Pop))
plt




findInt <- function(model, value) {
  function(x) {
    predict(model, data.frame(Treat=x), type="response") - value
  }
}

uniroot(findInt(modS, .5), range(datS$Treat))$root
uniroot(findInt(modS, .95), range(datS$Treat))$root
uniroot(findInt(modS, .05), range(datS$Treat))$root

