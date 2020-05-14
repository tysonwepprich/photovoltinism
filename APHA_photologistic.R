# Aphalara photoperiod logistic regression

library(ggplot2)
library(tidyverse)
library(broom)
library(MASS)
dat <- read.csv("data/Aphalara_diapause_R.csv") %>% 
  mutate(Total = Repro + Diap,
         perc_repro = Repro/Total,
         perc_diap = Diap/Total)

mod1 <- glm(perc_repro ~ Treat + Pop, weights = Total, data = dat, family = binomial(link = "probit"))
mod2 <- glm(perc_repro ~ Treat * Pop, weights = Total, data = dat, family = binomial(link = "probit"))
car::Anova(mod2, type = 2)

mods <- dat %>% 
  group_by(Pop) %>% 
  do(fit = glm(perc_diap ~ Treat, weights = Total, data = ., family = binomial(link = "probit")))

modcoefs <- bind_rows(lapply(mods$fit, FUN = function(x) data.frame(b0 = coef(x)[1], b1 = coef(x)[2]))) %>% 
  mutate(mods$Pop) %>% 
  mutate(cp_mean = -b0/b1,
         cp_sd = 1/abs(b1))

modcps <- bind_rows(lapply(mods$fit, FUN = function(x) {
  tmp <- MASS::dose.p(x, p = c(.05, .5, .95))
  data.frame(pred = as.vector(tmp),
             se = as.vector(attr(tmp, "SE")),
             p = as.vector(attr(tmp, "p")))
})) %>% 
  # mutate(mods$Pop) %>% 
  pivot_wider(names_from = p, values_from = c(pred, se))

# fakedat <- data.frame(Pop = c("North", "North", "South", "South"),
#                       Treat = c(12, 17, 12, 17),
#                       Repro = c(0, 1000, 0, 1000),
#                       Diap = c(1000, 0, 1000, 0),
#                       fraction = c(0, 1, 0, 1))
# dat <- rbind(dat, fakedat)

datN <- dat %>% filter(Pop == "North")
modN <- glm(cbind(Diap, Repro) ~ Treat, data = datN, family = binomial(link = "logit"))

predN <- data.frame(Treat = rep(seq(12.5, 16.5, length.out = 100), 1), Pop = "North")
predN$pred <- predict(modN, predN, type = "response")             


datS <- dat %>% filter(Pop == "South")
modS <- glm(cbind(Diap, Repro) ~ Treat, data = datS, family = binomial(link = "logit"))
# coefs <- coef(modS)
predS <- data.frame(Treat = rep(seq(12.5, 16.5, length.out = 100), 1), Pop = "South")
predS$pred <- predict(modS, predS, type = "response")       
# predS$pred2 <- exp(coefs[1]+coefs[2]*predS$Treat)/(1+exp(coefs[1]+coefs[2]*predS$Treat))

preds <- rbind(predN, predS)


# best model has additive Pop by AIC, not interaction or omitted
mod <- glm(cbind(Diap, Repro) ~ Treat + Pop - 1, data = dat, family = binomial(link = "probit"))
# coefs <- coef(modS)
pred <- data.frame(Treat = rep(seq(12.5, 16.5, length.out = 100), 2), Pop = c(rep("South", 100), rep("North", 100)))
pred$pred <- predict(mod, pred, type = "response")      

modcps <- bind_rows(lapply(c(3, 2), FUN = function(x) {
  tmp <- MASS::dose.p(mod, cf = c(x, 1), p = c(.05, .5, .95))
  data.frame(pred = as.vector(tmp),
             se = as.vector(attr(tmp, "SE")),
             p = as.vector(attr(tmp, "p")))
})) %>% 
  mutate(Pop = c(rep("S", 3), rep("N", 3))) %>%
  pivot_wider(names_from = p, values_from = c(pred, se))


plt <- ggplot(dat, aes(x = Treat, y = fraction, group = Pop, color = Pop)) +
  geom_point(size = 2) +
  geom_line(data = pred, aes(x = Treat, y = 1 - pred, group = Pop, color = Pop)) +
  xlab("Hours of daylight") +
  ylab("Percent reproductive") + 
  theme_bw(base_size = 16) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plt




findInt <- function(model, value) {
  function(x) {
    predict(model, data.frame(Treat=x), type="response") - value
  }
}

uniroot(findInt(modS, .5), range(datS$Treat))$root
uniroot(findInt(modS, .95), range(datS$Treat))$root
uniroot(findInt(modS, .05), range(datS$Treat))$root

