library(nlme)
data(Oxboys)
fm1 <- lmList(Oxboys)
fm1
fm2 <- lme(fm1)
fm2

# bug report from Arne.Mueller@sanofi-aventis.com
mod <- distance ~ age + Sex
fm3 <- lme(mod, Orthodont, random = ~ 1)
predict(fm3, Orthodont)
