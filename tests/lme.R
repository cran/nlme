library(nlme)
data(Oxboys)
fm1 <- lmList(Oxboys)
fm1
fm2 <- lme(fm1)
fm2
