library(nlme)
# Fit the null model for Table 4.1, p. 47
data(bdf)
fm1 <- lme(langPOST ~ 1, data = bdf, random = ~ 1 | schoolNR)
VarCorr(fm1)
-2*logLik(fm1)
