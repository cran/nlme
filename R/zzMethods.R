###   Miscellaneous methods that must be defined last in the library
###
### Copyright 1997-2003  Jose C. Pinheiro <Jose.Pinheiro@pharma.novartis.com>,
###                      Douglas M. Bates <bates@stat.wisc.edu>

## Note that  require( nls )  has already happened ...

BIC.lme <- BIC.lmList <- BIC.gls <- BIC.lm

comparePred.lme <- comparePred.lmList <- comparePred.gls

getData.nlme <- getData.gnls

getData.lme <- getData.gls <- getData.nls

qqnorm.gls <- qqnorm.lm <- qqnorm.nls

plot.lme <- plot.nls

fitted.gnls <- fitted.gls

residuals.gnls <- residuals.gls

vcov.gls <- function (object, ...) object$varBeta

vcov.lme <- function (object, ...) object$varFix
