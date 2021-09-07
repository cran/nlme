###   Miscellaneous methods that must be defined last in the library
###
### Copyright 2007-2018 The R Core team
### Copyright 1997-2003  Jose C. Pinheiro,
###                      Douglas M. Bates <bates@stat.wisc.edu>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

`%||%` <- function(x, y) if(is.null(x)) y else x

## Note that  require( stats )  has already happened ...

comparePred.lme <- comparePred.lmList <- comparePred.gls

getData.nlme <- getData.gnls

getData.lme <- getData.gls <- getData.nls

qqnorm.gls <- qqnorm.lm <- qqnorm.nls

plot.lme <- plot.nls

fitted.gnls <- fitted.gls

residuals.gnls <- residuals.gls

vcov.gls <- function (object, ...) object$varBeta

vcov.lme <- function (object, ...) object$varFix

deviance.gls <- deviance.lme <- function(object, ...) {
    if(object$method == "ML")
	-2 * logLik(object)
    else {
	warning("deviance undefined for REML fit")
	NULL
    }
}

## From MASS/R/stepAIC.R :
extractAIC.gls <- extractAIC.lme <- function(fit, scale, k = 2, ...)
{
    if(fit$method != "ML") stop("AIC undefined for REML fit")
    res <- logLik(fit)
    edf <- attr(res, "df")
    c(edf,  -2*res + k * edf)
}

terms.gls <- function(x, ...) terms(formula(x), ...)
if(FALSE)## Not needed, because 'lme' object has "terms" attribute:
    terms.lme <- function(x, ...) terms(formula(x), ...)
## end{from MASS}


sigma.gls <- sigma.lme <- function(object, ...) object$sigma

## also works for "nlsList"
sigma.lmList <- function(object, ...) vapply(object, sigma, 1, ...)

## confint() works for "gls" via confint.default() !
confint.lme <- function(object, ...)
    stop("not (yet) implemented.  Contributions are welcome; use intervals() instead (for now)")

confint.lmList <- function(object, ...) sapply(object, confint, ..., simplify=FALSE)
confint.nlsList <- function(object, ...) {
    sapply(object, function(ob) tryCatch(confint(ob, ...), error = function(e)
	structure(c(NA,NA), errMsg = conditionMessage(e))),
	simplify=FALSE)
}

.ns <- environment() # == asNamespace("nlme")

##  at the very end : ---------------------------
.onUnload <- function(libpath)
    library.dynam.unload("nlme", libpath)
