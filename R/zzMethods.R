###   Miscellaneous methods that must be defined last in the library
###
### Copyright 1997-2003  Jose C. Pinheiro,
###                      Douglas M. Bates <bates@stat.wisc.edu>
# Copyright 2007-2016 The R Core team
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

## Because of the conditional pre-/post 3.3.0 behavior with 'sigma' :
.onLoad <- function(libname, pkgname) {
    ## pInfo <- readRDS(attr(packageDescription(pkgname), "file"))
    ## built.R.ver <- pInfo$Built$R
    ## First case: does signal Error : object 'sigma' is not exported by 'namespace:stats'
    ##             so we don't need the following:
    ## if(getRversion() < "3.3") { ## check that nlme was not *installed* with R >= 3.3.x
    ##     if(built.R.ver >= "3.3")
    ##         warning("Package ", dQuote(pkgname), " installed with R version ", built.R.ver,
    ##     	 " can not safely be used with old R version ", getRversion())
    ##
    ## second case:
    if(getRversion() >= "3.3") {
        pInfo <- readRDS(attr(packageDescription(pkgname), "file"))
        built.R.ver <- pInfo$Built$R
	if(built.R.ver < "3.3")## installed with R < 3.3.x :
	    warning("Package ", dQuote(pkgname), " installed with old R version ",
		    built.R.ver, " should not be used with R version ", getRversion(),
		    "\n  Rather re-install it with this version of R.")
    }
}

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

deviance.gls <- deviance.lme <- function(object, ...)
    stop("not yet implemented.  Contributions are welcome (should be compatible with lme4's)")

if(getRversion() < "3.3") {
    sigma <- function(object, ...) UseMethod("sigma")
}

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

##  at the very end : ---------------------------
.onUnload <- function(libpath)
    library.dynam.unload("nlme", libpath)

