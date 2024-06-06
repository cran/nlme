## Contributed by Mary Lindstrom <lindstro@biostat.wisc.edu>

# Copyright 2007-2024 The R Core team
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

getVarCov <- function(obj, ...) UseMethod("getVarCov")

getVarCov.lme <-
    function(obj,
             individuals,
             type= c("random.effects","conditional","marginal"), ...)
{
    type  <-  match.arg(type)
    if(inherits(obj, "nlme"))
        stop("not implemented for \"nlme\" objects")
    if(length(obj$groups) > 1)
        stop("not implemented for multiple levels of nesting")
    sigma <- obj$sigma
    D <- as.matrix(obj$modelStruct$reStruct[[1]]) * sigma^2
    if (type=="random.effects")
    {
        result  <-  D
    }
    else
    {
        result <- list()
        groups  <-  obj$groups[[1]]
        ugroups  <-  unique(groups)
        ## CAVE: both default and numeric 'individuals' are undocumented,
        ##       but unlike getVarCov.gls and perhaps unintentionally, these
        ##       index the groups as they occur in the data, not the levels!
        if (missing(individuals)) {
            individuals  <-  ugroups[1]
        } else if (is.numeric(individuals)) {
            individuals  <-  ugroups[individuals]
        }
        for (individ in as.character(individuals))
        {
            ind <- groups == individ
            ni <- sum(ind, na.rm = TRUE)
            if (ni == 0)
                stop(gettextf("individual %s was not used in the fit",
                              sQuote(individ)), domain = NA)
            if(!is.null(csT <- obj$modelStruct$corStruct)
               && ni > 1) { # corMatrix.corSpatial() excludes 1-obs groups (PR#16806)
                V <- corMatrix(csT)[[individ]]
            } else
                V <- diag(ni)
            if(!is.null(obj$modelStruct$varStruct)) {
                ## CAVE: stored weights are based on internally reordered data,
                ##       so cannot be indexed via obj$groups
                grp <- if(!is.null(csT)) getGroups(csT) else groups[order(groups)]
                sds <- 1/varWeights(obj$modelStruct$varStruct)[grp == individ]
            } else
                sds <- rep(1, ni)
            sds <- obj$sigma * sds
            cond.var <- t(V * sds) * sds
            dimnames(cond.var)  <-  list(1:nrow(cond.var),1:ncol(cond.var))
            if (type=="conditional")
                result[[individ]] <- cond.var
            else
            {
                Z <- model.matrix(obj$modelStruct$reStruct,
                                  getData(obj))[ind, , drop = FALSE]
                result[[individ]] <-
                    cond.var + Z %*% D %*% t(Z)
            }
        }
    }
    class(result)  <-  c(type,"VarCov")
    attr(result,"group.levels")  <-  names(obj$groups)
    result
}

getVarCov.gls <-
    function(obj, individual = 1, ...)
{
    if (is.null(csT <- obj$modelStruct$corStruct))
        stop("not implemented for uncorrelated errors")
    if (is.null(grp <- getGroups(csT)))
        stop("not implemented for correlation structures without a grouping factor")
    ind <- if (is.numeric(individual)) {
               as.integer(grp) == individual
           } else         grp  == individual
    ni <- sum(ind, na.rm = TRUE)
    if (ni == 0)
        stop(gettextf("individual %s was not used in the fit",
                      sQuote(individual)), domain = NA)
    S <- if (ni > 1) corMatrix(csT)[[individual]] else diag(1) # PR#16806
    if (!is.null( obj$modelStruct$varStruct))
    {
        vw  <-  1/varWeights(obj$modelStruct$varStruct)[ind]
    }
    else vw  <-  rep(1, ni)
    vars  <-  (obj$sigma * vw)^2
    result  <-  t(S * sqrt(vars))*sqrt(vars)
    class(result)  <-  c("marginal","VarCov")
    attr(result,"group.levels")  <-  names(obj$groups)
    result
}

print.VarCov <-
    function(x, corr = FALSE, stdevs = TRUE, digits = 5, ...)
{
    pvc  <-  function(x, type, corr, stdevs, digits) {
        cat(c("Random effects","Conditional",
              "Marginal")[match(type,
                                c("random.effects","conditional",
                                  "marginal"))], " ", sep = "")
        x  <-  as.matrix(x)
        class(x)  <-  NULL
        attr(x,"group.levels")  <-  NULL
        if (corr)
        {
            cat("correlation matrix\n")
            sds <- sqrt(diag(x))
            print(signif(t(x/sds)/sds,digits))
        }
        else
        {
            cat("variance covariance matrix\n")
            print(signif(x,digits))
            if(stdevs)
                sds <- sqrt(diag(x))
        }
        if (stdevs) cat("  Standard Deviations:",signif(sds,digits),"\n")
    }
    if (!is.list(x))
        pvc(x,class(x)[1],corr,stdevs,digits)
    else
    {
        for (nm in names(x))
        {
            cat(attr(x,"group.levels"),nm,"\n")
            pvc(x[[nm]],class(x)[1],corr,stdevs,digits)
        }
    }
    invisible(x)
}

