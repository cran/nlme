### $Id: nlsList.R,v 1.4.2.2 2003/01/18 20:49:28 saikat Exp $
###
###                  Create a list of nls objects
###
### Copyright 1997-2001  Jose C. Pinheiro <jcp@research.bell-labs.com>,
###                      Douglas M. Bates <bates@stat.wisc.edu>
###
### This file is part of the nlme library for S and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

nlsList <-
  ## A list of nls objects
  function(model, data, start, control, level, subset, na.action = na.fail,
           pool = TRUE) UseMethod("nlsList")

nlsList.selfStart <-
  function (model, data, start, control, level, subset, na.action = na.fail,
            pool = TRUE)
{
  mCall <- as.list(match.call())[-1]
  if (!inherits(data, "groupedData")) {
    stop("second argument must be a groupedData object")
  }
  marg <- substitute(model)
  if (mode(marg) != "name") {
    stop("cannot use an anonymous function for the model")
  }
					# Build up a call to the model function
  m <- call(as.character(marg))
  args <- lapply(names(formals(eval(marg))), as.name)
  args[[1]] <- getCovariateFormula(data)[[2]]
  m[1 + seq(along = args)] <- args
  form <- formula(data)
  form[[3]][[2]] <- m
  mCall$model <- form
  do.call("nlsList.formula", mCall)
}

nlsList.formula <-
  function(model, data, start = NULL, control, level, subset,
           na.action = na.fail, pool = TRUE)
{
  Call <- match.call()
  if (!missing(subset)) {
    data <-
      data[eval(asOneSidedFormula(Call[["subset"]])[[2]], data),, drop = FALSE]
  }
  if (!inherits(data, "data.frame")) data <- as.data.frame(data)
  data <- na.action(data)
  if (is.null(grpForm <- getGroupsFormula(model))) {
    if (inherits(data, "groupedData")) {
      if (missing(level)) level <- length(getGroupsFormula(data, asList = TRUE))
      else if (length(level) > 1) {
	stop("Multiple levels not allowed")
      }
      groups <- pruneLevels(getGroups(data, level = level))
      grpForm <- getGroupsFormula(data)
    } else {
      stop (paste("Data must be a groupedData object if formula",
                  "does not include groups"))
    }
  } else {
    if (missing(level)) {
      level <- length(getGroupsFormula(model, asList = TRUE))
    } else if (length(level) > 1) {
      stop("Multiple levels not allowed")
    }
    model <- eval(parse(text = paste(deparse(model[[2]]),
                        deparse(getCovariateFormula(model)[[2]]), sep = "~")))
    groups <- pruneLevels(getGroups(data, form = grpForm, level = level))
  }
  if (is.null(start) && is.null(attr(data, "parameters"))) {
    ## no starting values
    ## checking for old-style selfStart functions
    FUN <- eval(model[[3]][[1]])
    if (is.function(FUN) && class(FUN) != "selfStart" &&
        !is.null(attr(FUN, "initial"))) {
      stop(paste("old-style self-starting model functions\n",
                 "are no longer supported.\n",
                 "New selfStart functions are available.\n",
                 "Use\n",
                 "  SSfpl instead of fpl,\n",
                 "  SSfol instead of first.order.log,\n",
                 "  SSbiexp instead of biexp,\n",
                 "  SSlogis instead of logistic.\n",
                 "If writing your own selfStart model, see\n",
                 "  \"help(selfStart)\"\n",
                 "for the new form of the \"initial\" attribute.", sep="" ))
    }
  }

  controlvals <- nls.control()
  if(!missing(control)) {
    controlvals[names(control)] <- control
  }
  val <- lapply(split(data, groups),
		function(dat, formula, start, control, first = TRUE)
		{
                  ans <- try({
                    data <- as.data.frame(dat)
                    if (is.null(start)) {
                      nls(formula = formula, data = data, control = control)
                    } else {
                      nls(formula = formula, data = data, start = start,
                          control = control)
                    }
                  })
                  if (inherits(ans, "try-error"))
                    NULL
                  else ans
		}, formula = model, start = start, control = controlvals)
  if (inherits(data, "groupedData")) {
    ## saving labels and units for plots
    attr(val, "units") <- attr(data, "units")
    attr(val, "labels") <- attr(data, "labels")
    attr(val, "outer") <- attr(data, "outer")
  }
  attr(val, "dims") <- list(N = nrow(data), M = length(val))
  attr(val, "call") <- Call
  attr(val,"groups") <- ordered(groups, levels = names(val))
  attr(val, "origOrder") <- match(unique(as.character(groups)), names(val))
  attr(val, "pool") <- pool
  attr(val, "groupsForm") <- grpForm
  class(val) <- c("nlsList", "lmList")
  val
}

###*# Methods for standard generics

coef.summary.nlsList <-
  function(object, ...) object$parameters

formula.nlsList <-
  function(x, ...) eval(attr(x, "call")[["model"]])

summary.nlsList <-
  function(object, ...)
{
  val <- NextMethod("summary")
  class(val) <- c("summary.nlsList", class(val))
  val
}

update.nlsList <-
  function(object, model, data, start, control, level, subset, na.action,
	   pool, ...)
{
  thisCall <- as.list(match.call())[-(1:2)]
  if (!missing(model)) {
    names(thisCall)[match(names(thisCall), "model")] <- "object"
  }
  nextCall <- as.list(attr(object, "call")[-1])
  nextCall[names(thisCall)] <- thisCall
  do.call("nlsList", nextCall)
}

### Local variables:
### mode: S
### End:
