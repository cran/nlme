% $Id: nlsList.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{nlsList}
\title{List of nls Objects with a Common Model}
\usage{
nlsList(model, data, start, control, level, na.action, pool)
}
\alias{nlsList}
\alias{nlsList.formula}
\arguments{
  \item{model}{either a nonlinear model formula, with the response on
    the left of a \code{~} operator and an expression involving
    parameters, covariates, and a grouping factor separated by the
    \code{|} operator on the right, or a \code{selfStart} function.  The
    method function \code{nlsList.selfStart} is documented separately.
  }
  \item{data}{
    a data frame in which to interpret the variables named in
    \code{model}.
  }
  \item{start}{an optional named list with initial values for the
    parameters to be estimated in \code{model}. It is passed as the
    \code{start} argument to each \code{nls} call and is required when
    the nonlinear function in \code{model} does not inherit from class
    \code{selfStart}.
  }
  \item{control}{a list of control values passed as the \code{control}
    argument to \code{nls}. Defaults to an empty list.
  }
  \item{level}{
    an optional integer specifying the level of grouping to be used when 
    multiple nested levels of grouping are present.
  }
  \item{na.action}{a function that indicates what should happen when the
    data contain \code{NA}s.  The default action (\code{na.fail}) causes
    \code{nlsList} to print an error message and terminate if there are any
    incomplete observations.
  }
  \item{pool}{
    an optional logical value that is preserved as an attribute of the
    returned value.  This will be used as the default for \code{pool} in
    calculations of standard deviations or standard errors for summaries.
  }
}
\description{
  \code{Data} is partitioned according to the levels of the grouping
  factor defined in \code{model} and individual \code{nls} fits are
  obtained for each \code{data} partition, using the model defined in
  \code{model}.
}
\value{
  a list of \code{nls} objects with as many components as the number of
  groups defined by the grouping factor. Generic functions such as
  \code{coef}, \code{fixed.effects}, \code{lme}, \code{pairs},
  \code{plot}, \code{predict}, \code{random.effects}, \code{summary},
  and \code{update} have methods that can be applied to an \code{nlsList}
  object. 
}
\seealso{
\code{\link{nls}}, \code{\link{nlme.nlsList}}.
}
\examples{
data(CO2)
fm1 <- nlsList(uptake ~ SSasympOff(conc, Asym, lrc, c0),
   data = CO2, start = c(Asym = 30, lrc = -4.5, c0 = 52))
fm1
}
\keyword{models}

