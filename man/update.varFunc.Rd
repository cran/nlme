% $Id: update.varFunc.Rd,v 1.1 2000/03/17 22:21:24 saikat Exp $
\name{update.varFunc}
\title{Update varFunc Object}
\usage{
update(object, data)
}
\alias{update.varFunc}
\alias{update.varComb}
\alias{update.varConstPower}
\alias{update.varExpon}
\alias{update.varPower}

\arguments{
 \item{object}{an object inheriting from class \code{varFunc},
   representing a variance function structure.} 
 \item{data}{a list with a component named \code{"."} with the current
   version of the fitted object (from which fitted values, coefficients,
   and residuals can be extracted) and, if necessary, other variables
   used to evaluate the variance covariate(s).}
}
\description{
  If the \code{formula(object)} includes a \code{"."} term, representing
  a fitted object, the variance covariate needs to be updated upon
  completion of an optimization cycle (in which the variance function
  weights are kept fixed). This method function allows a reevaluation of
  the variance covariate using the current fitted object and,
  optionally, other variables in the original data.
}
\value{
  if \code{formula(object)} includes a \code{"."} term, an
  \code{varFunc} object similar to \code{object}, but with the 
  variance covariate reevaluated at the current fitted object value;
  else \code{object} is returned unchanged.
}

\author{Jose Pinheiro and Douglas Bates}


\seealso{\code{\link{needUpdate}}, \code{\link{covariate<-.varFunc}}}

\keyword{models}
