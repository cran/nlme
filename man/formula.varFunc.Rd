% $Id: formula.varFunc.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{formula.varFunc}
\title{Extract varFunc Formula}
\usage{
\method{formula}{varFunc}(object)
}
\alias{formula.varFunc}
\alias{formula.varComb}
\arguments{
 \item{object}{an object inheriting from class \code{varFunc},
   representing a variance function structure.} 
}
\description{
  This method function extracts the formula associated with a
  \code{varFunc} object, in which covariates and grouping factors are
  specified. 
}
\value{
  if \code{object} has a \code{formula} attribute, its value is
  returned; else \code{NULL} is returned.
}

\author{Jose Pinheiro and Douglas Bates}

\examples{
formula(varPower(form = ~ fitted(.) | Sex))
}
\keyword{models}
