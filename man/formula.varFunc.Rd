% $Id: formula.varFunc.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{formula.varFunc}
\title{Extract varFunc Formula}
\usage{
formula(object)
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
library(nlme)
formula(varPower(form = ~ fitted(.) | Sex))
}
\keyword{models}
