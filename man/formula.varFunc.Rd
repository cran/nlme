% $Id: formula.varFunc.Rd,v 1.6 2001/10/30 20:51:15 bates Exp $
\name{formula.varFunc}
\title{Extract varFunc Formula}
\usage{
\method{formula}{varFunc}(x, \dots)
}
\alias{formula.varFunc}
\alias{formula.varComb}
\arguments{
 \item{x}{an object inheriting from class \code{varFunc},
   representing a variance function structure.} 
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  This method function extracts the formula associated with a
  \code{varFunc} object, in which covariates and grouping factors are
  specified. 
}
\value{
  if \code{x} has a \code{formula} attribute, its value is
  returned; else \code{NULL} is returned.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\examples{
formula(varPower(form = ~ fitted(.) | Sex))
}
\keyword{models}
