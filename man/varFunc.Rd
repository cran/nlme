% $Id: varFunc.Rd,v 1.2 1999/11/04 16:39:00 saikat Exp $
\name{varFunc}
\title{Variance Function Structure}
\usage{
varFunc(object)
}
\alias{varFunc}
\arguments{
 \item{object}{either an one-sided formula specifying a variance
   covariate, or an object inheriting from class \code{varFunc},
   representing a variance function structure.} 
}
\description{
  If \code{object} is a one-sided formula, it is used as the argument to
  \code{varFixed} and the resulting object is returned. Else, if
  \code{object} inherits from class \code{varFunc}, it is returned
  unchanged. 
}
\value{
  an object from class \code{varFunc}, representing a variance function
  structure. 
}

\author{Jose Pinheiro and Douglas Bates}


\seealso{\code{\link{varFixed}}, \code{\link{varWeights.varFunc}},
  \code{\link{coef.varFunc}}} 

\examples{
library(nlme)
vf1 <- varFunc(~age)
}
\keyword{models}
