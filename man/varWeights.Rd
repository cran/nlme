% $Id: varWeights.Rd,v 1.2 1999/11/04 16:39:00 saikat Exp $
\name{varWeights}
\title{Extract Variance Function Weights}
\usage{
varWeights(object)
}
\alias{varWeights}
\alias{varWeights.varComb}
\alias{varWeights.varFunc}
\arguments{
 \item{object}{an object inheriting from class \code{varFunc},
   representing a variance function structure.} 
}
\description{
The inverse of the standard deviations corresponding to the variance
function structure represented by \code{object} are returned.
}
\value{
  if \code{object} has a \code{weights} attribute, its value is
  returned; else \code{NULL} is returned.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{logLik.varFunc}}}

\examples{
library(nlme)
data(Orthodont)
vf1 <- varPower(form=~age)
vf1 <- initialize(vf1, Orthodont)
coef(vf1) <- 0.3
varWeights(vf1)[1:10]
}
\keyword{models}
