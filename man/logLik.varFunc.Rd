% $Id: logLik.varFunc.Rd,v 1.6 2001/10/30 20:51:15 bates Exp $
\name{logLik.varFunc}
\title{Extract varFunc logLik}
\usage{
\method{logLik}{varFunc}(object, data, \dots)
}
\alias{logLik.varFunc}
\alias{logLik.varComb}
\arguments{
 \item{object}{an object inheriting from class \code{varFunc},
   representing a variance function structure.} 
 \item{data}{this argument is included to make this method function
   compatible with other \code{logLik} methods and will be ignored.} 
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  This method function extracts the component of a Gaussian
  log-likelihood associated with the variance function structure
  represented by \code{object}, which is equal to the sum of the
  logarithms of the corresponding weights. 
}
\value{
  the sum of the logarithms of the weights corresponding to the variance
  function structure represented by \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}


\examples{
data(Orthodont)
vf1 <- varPower(form = ~age)
vf1 <- initialize(vf1, Orthodont)
coef(vf1) <- 0.1
logLik(vf1)
}
\keyword{models}
