% $Id: logLik.corStruct.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{logLik.corStruct}
\title{Extract corStruct Log-Likelihood}
\usage{
logLik(object, data)
}
\alias{logLik.corStruct}
\arguments{
 \item{object}{an object inheriting from class \code{corStruct},
   representing a correlation structure.}
 \item{data}{this argument is included to make this method function
   compatible with other \code{logLik} methods and will be ignored.} 
}
\description{
  This method function extracts the component of a Gaussian
  log-likelihood associated with the correlation structure, which is
  equal to the negative of the logarithm of the determinant (or sum of
  the logarithms of the determinants) of the matrix (or matrices)
  represented by \code{object}.
}
\value{
  the negative of the logarithm of the determinant (or sum of
  the logarithms of the determinants) of the correlation matrix (or
  matrices) represented by \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{logDet.corStruct}} }

\examples{
library(nlme)
data(Orthodont)
cs1 <- corAR1(0.2)
cs1 <- initialize(cs1, data = Orthodont)
logLik(cs1)
}
\keyword{models}
