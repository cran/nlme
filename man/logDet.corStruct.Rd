% $Id: logDet.corStruct.Rd,v 1.4 2000/07/03 18:22:47 bates Exp $
\name{logDet.corStruct}
\title{Extract corStruct Log-Determinant}
\usage{
logDet(object, covariate)
}
\alias{logDet.corStruct}
\alias{logDet.corIdent}
\arguments{
 \item{object}{an object inheriting from class \code{corStruct},
   representing a correlation structure.}
 \item{covariate}{an optional covariate vector (matrix), or list of
   covariate vectors (matrices), at which values the correlation matrix,
   or list of correlation  matrices, are to be evaluated. Defaults to
   \code{getCovariate(object)}.} 
}
\description{
  This method function extracts the logarithm of the determinant of a
  square-root factor of the correlation matrix associated with
  \code{object}, or the sum of the log-determinants of square-root
  factors of the list of correlation matrices associated with
  \code{object}.  
}
\value{
  the log-determinant of a square-root factor of the correlation matrix
  associated with \code{object}, or the sum of the log-determinants of
  square-root factors of the list of correlation matrices associated
  with \code{object}. 
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{logLik.corStruct}}, \code{\link{corMatrix.corStruct}}}

\examples{
cs1 <- corAR1(0.3)
logDet(cs1, covariate = 1:4)
}
\keyword{models}
