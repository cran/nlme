% $Id: logDet.reStruct.Rd,v 1.4 2000/07/03 18:22:47 bates Exp $
\name{logDet.reStruct}
\title{Extract reStruct Log-Determinants}
\usage{
logDet(object)
}
\alias{logDet.reStruct}
\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
}
\description{
  Calculates, for each of the \code{pdMat} components of \code{object},
  the logarithm of the determinant of a square-root factor.  
}
\value{
  a vector with the log-determinants of square-root factors of the
  \code{pdMat} components of \code{object}.
}

\author{Jose Pinheiro}

\seealso{\code{\link{reStruct}}, \code{\link{pdMat}}}

\examples{
rs1 <- reStruct(list(A = pdSymm(diag(1:3), form = ~Score),
  B = pdDiag(2 * diag(4), form = ~Educ)))
logDet(rs1)
}
\keyword{models}
