% $Id: logDet.Rd,v 1.3 2000/07/03 18:22:47 bates Exp $
\name{logDet}
\title{Extract the Logarithm of the Determinant}
\usage{
logDet(object, ...)
}
\alias{logDet}
\arguments{
 \item{object}{any object from which a matrix, or list of matrices, can
   be extracted}
 \item{\dots}{some methods for this generic function require additional
   arguments.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include: \code{corStruct}, several \code{pdMat} classes,
  and \code{reStruct}. 
}
\value{
  will depend on the method function used; see the appropriate
  documentation. 
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{logLik}}}
\examples{
## see the method function documentation
}
\keyword{models}
