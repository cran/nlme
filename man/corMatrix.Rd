% $Id: corMatrix.Rd,v 1.1 2000/03/17 22:21:22 saikat Exp $
\name{corMatrix}
\title{Extract Correlation Matrix}
\usage{
corMatrix(object, ...)
}
\alias{corMatrix}
\arguments{
 \item{object}{an object for which a correlation matrix can be
   extracted.} 
 \item{\dots}{some methods for this generic function require additional
   arguments.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include all \code{corStruct} classes.
}
\value{
  will depend on the method function used; see the appropriate
  documentation. 
}
\author{Jose Pinheiro and Douglas Bates}
\examples{
## see the method function documentation
}
\keyword{models}
