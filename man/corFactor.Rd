% $Id: corFactor.Rd,v 1.1 1999/10/13 00:50:11 saikat Exp $
\name{corFactor}
\title{Factor of a Correlation Matrix}
\usage{
corFactor(object, ...)
}
\alias{corFactor}
\arguments{
 \item{object}{an object from which a correlation matrix can be
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
