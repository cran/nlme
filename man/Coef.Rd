% $Id: Coef.Rd,v 1.3 2000/04/05 18:43:21 bates Exp $
\name{Coef}
\title{Assign Values to Coefficients}
\usage{
coef(object, ...) <-  value
}
\alias{coef<-}
\alias{coefficients<-}
\arguments{
 \item{object}{any object representing a fitted model, or, by default,
   any object with a \code{coef} component.}
 \item{\dots}{some methods for this generic function may require
   additional arguments.}
 \item{value}{a value to be assigned to the coefficients associated with
 \code{object}.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include all \code{pdMat}, \code{corStruct}, and 
  \code{varFunc} classes, \code{reStruct}, and \code{modelStruct}.
}
\value{
  will depend on the method function; see the appropriate documentation.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{coef}}}

\examples{
## see the method function documentation
}
\keyword{models}
