% $Id: covariate<-.Rd,v 1.1 1999/10/13 00:50:11 saikat Exp $
\name{covariate<-}
\title{Assign Covariate Values}
\usage{
covariate(object) <- value
}
\alias{covariate<-}
\arguments{
 \item{object}{any object with a \code{covariate} component.}
 \item{value}{a value to be assigned to the covariate associated with
   \code{object}.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include all \code{varFunc} classes.
}
\value{
  will depend on the method function; see the appropriate documentation.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{getCovariate}}}

\examples{
## see the method function documentation
}
\keyword{models}
