% File nlme/man/Coef.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{Coef}
\title{Assign Values to Coefficients}
\usage{
coef(object, \dots) <-  value
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
  this function include all \code{"\link{pdMat}"},
  \code{"\link{corStruct}"} and \code{"\link{varFunc}"} classes,
  \code{"\link{reStruct}"}, and \code{"modelStruct"}.
}
\value{
  will depend on the method function; see the appropriate documentation.
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{coef}}}

\examples{
## see the method function documentation
}
\keyword{models}
