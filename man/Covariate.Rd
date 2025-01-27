% File nlme/man/Covariate.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{Covariate}
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
  this function include all \code{"\link{varFunc}"} classes,
  see \sQuote{\link{varClasses}}.
}
\value{
  will depend on the method function; see the appropriate documentation.
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{covariate<-.varFunc}}, \code{\link{getCovariate}}}

\keyword{models}
