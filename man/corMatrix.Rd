% File nlme/man/corMatrix.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{corMatrix}
\title{Extract Correlation Matrix}
\usage{
corMatrix(object, \dots)
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
  this function include all \code{"corStruct"} classes,
  see \sQuote{\link{corClasses}}.
}
\value{
  will depend on the method function used; see the appropriate
  documentation. 
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{
  \code{\link{corMatrix.corStruct}},
  \code{\link{corMatrix.pdMat}},
  \code{\link{corMatrix.reStruct}}
}
\keyword{models}
