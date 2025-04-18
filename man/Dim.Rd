% File nlme/man/Dim.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{Dim}
\title{Extract Dimensions from an Object}
\usage{
Dim(object, \dots)
}
\alias{Dim}
\alias{Dim.default}
\arguments{
 \item{object}{any object for which dimensions can be extracted.}
 \item{\dots}{some methods for this generic function require additional
   arguments.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include: \code{"\link{corSpatial}"}, \code{"\link{corStruct}"},
  \code{"pdCompSymm"}, \code{"pdDiag"}, \code{"pdIdent"}, \code{"pdMat"},
  and \code{"pdSymm"}.
}
\value{
  will depend on the method function used; see the appropriate documentation.
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}
\note{
  If \code{dim} allowed more than one argument, there would be no
  need for this generic function.
}
\seealso{
  \code{\link{Dim.corSpatial}},
  \code{\link{Dim.pdMat}}, \code{\link{Dim.corStruct}}
}
\keyword{models}
