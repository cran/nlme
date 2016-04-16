% File nlme/man/Matrix.pdMat.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{Matrix.pdMat}
\title{Assign Matrix to a pdMat or pdBlocked Object}
\alias{matrix<-.pdMat}
\alias{matrix<-.pdBlocked}
\usage{
\method{matrix}{pdMat}(object) <- value
\method{matrix}{pdBlocked}(object) <- value
}
\arguments{
 \item{object}{an object inheriting from class \code{"\link{pdMat}"},
    representing a positive definite matrix.}
 \item{value}{a matrix with the new values to be assigned to the
   positive-definite matrix represented by \code{object}.  Must have the
   same dimensions as \code{as.matrix(object)}.}
}
\description{
  The positive-definite matrix represented by \code{object} is replaced
  by \code{value}.  If the original matrix had row and/or column names,
  the corresponding names for \code{value} can either be \code{NULL}, or
  a permutation of the original names.
}
\value{
  a \code{pdMat} or \code{pdBlocked} object similar to \code{object},
  but with its coefficients modified to produce the matrix in
  \code{value}.
}

\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{pdMat}}, \code{"\link{matrix<-}"}
}

\examples{
class(pd1 <- pdSymm(diag(3))) # "pdSymm" "pdMat"
matrix(pd1) <- diag(1:3)
pd1
}
\keyword{models}
