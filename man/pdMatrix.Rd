% $Id: pdMatrix.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{pdMatrix}
\title{Extract Matrix or Square-Root Factor from a pdMat Object}
\usage{
pdMatrix(object, fact, ...)
}
\alias{pdMatrix}
\alias{pdMatrix.pdBlocked}
\alias{pdMatrix.pdCompSymm}
\alias{pdMatrix.pdDiag}
\alias{pdMatrix.pdIdent}
\alias{pdMatrix.pdMat}
\alias{pdMatrix.pdSymm}

\arguments{
 \item{object}{an object inheriting from class \code{pdMat}, representing
   a positive definite matrix.}
 \item{fact}{an optional logical value. If \code{TRUE}, a square-root
   factor of the positive-definite matrix represented by \code{object}
   is returned; else, if \code{FALSE}, the positive-definite matrix is
   returned. Defaults to \code{FALSE}.}
 \item{\dots}{some methods for the generic function may require
   additional arguments.}
}

\description{
  The positive-definite matrix represented by \code{object}, or a
  square-root factor of it is obtained. Letting \eqn{\Sigma}{S} denote a
  positive-definite matrix, a square-root factor of \eqn{\Sigma}{S} is
  any square matrix \eqn{L}{L} such that \eqn{\Sigma = L'L}{S =
    L'L}. This function extracts \eqn{S} or \eqn{L}.
}
\value{
  if \code{fact} is \code{FALSE} the positive-definite matrix
  represented by \code{object} is returned; else a square-root of the
  positive-definite matrix is returned.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{as.matrix.pdMat}}, \code{\link{pdFactor}},
  \code{\link{corMatrix}}}

\examples{
library(nlme)
pd1 <- pdSymm(diag(1:4))
pdMatrix(pd1)
}
\keyword{models}
