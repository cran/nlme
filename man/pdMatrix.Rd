% $Id: pdMatrix.Rd,v 1.6 2001/06/18 21:16:50 bates Exp $
\name{pdMatrix}
\title{Extract Matrix or Square-Root Factor from a pdMat Object}
\usage{
pdMatrix(object, factor)
}
\alias{pdMatrix}
\alias{pdMatrix.pdBlocked}
\alias{pdMatrix.pdCompSymm}
\alias{pdMatrix.pdDiag}
\alias{pdMatrix.pdIdent}
\alias{pdMatrix.pdMat}
\alias{pdMatrix.pdSymm}
\alias{pdMatrix.pdNatural}
\arguments{
 \item{object}{an object inheriting from class \code{pdMat}, representing
   a positive definite matrix.}
 \item{factor}{an optional logical value. If \code{TRUE}, a square-root
   factor of the positive-definite matrix represented by \code{object}
   is returned; else, if \code{FALSE}, the positive-definite matrix is
   returned. Defaults to \code{FALSE}.}
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
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{as.matrix.pdMat}}, \code{\link{pdFactor}},
  \code{\link{corMatrix}}}

\examples{
pd1 <- pdSymm(diag(1:4))
pdMatrix(pd1)
}
\keyword{models}
