% $Id: solve.pdMat.Rd,v 1.4 2000/07/03 18:22:48 bates Exp $
\name{solve.pdMat}
\title{Calculate Inverse of a Positive-Definite Matrix}
\usage{
solve(a, b, tol)
}
\alias{solve.pdMat}
\alias{solve.pdBlocked}
\alias{solve.pdDiag}
\alias{solve.pdIdent}
\alias{solve.pdNatural}

\arguments{
 \item{a}{an object inheriting from class \code{pdMat}, representing
   a positive definite matrix.} 
 \item{b}{this argument is only included for consistency with the generic
   function and is not used in this method function.}
 \item{tol}{an optional numeric value for the tolerance used in the
   numerical algorithm. Defaults to \code{1e-7}.}
}
\description{
  The positive-definite matrix represented by \code{a} is inverted and
  assigned to \code{a}.
}
\value{
  a \code{pdMat} object similar to \code{a}, but with coefficients
  corresponding to the inverse of the positive-definite matrix
  represented by \code{a}.
}

\author{Jose Pinheiro and Douglas Bates}


\seealso{\code{\link{pdMat}}}

\examples{
pd1 <- pdCompSymm(3 * diag(3) + 1)
solve(pd1)
}
\keyword{models}
