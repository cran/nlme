% $Id: print.pdMat.Rd,v 1.4 2000/07/03 18:22:48 bates Exp $
\name{print.pdMat}
\title{Print a pdMat Object}
\usage{
print(x, ...)
}
\alias{print.pdMat}
\arguments{
 \item{object}{an object inheriting from class \code{pdMat}, representing
   a positive definite matrix.} 
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  Print the standard deviations and correlations (if any) associated the
  positive-definite matrix represented by \code{x} (considered as a
  variance-covariance matrix).
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{print.summary.pdMat}}}

\examples{
pd1 <- pdSymm(diag(1:3), nam = c("A","B","C"))
print(pd1)
}
\keyword{models}
