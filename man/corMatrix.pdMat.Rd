% $Id: corMatrix.pdMat.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{corMatrix.pdMat}
\title{Extract Correlation Matrix from a pdMat Object}
\usage{
\method{corMatrix}{pdMat}(object)
}
\alias{corMatrix.pdBlocked}
\alias{corMatrix.pdCompSymm}
\alias{corMatrix.pdDiag}
\alias{corMatrix.pdIdent}
\alias{corMatrix.pdMat}
\alias{corMatrix.pdSymm}

\arguments{
 \item{object}{an object inheriting from class \code{pdMat}, representing
   a positive definite matrix.}
}

\description{
  The correlation matrix corresponding to the positive-definite matrix
  represented by \code{object} is obtained. 
}
\value{
  the correlation matrix corresponding to the positive-definite matrix
  represented by \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{as.matrix.pdMat}}, \code{\link{pdMatrix}}}

\examples{
pd1 <- pdSymm(diag(1:4))
corMatrix(pd1)
}
\keyword{models}
