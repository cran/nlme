% $Id: as.matrix.pdMat.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{as.matrix.pdMat}
\title{Matrix of a pdMat Object}
\usage{
\method{as.matrix}{pdMat}(x)
}
\alias{as.matrix.pdMat}
\arguments{
 \item{x}{an object inheriting from class \code{pdMat}, representing a
   positive-definite matrix.}
}
\description{
  This method function extracts the positive-definite matrix represented
  by \code{x}.
}
\value{
  a matrix corresponding to the positive-definite matrix represented by
  \code{x}.
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{pdMat}}, \code{\link{corMatrix}}}

\examples{
as.matrix(pdSymm(diag(4)))
}

\keyword{models}
