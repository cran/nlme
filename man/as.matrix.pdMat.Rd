% $Id: as.matrix.pdMat.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{as.matrix.pdMat}
\title{Matrix of a pdMat Object}
\usage{
as.matrix(x)
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
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{pdMat}}, \code{\link{corMatrix}}}

\examples{
library(nlme)
as.matrix(pdSymm(diag(4)))
}

\keyword{models}
