% $Id: Dim.pdMat.Rd,v 1.2 1999/11/04 16:38:58 saikat Exp $
\name{Dim.pdMat}
\title{Dimensions of a pdMat Object}
\usage{
Dim(object)
}
\alias{Dim.pdMat}
\alias{Dim.pdCompSymm}
\alias{Dim.pdDiag}
\alias{Dim.pdIdent}
\arguments{
 \item{object}{an object inheriting from class \code{pdMat},
   representing a positive-definite matrix.}
}
\description{
  This method function returns the dimensions of the matrix represented
  by \code{object}. 
}
\value{
  an integer vector with the number of rows and columns of the
  matrix represented by \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{Dim}}}
}
\examples{
library(nlme)
Dim(pdSymm(diag(3)))
}
\keyword{models}
