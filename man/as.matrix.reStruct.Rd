% $Id: as.matrix.reStruct.Rd,v 1.6 2001/10/30 20:51:14 bates Exp $
\name{as.matrix.reStruct}
\title{Matrices of an reStruct Object}
\usage{
\method{as.matrix}{reStruct}(x)
}
\alias{as.matrix.reStruct}
\arguments{
 \item{x}{an object inheriting from class \code{reStruct}, representing
   a random effects structure and consisting of a list of \code{pdMat}
   objects.}
}
\description{
  This method function extracts the positive-definite matrices
  corresponding to the \code{pdMat} elements of \code{object}. 
}
\value{
  a list with components given by the positive-definite matrices
  corresponding to the elements of \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{as.matrix.pdMat}}, \code{\link{reStruct}},
  \code{\link{pdMat}}}

\examples{
data(Orthodont)
rs1 <- reStruct(pdSymm(diag(3), ~age+Sex, data = Orthodont))
as.matrix(rs1)
}
\keyword{models}
