% $Id: corMatrix.reStruct.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{corMatrix.reStruct}
\title{Extract Correlation Matrix from Components of an reStruct Object}
\usage{
corMatrix(object)
}
\alias{corMatrix.reStruct}

\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
}
\description{
  This method function extracts the correlation matrices
  corresponding to the \code{pdMat} elements of \code{object}. 
}
\value{
  a list with components given by the correlation matrices
  corresponding to the elements of \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{as.matrix.reStruct}}, \code{\link{reStruct}},
  \code{\link{pdMat}}}

\examples{
data(Orthodont)
rs1 <- reStruct(pdSymm(diag(3), ~age+Sex, data = Orthodont))
corMatrix(rs1)
}
\keyword{models}
