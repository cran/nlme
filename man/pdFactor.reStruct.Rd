% $Id: pdFactor.reStruct.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{pdFactor.reStruct}
\title{Extract Square-Root Factor from Components of an reStruct Object}
\usage{
pdFactor(object)
}
\alias{pdFactor.reStruct}

\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
}
\description{
  This method function extracts square-root factors of the
  positive-definite matrices corresponding to the \code{pdMat} elements
  of \code{object}.  
}
\value{
  a vector with square-root factors of the positive-definite matrices
  corresponding to the elements of \code{object} stacked column-wise.
}
\author{Jose Pinheiro and Douglas Bates}
\note{
  This function is used intensively in optimization
  algorithms and its value is returned as a vector for efficiency
  reasons. The \code{pdMatrix} function can be used to obtain
  square-root factors in matrix form.
}
\seealso{\code{\link{pdMatrix.reStruct}}, \code{\link{pdFactor.pdMat}}}
\examples{
library(nlme)
data(Orthodont)
rs1 <- reStruct(pdSymm(diag(3), ~age+Sex, data = Orthodont))
pdFactor(rs1)
}
\keyword{models}
