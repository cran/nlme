% $Id: corMatrix.reStruct.Rd,v 1.6 2001/10/30 20:51:15 bates Exp $
\name{corMatrix.reStruct}
\title{Extract Correlation Matrix from Components of an reStruct Object}
\usage{
\method{corMatrix}{reStruct}(object, \dots)
}
\alias{corMatrix.reStruct}

\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  This method function extracts the correlation matrices
  corresponding to the \code{pdMat} elements of \code{object}. 
}
\value{
  a list with components given by the correlation matrices
  corresponding to the elements of \code{object}.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{as.matrix.reStruct}}, \code{\link{reStruct}},
  \code{\link{pdMat}}}

\examples{
data(Orthodont)
rs1 <- reStruct(pdSymm(diag(3), ~age+Sex, data = Orthodont))
corMatrix(rs1)
}
\keyword{models}
