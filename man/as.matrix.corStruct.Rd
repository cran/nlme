% $Id: as.matrix.corStruct.Rd,v 1.2 2000/04/17 19:03:38 bates Exp $
\name{as.matrix.corStruct}
\title{Matrix of a corStruct Object}
\usage{
as.matrix(x)
}
\alias{as.matrix.corStruct}
\arguments{
 \item{x}{an object inheriting from class \code{corStruct},
   representing a correlation structure.}
}
\description{
  This method function extracts the correlation matrix, or list of
  correlation matrices, associated with \code{object}.
}
\value{
  If the correlation structure includes a grouping factor, the returned
  value will be a list with components given by the correlation
  matrices for each group. Otherwise, the returned value will be a
  matrix representing the correlation structure associated with
  \code{object}. 
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{corClasses}}, \code{\link{corMatrix}}}
\examples{
data(Orthodont)
cst1 <- corAR1(form = ~1|Subject)
cst1 <- initialize(cst1, data = Orthodont)
as.matrix(cst1)
}
\keyword{models}
