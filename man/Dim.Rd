% $Id: Dim.Rd,v 1.1 1999/10/13 00:50:10 saikat Exp $
\name{Dim}
\title{Extract Dimensions from an Object}
\usage{
Dim(object, ...)
}
\alias{Dim}
\arguments{
 \item{object}{any object for which dimensions can be extracted.}
 \item{\dots}{some methods for this generic function require additional
   arguments.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include: \code{corSpatial}, \code{corStruct},
  \code{pdCompSymm}, \code{pdDiag}, \code{pdIdent}, \code{pdMat},
  and \code{pdSymm}.
}
\value{
  will depend on the method function used; see the appropriate documentation.
}
\author{Jose Pinheiro and Douglas Bates}
\note{
  If \code{dim} allowed more than one argument, there would be no
  need for this generic function.
}
\seealso{
  \code{\link{Dim.pdMat}}, \code{\link{Dim.corStruct}}
}
\examples{
## see the method function documentation
}
\keyword{models}
