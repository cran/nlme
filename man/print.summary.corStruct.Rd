% $Id: print.summary.corStruct.Rd,v 1.2 2000/04/17 19:03:40 bates Exp $
\name{print.summary.corStruct}
\title{Print a summary.corStruct Object}
\usage{
print(x, ...)
}
\alias{print.summary.corStruct}
\arguments{
 \item{x}{an object inheriting from class \code{summary.corStruct},
   generally resulting from applying \code{summary} to an object
   inheriting from class \code{corStruct}.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  This method function prints the constrained coefficients of
  an initialized \code{corStruct} object, with a header specifying the
  type of correlation structure associated with the object.
}
\value{
  the printed coefficients of \code{x} in constrained form, with a
  header specifying the associated correlation structure type.
}
\author{Jose Pinheiro and Douglas Bates }

\seealso{\code{\link{summary.corStruct}} }

\examples{
cs1 <- corAR1(0.3)
print(summary(cs1))
}
\keyword{models}
