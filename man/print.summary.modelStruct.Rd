% $Id: print.summary.modelStruct.Rd,v 1.2 2000/04/17 19:03:40 bates Exp $
\name{print.summary.modelStruct}
%\title{Use \code{print} on a \code{summary.modelStruct} Object}
\title{Print a summary.modelStruct Object}
\usage{
print(x, ...)
}
\alias{print.summary.modelStruct}
\arguments{
 \item{x}{an object inheriting from class \code{summary.modelStruct},
   generally resulting from applying \code{summary} to an object
   inheriting from class \code{modelStruct}.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  This method function prints the constrained coefficients of
  an initialized \code{modelStruct} object, with a header specifying the
  type of correlation structure associated with the object.
}
\value{
  the printed coefficients of \code{x} in constrained form, with a
  header specifying the associated correlation structure type.
}
\author{Jose Pinheiro and Douglas Bates }

\seealso{\code{\link{summary.modelStruct}} }

\examples{
cs1 <- corAR1(0.3)
print(summary(cs1))
}
\keyword{models}
