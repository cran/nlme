% $Id: print.modelStruct.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{print.modelStruct}
%\title{Use \code{print} on a \code{modelStruct} Object}
\title{Print a modelStruct Object}
\usage{
print(x, ...)
}
\alias{print.modelStruct}
\arguments{
 \item{x}{an object inheriting from class \code{modelStruct},
   representing a list of model components, such as \code{corStruct} and
   \code{varFunc} objects.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  This method function applies \code{print} to each element of
  \code{object}.
}
\value{
  the printed elements of \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{print}}}

\examples{
library(nlme)
lms1 <- lmeStruct(reStruct = reStruct(pdDiag(diag(2), ~age)),
   corStruct = corAR1(0.3))
print(lms1)
}
\keyword{models}
