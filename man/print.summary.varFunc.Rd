% $Id: print.summary.varFunc.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{print.summary.varFunc}
\title{Print a summary.varFunc Object}
\usage{
print(x, header, ...)
}
\alias{print.summary.varFunc}
\alias{print.summary.varComb}
\alias{print.summary.varFixed}

\arguments{
 \item{x}{an object inheriting from class \code{varFunc}, representing a
 variance function structure.}
 \item{header}{an optional logical value controlling whether a header
   should be included with the rest of the output. Defaults to
   \code{TRUE}.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  The variance function structure description, the formula and the
  coefficients  associated with \code{x} are printed.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{summary.varFunc}}}

\examples{
library(nlme)
data(Orthodont)
vf1 <- varPower(0.3, form = ~age)
vf1 <- initialize(vf1, Orthodont)
print(summary(vf1))
}
\keyword{models}
