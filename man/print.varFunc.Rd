% File nlme/man/print.varFunc.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{print.varFunc}
\title{Print a varFunc Object}
\usage{
\method{print}{varFunc}(x, \dots)
}
\alias{print.varFunc}
\alias{print.varComb}
\arguments{
  \item{x}{an object inheriting from class \code{"\link{varFunc}"},
    representing a variance function structure.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  The class and the coefficients associated with \code{x} are printed.
}

\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{summary.varFunc}}}

\examples{
vf1 <- varPower(0.3, form = ~age)
vf1 <- Initialize(vf1, Orthodont)
print(vf1)
}
\keyword{models}
