% $Id: print.gls.Rd,v 1.5 2001/01/10 19:21:05 bates Exp $
\name{print.gls}
\title{Print a gls Object}
\usage{
print(x, ...)
}
\alias{print.gls}
\arguments{
 \item{x}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  Information describing the fitted linear model represented by \code{x}
  is printed. This includes the coefficients, correlation and variance
  function parameters, if any are present, and the residual standard
  error. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gls}}, \code{\link{print.summary.gls}} }

\examples{
data(Ovary)
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
print(fm1)
}
\keyword{models}
\keyword{print}
