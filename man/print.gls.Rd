% $Id: print.gls.Rd,v 1.6 2001/06/18 21:16:50 bates Exp $
\name{print.gls}
\title{Print a gls Object}
\usage{
\method{print}{gls}(x, ...)
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

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{gls}}, \code{\link{print.summary.gls}} }

\examples{
data(Ovary)
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
print(fm1)
}
\keyword{models}
\keyword{print}
