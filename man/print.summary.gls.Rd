% $Id: print.summary.gls.Rd,v 1.6 2001/06/18 21:16:50 bates Exp $
\name{print.summary.gls}
\title{Print a summary.gls Object}
\usage{
\method{print}{summary.gls}(x, verbose, digits, ...)
}
\alias{print.summary.gls}
\arguments{
 \item{x}{an object inheriting from class \code{summary.gls},
   representing a summarized \code{gls} object.}
 \item{verbose}{an optional logical value used to control the amount of
   printed output. Defaults to \code{FALSE}.}
 \item{digits}{an optional integer value specifying the number of
   significant digits to print for numeric values.  Defaults to
   \code{options("digits").  Probably not used.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  Information summarizing the fitted linear model represented by
  \code{x} is printed. This includes the AIC, BIC, and 
  log-likelihood at convergence, the coefficient estimates and
  their respective standard errors,  correlation and variance
  function parameters, if any are present, and the residual standard
  error. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{summary.gls}}, \code{\link{gls}} }

\examples{
data(Ovary)
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
print(summary(fm1))
}
\keyword{models}
