% $Id: logLik.gls.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{logLik.gls}
\title{Log-Likelihood of a gls Object}
\usage{
logLik(object, REML)
}
\alias{logLik.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
 \item{REML}{an optional logical value. If \code{TRUE} the restricted
   log-likelihood is returned, else, if \code{FALSE}, the log-likelihood
   is returned. Defaults to \code{FALSE}. }
}
\description{
  If \code{REML=FALSE}, returns the log-likelihood value of the linear
  model represented by \code{object} evaluated at the estimated
  coefficients; else, the restricted log-likelihood evaluated at the
  estimated coefficients is returned. 
}
\value{
  the (restricted) log-likelihood of the linear model represented by
  \code{object} evaluated at the estimated coefficients.  
}
\references{
  Harville, D.A. (1974) "Bayesian Inference for Variance Components
  Using Only Error Contrasts", Biometrika, 61, 383-385.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gls}}}

\examples{
data(Ovary)
%\dontrun{
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
logLik(fm1)
logLik(fm1, REML = FALSE)
%}
}
\keyword{models}
