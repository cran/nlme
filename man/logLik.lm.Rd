% $Id: logLik.lm.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{logLik.lm}
\title{Extract Log-Likelihood from an lm Object}
\usage{
logLik(object, REML)
}
\alias{logLik.lm}
\arguments{
 \item{object}{an object inheriting from class \code{lm}.}
 \item{REML}{an optional logical value. If \code{TRUE} the restricted
   log-likelihood is returned, else, if \code{FALSE}, the log-likelihood
   is returned. Defaults to \code{FALSE}.}
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

\seealso{\code{\link{lm}}}

\examples{
library(nlme)
data(Orthodont)
fm1 <- lm(distance ~ Sex * age, Orthodont)
logLik(fm1)
logLik(fm1, REML = TRUE)
}
\keyword{models}
