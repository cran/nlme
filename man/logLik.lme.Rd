% $Id: logLik.lme.Rd,v 1.7.2.1 2003/04/17 22:28:45 bates Exp $
\name{logLik.lme}
\title{Log-Likelihood of an lme Object}
\usage{
\method{logLik}{lme}(object, REML, \dots)
}
\alias{logLik.lme}
\alias{logLik.gls}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
 \item{REML}{an optional logical value. If \code{TRUE} the restricted
   log-likelihood is returned, else, if \code{FALSE}, the log-likelihood
   is returned. Defaults to \code{FALSE}. }
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  If \code{REML=FALSE}, returns the log-likelihood value of the linear
  mixed-effects model represented by \code{object} evaluated at the
  estimated  coefficients; else, the restricted log-likelihood evaluated
  at the estimated coefficients is returned.
}
\value{
  the (restricted) log-likelihood of the model represented by
  \code{object} evaluated at the estimated coefficients.
}
\references{
  Harville, D.A. (1974) ``Bayesian Inference for Variance Components
  Using Only Error Contrasts'', \emph{Biometrika}, \bold{61}, 383--385.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{lme}},\code{\link{gls}}}
\examples{
fm1 <- lme(distance ~ Sex * age, Orthodont, random = ~ age, method = "ML")
logLik(fm1)
logLik(fm1, REML = TRUE)
}
\keyword{models}
