% $Id: logLik.gnls.Rd,v 1.8 2002/03/05 14:59:39 bates Exp $
\name{logLik.gnls}
\title{Log-Likelihood of a gnls Object}
\usage{
\method{logLik}{gnls}(object, REML, \dots)
}
\alias{logLik.gnls}
\arguments{
 \item{object}{an object inheriting from class \code{gnls}, representing
   a generalized nonlinear least squares fitted model.}
 \item{REML}{an optional logical value. If \code{TRUE} the restricted
   log-likelihood is returned, else, if \code{FALSE}, the log-likelihood
   is returned. Defaults to \code{FALSE}. }
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  Returns the log-likelihood value of the nonlinear model represented by
  \code{object} evaluated at the estimated coefficients.  
}
\value{
  the log-likelihood of the linear model represented by
  \code{object} evaluated at the estimated coefficients.  
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{gnls}},
  \code{\link{logLik.lme}}
}
\examples{
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
logLik(fm1)
}
\keyword{models}
