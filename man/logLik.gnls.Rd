% $Id: logLik.gnls.Rd,v 1.6 2001/06/18 21:16:50 bates Exp $
\name{logLik.gnls}
\title{Log-Likelihood of a gnls Object}
\usage{
\method{logLik}{gnls}(object, REML)
}
\alias{logLik.gnls}
\arguments{
 \item{object}{an object inheriting from class \code{gnls}, representing
   a generalized nonlinear least squares fitted model.}
 \item{REML}{an optional logical value. If \code{TRUE} the restricted
   log-likelihood is returned, else, if \code{FALSE}, the log-likelihood
   is returned. Defaults to \code{FALSE}. }
}
\description{
  Returns the log-likelihood value of the nonlinear model represented by
  \code{object} evaluated at the estimated coefficients.  
}
\value{
  the log-likelihood of the linear model represented by
  \code{object} evaluated at the estimated coefficients.  
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{gnls}}}
\examples{
data(Soybean)
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
logLik(fm1)
}
\keyword{models}
