% $Id: BIC.logLik.Rd,v 1.2 2000/04/17 19:03:38 bates Exp $
\name{BIC.logLik}
\title{BIC of a logLik Object}
\usage{
BIC(object)
}
\alias{BIC.logLik}
\arguments{
 \item{object}{an object inheriting from class \code{logLik}, usually
   resulting from applying a \code{logLik} method to a fitted model
   object.}
}
\description{
  This function calculates the Bayesian information criterion,
  also known as Schwarz's Bayesian criterion (SBC) for an object
  inheriting from class \code{logLik}, according to the formula 
  \eqn{-2 \mbox{log-likelihood} + n_{par} \log(n_{obs})}{log-likelihood +
  npar*log(nobs)}, where \eqn{n_{par}}{npar} represents the number of
  parameters  and \eqn{n_{obs}}{nobs} the number of  observations in the
  fitted model. When comparing fitted objects, the smaller the BIC, the
  better the fit. 
}
\value{
  a numeric value with the corresponding BIC.
}
\references{
  Schwarz, G. (1978) "Estimating the Dimension of a Model", Annals of
  Statistics, 6, 461-464.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{BIC}}, \code{\link{logLik}}, \code{\link{AIC}}}.
\examples{
data(Orthodont)
fm1 <- lm(distance ~ age, data = Orthodont) 
BIC(logLik(fm1))
}
\keyword{models}
