% $Id: BIC.Rd,v 1.4 2000/07/03 18:22:45 bates Exp $
\name{BIC}
\title{Bayesian Information Criterion}
\usage{
BIC(object, ...)
}
\alias{BIC}
\alias{BIC.lm}
\alias{BIC.lmList}
\alias{BIC.gls}
\alias{BIC.lme}
\alias{BIC.nls}
\alias{BIC.nlsList}
\alias{BIC.nlme}
\arguments{
 \item{object}{a fitted model object, for which there exists a
   \code{logLik} method to extract the corresponding log-likelihood, or
   an object inheriting from class \code{logLik}.}
 \item{\dots}{optional fitted model objects.}
}
\description{
  This generic function calculates the Bayesian information criterion,
  also known as Schwarz's Bayesian criterion (SBC), for one or several
  fitted model objects for which a log-likelihood value can be obtained,
  according to the formula \eqn{-2 \mbox{log-likelihood} + n_{par}
    \log(n_{obs})}{-2*log-likelihood + npar*log(nobs)}, where
  \eqn{n_{par}}{npar}  represents the
  number of parameters and \eqn{n_{obs}}{nobs} the number of
  observations in the fitted model. 
}
\value{
 if just one object is provided, returns a numeric value with the
 corresponding BIC; if more than one object are provided, returns a
 \code{data.frame} with rows corresponding to the objects and columns
 representing the number of parameters in the model (\code{df}) and the
 BIC.
}
\references{
  Schwarz, G. (1978) "Estimating the Dimension of a Model", Annals of
  Statistics, 6, 461-464.
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{logLik}}, \code{\link{AIC}}, \code{\link{BIC.logLik}}}
\examples{
data(Orthodont)
fm1 <- lm(distance ~ age, data = Orthodont) # no random effects
BIC(fm1)
%\dontrun{
fm2 <- lme(distance ~ age, data = Orthodont) # random is ~age
BIC(fm1, fm2)
%}
}
\keyword{models}
