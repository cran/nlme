% $Id: AIC.lmList.Rd,v 1.1 2001/10/30 20:51:14 bates Exp $
\name{AIC}
\title{Akaike Information Criterion}
\usage{
\method{AIC}{lmList}(object, \dots, k)
}
\alias{AIC.nls}
\alias{AIC.lmList}
\arguments{
 \item{object}{a fitted model object, for which there exists a
   \code{logLik} method to extract the corresponding log-likelihood, or
   an object inheriting from class \code{logLik}.}
 \item{\dots}{optional fitted model objects.}
 \item{k}{numeric, the ``penalty'' per parameter to be used; the
   default \code{k = 2} is the classical AIC.}
}
\description{
This generic function calculates the Akaike information criterion for
one or several fitted model objects for which a log-likelihood value
can be obtained, according to the formula \eqn{-2 \mbox{log-likelihood}
  + 2 n_{par}}{-2*log-likelihood + 2*npar}, where \eqn{n_{par}}{npar} 
represents the number of parameters in the fitted model. When comparing
fitted objects, the smaller the AIC, the better the fit.
}
\value{
 if just one object is provided, returns a numeric value
 with the corresponding AIC; if more than one object are provided,
 returns a \code{data.frame} with rows corresponding to the objects and
 columns representing the number of parameters in the model
 (\code{df}) and the AIC.
}
\references{
  Sakamoto, Y., Ishiguro, M., and Kitagawa G. (1986) "Akaike
  Information Criterion Statistics", D. Reidel Publishing Company.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{logLik}}, \code{\link{BIC}}, \code{\link{AIC.logLik}}}.
\examples{
data(Orthodont)
fm1 <- lm(distance ~ age, data = Orthodont) # no random effects
AIC(fm1)
%\dontrun{
fm2 <- lme(distance ~ age, data = Orthodont) # random is ~age
AIC(fm1, fm2)
%}
}
\keyword{models}
