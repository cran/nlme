% $Id: AIC.logLik.Rd,v 1.2 1999/11/04 16:38:58 saikat Exp $
\name{AIC.logLik}
\title{AIC of a logLik Object}
\usage{
AIC(object)
}
\alias{AIC.logLik}
\arguments{
 \item{object}{an object inheriting from class \code{logLik}, usually
   resulting from applying a \code{logLik} method to a fitted model
   object.}
}
\description{
This function calculates the Akaike information criterion for
an object inheriting from class \code{logLik}, according to the formula
\eqn{-2 \mbox{log-likelihood} + 2 n_{par}}{2*log-likelihood + 2*npar}, where
\eqn{n_{par}}{npar} represents the number of parameters in the fitted 
model. When comparing fitted objects, the smaller the AIC, the better
the fit. 
}
\value{
  a numeric value with the corresponding AIC.
}
\references{
  Sakamoto, Y., Ishiguro, M., and Kitagawa G. (1986) "Akaike
  Information Criterion Statistics", D. Reidel Publishing Company.
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{AIC}}, \code{\link{logLik}}, \code{\link{BIC}}}.
\examples{
library(nlme)
data(Orthodont)
fm1 <- lm(distance ~ age, data = Orthodont) 
AIC(logLik(fm1))
}
\keyword{models}
