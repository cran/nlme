% $Id: print.summary.lme.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{print.summary.lme}
\title{Print a summary.lme Object}
\usage{
print(x, verbose, ...)
}
\alias{print.summary.lme}
\arguments{
 \item{x}{an object inheriting from class \code{summary.lme},
   representing a summarized \code{lme} object.}
 \item{verbose}{an optional logical value used to control the amount of
   printed output. Defaults to \code{FALSE}.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  Information summarizing the fitted linear mixed-effects model
  represented by \code{x} is printed. This includes the AIC, BIC, and
  log-likelihood at convergence, the fixed effects estimates and
  respective standard errors, the standard deviations and correlations
  for the random effects, the within-group correlation and variance
  function parameters, if any are present, and the within-group standard
  deviation. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{summary.lme}}, \code{\link{lme}} }

\examples{
library(nlme)
data(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
print(summary((fm1)))
}
\keyword{models}
