% $Id: print.lme.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{print.lme}
\title{Print an lme Object}
\usage{
\method{print}{lme}(x, ...)
}
\alias{print.lme}
\arguments{
 \item{x}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  Information describing the fitted linear mixed-effects model
  represented by \code{x} is printed. This includes the fixed
  effects, the standard deviations and correlations for the
  random effects, the within-group correlation and variance
  function parameters, if any are present, and the within-group standard
  deviation. 
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{lme}}, \code{\link{print.summary.lme}} }

\examples{
data(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
print(fm1)
}
\keyword{models}
