% File nlme/man/summary.lme.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{summary.lme}
\title{Summarize an lme Object}
\alias{summary.lme}
%\alias{coef.summary.lme}
\alias{print.summary.lme}
\usage{
\method{summary}{lme}(object, adjustSigma, verbose, \dots)
\method{print}{summary.lme}(x, verbose = FALSE, \dots)
}
\arguments{
  \item{object}{an object inheriting from class \code{"\link{lme}"}, representing
    a fitted linear mixed-effects model.}
  \item{adjustSigma}{an optional logical value.  If \code{TRUE} and the
    estimation method used to obtain \code{object} was maximum
    likelihood, the residual standard error is multiplied by
    \eqn{\sqrt{n_{obs}/(n_{obs} - n_{par})}}{sqrt(nobs/(nobs - npar))},
    converting it to a REML-like estimate.  This argument is only used
    when a single fitted object is passed to the function.  Default is
    \code{TRUE}.}
  \item{verbose}{an optional logical value used to control the amount of
    output in the \code{print.summary.lme} method.  Defaults to
    \code{FALSE}.}
  \item{\dots}{additional optional arguments passed to methods, mainly
    for the \code{\link{print}} method.}
  \item{x}{a \code{"summary.lme"} object.}
!}
\description{
  Additional information about the linear mixed-effects fit represented
  by \code{object} is extracted and included as components of
  \code{object}.  The returned object has a \code{\link{print}} and a
  \code{\link{coef}} method, the latter returning the coefficient's
  \code{tTtable}.
}
\value{
  an object inheriting from class \code{summary.lme} with all components
  included in \code{object} (see \code{\link{lmeObject}} for a full
  description of the components) plus the following components:
  \item{corFixed}{approximate correlation matrix for the fixed effects
    estimates.}
  \item{tTable}{a matrix with columns named \code{Value},
    \code{Std. Error}, \code{DF}, \code{t-value}, and \code{p-value}
    representing respectively the fixed effects estimates, their
    approximate standard errors, the denominator degrees of freedom, the
    ratios between the estimates and their standard errors, and the
    associated p-value from a t distribution.  Rows correspond to the
    different fixed effects.}
  \item{residuals}{if more than five observations are used in the
    \code{lme} fit, a vector with the minimum, first quartile, median, third
    quartile, and maximum of the innermost grouping level residuals
    distribution; else the innermost grouping level residuals.}
  \item{AIC}{the Akaike Information Criterion corresponding to
    \code{object}.}
  \item{BIC}{the Bayesian Information Criterion corresponding to
    \code{object}.}
}

\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{
  \code{\link{AIC}}, \code{\link{BIC}},
  \code{\link{lme}}.
}
\examples{
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
(s1 <- summary(fm1))
\dontrun{coef(s1) # the (coef | Std.E | t | P-v ) matrix}
\dontshow{stopifnot(is.matrix(coef(s1)))}
}
\keyword{models}
