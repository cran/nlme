\name{MEestimate}
\alias{MEestimate}
\title{Estimate mixed-effects parameters}
\description{
  Estimate the parameters in a linear mixed-effects model.
}
\usage{
MEestimate(object, groups, conLin=attr(object, "conLin"))
}
\arguments{
  \item{object}{an object that inherits from class \code{lmeStruct},
    describing a linear mixed-effects model or a nonlinear mixed-effects model.}
  \item{groups}{a factor or list of factors giving the grouping of the observations}
  \item{conLin}{a condensed linear model object.  The default is
    \code{attr(object, "conLin")}.}
}
\details{
  This is a utility function used by \code{lme} and \code{nlme}.  It is
  not intended for direct use.}
\value{
  A list with components
  \item{logLik}{the log-likelihood at the estimate.}
  \item{b}{the estimated BLUP's (Best Linear Unbiased Predictors) of the
  random effects}
  \item{beta}{the estimates of the fixed effects}
  \item{sigma}{the estimated standard deviation of the within-subject
    noise term.}
  \item{varFix}{the estimated variance-covariance matrix of the fixed effects}
}
\author{Jose Pinheiro \email{jose.pinheiro@pharma.novartis.com}, Douglas Bates
  \email{bates@stat.wisc.edu}} 
\seealso{\code{\link{MEEM}}, \code{\link{MEdecomp}}, \code{\link{lme}}, \code{\link{nlme}}}
\keyword{models}
