% $Id: logLik.glsStruct.Rd,v 1.4 2001/06/18 21:16:50 bates Exp $
\name{logLik.glsStruct}
\title{Log-Likelihood of a glsStruct Object}
\usage{
\method{logLik}{glsStruct}(object, Pars, conLin)
}
\alias{logLik.glsStruct}

\arguments{
 \item{object}{an object inheriting from class \code{glsStruct},
   representing a list of linear model components, such as
   \code{corStruct} and \code{varFunc} objects.}
 \item{Pars}{the parameter values at which the (restricted)
   log-likelihood is to be evaluated.}
 \item{conLin}{an optional condensed linear model object, consisting of
   a list with components \code{"Xy"}, corresponding to a regression
   matrix (\code{X}) combined with a response vector (\code{y}), and 
   \code{"logLik"}, corresponding to the log-likelihood of the
   underlying linear model. Defaults to \code{attr(object, "conLin")}.}
}
\description{
  \code{Pars} is used to update the coefficients of the model components
  of \code{object} and the individual (restricted) log-likelihood
  contributions of each component are added together. The type of
  log-likelihood (restricted or not) is determined by the
  \code{settings} attribute of \code{object}.
}
\value{
  the (restricted) log-likelihood for the linear model described by
  \code{object}, evaluated at \code{Pars}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gls}}, \code{\link{glsStruct}}}

\keyword{models}
