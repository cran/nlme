% $Id: logLik.gnlsStruct.Rd,v 1.3 2000/07/03 18:22:47 bates Exp $
\name{logLik.gnlsStruct}
\title{Log-Likelihood of a gnlsStruct Object}
\usage{
logLik(object, Pars, conLin)
}
\alias{logLik.gnlsStruct}

\arguments{
 \item{object}{an object inheriting from class \code{gnlsStruct},
   representing a list of model components, such as
   \code{corStruct} and \code{varFunc} objects, and attributes
   specifying the underlying nonlinear model and the response variable.}
 \item{Pars}{the parameter values at which the log-likelihood is to be
   evaluated.} 
 \item{conLin}{an optional condensed linear model object, consisting of
   a list with components \code{"Xy"}, corresponding to a regression
   matrix (\code{X}) combined with a response vector (\code{y}), and 
   \code{"logLik"}, corresponding to the log-likelihood of the
   underlying nonlinear model. Defaults to \code{attr(object,
     "conLin")}.} 
}
\description{
  \code{Pars} is used to update the coefficients of the model components
  of \code{object} and the individual log-likelihood
  contributions of each component are added together.
}
\value{
  the log-likelihood for the linear model described by \code{object},
  evaluated at \code{Pars}.  
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gnls}}, \code{\link{gnlsStruct}}}

\keyword{models}
