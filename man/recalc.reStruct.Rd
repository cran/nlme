% $Id: recalc.reStruct.Rd,v 1.3 2000/07/03 18:22:48 bates Exp $
\name{recalc.reStruct}
\title{Recalculate for an reStruct Object}
\usage{
recalc(object, conLin)
}
\alias{recalc.reStruct}
\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
 \item{conLin}{a condensed linear model object, consisting of a list
   with components \code{"Xy"}, corresponding to a regression matrix
   (\code{X}) combined with a response vector (\code{y}), and
   \code{"logLik"}, corresponding to the log-likelihood of the
   underlying model.} 
}
\description{
  The log-likelihood, or restricted log-likelihood, of the
  Gaussian linear mixed-effects model represented by \code{object} and
  \code{conLin} (assuming spherical within-group covariance structure),
  evaluated at \code{coef(object)} is calculated and added to the
  \code{logLik} component of \code{conLin}. The \code{settings}
  attribute of \code{object} determines whether the log-likelihood, or
  the restricted log-likelihood, is to be calculated. The computational
  methods for the (restricted) log-likelihood calculations are described
  in Bates and Pinheiro (1998). 
}
\value{
  the condensed linear model with its \code{logLik} component updated.
}

\references{
  Bates, D.M. and Pinheiro, J.C. (1998) "Computational methods for
  multilevel models" available in PostScript or PDF formats at
  http://franz.stat.wisc.edu/pub/NLME/
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{reStruct}}, \code{\link{logLik}},
  \code{\link{lme}}} 

\keyword{models}
