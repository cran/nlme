% $Id: initialize.reStruct.Rd,v 1.3 2000/07/03 18:22:47 bates Exp $
\name{initialize.reStruct}
\title{Initialize reStruct Object}
\usage{
initialize(object, data, conLin, control)
}
\alias{initialize.reStruct}
\arguments{
 \item{object}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
 \item{data}{a data frame in which to evaluate the variables defined in
   \code{formula(object)}.}
 \item{conLin}{a condensed linear model object, consisting of a list
   with components \code{"Xy"}, corresponding to a regression matrix
   (\code{X}) combined with a response vector (\code{y}), and
   \code{"logLik"}, corresponding to the log-likelihood of the
   underlying model.} 
 \item{control}{an optional list with a single component \code{niterEM}
   controlling the number of iterations for the EM algorithm used to
   refine initial parameter estimates. It is given as a list for
   compatibility with other \code{initialize} methods. Defaults to
   \code{list(niterEM = 20)}.}
}
\description{
  Initial estimates for the parameters in the \code{pdMat} objects
  forming \code{object}, which have not yet been initialized, are
  obtained using the methodology described in Bates and Pinheiro
  (1998). These estimates may be refined using a series of EM
  iterations, as described in Bates and Pinheiro (1998). The number of
  EM iterations to be used is defined in \code{control}. 
}
\value{
  an \code{reStruct} object similar to \code{object}, but with all
  \code{pdMat} components initialized.
}
\references{
  Bates, D.M. and Pinheiro, J.C. (1998) "Computational methods for
  multilevel models" available in PostScript or PDF formats at
  http://franz.stat.wisc.edu/pub/NLME/
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{reStruct}}, \code{\link{pdMat}}}

\keyword{models}
