% $Id: initialize.lmeStruct.Rd,v 1.1 1999/10/13 00:50:11 saikat Exp $
\name{initialize.lmeStruct}
\title{Initialize an lmeStruct Object}
\usage{
initialize(object, data, groups, conLin, control)
}
\alias{initialize.lmeStruct}
\arguments{
 \item{object}{an object inheriting from class \code{lmeStruct},
   representing a list of linear mixed-effects model components, such as
   \code{reStruct}, \code{corStruct}, and \code{varFunc} objects.}
 \item{data}{a data frame in which to evaluate the variables defined in
   \code{formula(object)}.}
 \item{groups}{a data frame with the grouping factors corresponding to
   the lme model associated with \code{object} as columns, sorted from
   innermost to outermost grouping level.}
 \item{conLin}{an optional condensed linear model object, consisting of
   a list with components \code{"Xy"}, corresponding to a regression
   matrix (\code{X}) combined with a response vector (\code{y}), and 
   \code{"logLik"}, corresponding to the log-likelihood of the
   underlying lme model. Defaults to \code{attr(object, "conLin")}.}
 \item{control}{an optional list with control parameters for the
   initialization and optimization algorithms used in
   \code{lme}. Defaults to \code{list(niterEM=20, gradHess=TRUE)},
   implying that 20 EM iterations are to be used in the derivation of
   initial estimates for the coefficients of the \code{reStruct}
   component of \code{object} and, if possible, numerical gradient
   vectors and Hessian matrices for the log-likelihood function are to
   be used in the optimization algorithm.}
}
\description{
  The individual linear mixed-effects model components of the
  \code{lmeStruct} list are initialized. 
}
\value{
  an \code{lmeStruct} object similar to \code{object}, but with
  initialized model components. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lme}}, \code{\link{initialize.reStruct}},
  \code{\link{initialize.corStruct}} , \code{\link{initialize.varFunc}}}

\keyword{models}
