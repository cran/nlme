% $Id: recalc.varFunc.Rd,v 1.1 2000/03/17 22:21:24 saikat Exp $
\name{recalc.varFunc}
\title{Recalculate for varFunc Object}
\usage{
recalc(object, conLin)
}
\alias{recalc.varFunc}
\arguments{
 \item{object}{an object inheriting from class \code{varFunc},
   representing a variance function structure.} 
 \item{conLin}{a condensed linear model object, consisting of a list
   with components \code{"Xy"}, corresponding to a regression matrix
   (\code{X}) combined with a response vector (\code{y}), and
   \code{"logLik"}, corresponding to the log-likelihood of the
   underlying model.} 
}
\description{
  This method function pre-multiples the \code{"Xy"} component of
  \code{conLin} by a diagonal matrix with diagonal elements given by the
  weights corresponding to the variance structure represented by
  \code{object}e and adds the log-likelihood contribution of
  \code{object}, given by \code{logLik(object)}, to the \code{"logLik"}
  component of \code{conLin}.  
}
\value{
  the recalculated condensed linear model object.
}
\author{Jose Pinheiro and Douglas Bates}
\note{This method function is only used inside model fitting functions,
    such as \code{lme} and \code{gls}, that allow heteroscedastic error
    terms.}
\seealso{\code{\link{varWeights}}, \code{\link{logLik.varFunc}}}
\keyword{models}
