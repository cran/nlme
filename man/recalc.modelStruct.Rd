% $Id: recalc.modelStruct.Rd,v 1.3 2000/07/03 18:22:48 bates Exp $
\name{recalc.modelStruct}
%\title{Recalculate Condensed Linear Model for a \code{modelStruct} Object}
\title{Recalculate for a modelStruct Object}
\usage{
recalc(object, conLin)
}
\alias{recalc.modelStruct}
\arguments{
 \item{object}{an object inheriting from class \code{modelStruct},
   representing a list of model components, such as \code{corStruct} and
   \code{varFunc} objects.}
 \item{conLin}{an optional  condensed linear model object, consisting of
   a list with components \code{"Xy"}, corresponding to a regression
   matrix (\code{X}) combined with a response vector (\code{y}), and 
   \code{"logLik"}, corresponding to the log-likelihood of the
   underlying model. Defaults to \code{attr(object, "conLin")}.}
}
\description{
  This method function recalculates the condensed linear model object
  using each element of \code{object} sequentially from last to first.
}
\value{
  the recalculated condensed linear model object.
}
\author{Jose Pinheiro and Douglas Bates}
\note{This method function is generally only used inside model fitting
  functions, such as \code{lme} and \code{gls}, that allow model
  components, such as correlated error terms and variance functions.
}
\seealso{\code{\link{recalc.corStruct}}, \code{\link{recalc.reStruct}},
    \code{\link{recalc.varFunc}}}
\keyword{models}
