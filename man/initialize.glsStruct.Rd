% File nlme/man/Initialize.glsStruct.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{Initialize.glsStruct}
\title{Initialize a glsStruct Object}
\usage{
\method{Initialize}{glsStruct}(object, data, control, \dots)
}
\alias{Initialize.glsStruct}
\arguments{
 \item{object}{an object inheriting from class \code{"\link{glsStruct}"},
   representing a list of linear model components, such as
   \code{corStruct} and \code{varFunc} objects.}
 \item{data}{a data frame in which to evaluate the variables defined in
   \code{formula(object)}.}
 \item{control}{an optional list with control parameters for the
   initialization and optimization algorithms used in
   \code{gls}. Defaults to \code{list(singular.ok = FALSE)}, implying that linear dependencies are not allowed in the model.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  The individual linear model components of the \code{glsStruct} list
  are initialized.  
}
\value{
  a \code{glsStruct} object similar to \code{object}, but with
  initialized model components. 
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}
\seealso{\code{\link{gls}},
  \code{\link{Initialize.corStruct}},
  \code{\link{Initialize.varFunc}},
  \code{\link{Initialize}}}
\keyword{models}
