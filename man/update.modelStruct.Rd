% $Id: update.modelStruct.Rd,v 1.1 2000/03/17 22:21:24 saikat Exp $
\name{update.modelStruct}
\title{Update a modelStruct Object}
\usage{
update(object, data)
}
\alias{update.modelStruct}
\arguments{
    \item{object}{an object inheriting from class \code{modelStruct},
	representing a list of model components, such as \code{corStruct} and
	\code{varFunc} objects.}
    \item{data}{a data frame in which to evaluate the variables needed for
	updating the elements of \code{object}.}
}
\description{
    This method function updates each element of \code{object}, allowing
    the access to \code{data}.
}
\value{
    an object similar to \code{object} (same class, length, and names),
    but with updated elements.
}
\author{Jose Pinheiro and Douglas Bates}
\note{This method function is primarily used within model fitting
  functions, such as \code{lme} and \code{gls}, that allow model
  components such as variance functions.}
\keyword{models}
