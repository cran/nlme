% $Id: varWeights.glsStruct.Rd,v 1.1 1999/10/13 00:50:13 saikat Exp $
\name{varWeights.glsStruct}
\title{Variance Weights for glsStruct Object}
\usage{
varWeights(object)
}
\alias{varWeights.glsStruct}
\arguments{
 \item{object}{an object inheriting from class \code{glsStruct},
   representing a list of linear model components, such as
   \code{corStruct} and \code{varFunc} objects.}
}
\description{
  If \code{object} includes a \code{varStruct} component, the inverse of
  the standard deviations of the variance function structure represented
  by the corresponding \code{varFunc} object are returned; else, a
  vector of ones of length equal to the number of observations in the
  data frame used to fit the associated linear model is returned. 
}
\value{
  if \code{object} includes a \code{varStruct} component, a vector with
  the corresponding variance weights; else, or a vector of ones.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{varWeights}}}

\keyword{models}
