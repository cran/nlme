% $Id: varWeights.lmeStruct.Rd,v 1.1 1999/10/13 00:50:13 saikat Exp $
\name{varWeights.lmeStruct}
\title{Variance Weights for lmeStruct Object}
\usage{
varWeights(object)
}
\alias{varWeights.lmeStruct}
\arguments{
  \item{object}{an object inheriting from class \code{lmeStruct},
    representing a list of linear mixed-effects model components, such as
    \code{reStruct}, \code{corStruct}, and \code{varFunc} objects.}
}
\description{
  If \code{object} includes a \code{varStruct} component, the inverse of
  the standard deviations of the variance function structure represented
  by the corresponding \code{varFunc} object are returned; else, a
  vector of ones of length equal to the number of observations in the
  data frame used to fit the associated linear mixed-effects model is
  returned. 
}
\value{
  if \code{object} includes a \code{varStruct} component, a vector with
  the corresponding variance weights; else, or a vector of ones.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{varWeights}}}

\keyword{models}
