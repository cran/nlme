% $Id: allCoef.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{allCoef}
\title{Extract Coefficients from a Set of Objects}
\usage{
allCoef(..., extract)
}
\alias{allCoef}
\arguments{
 \item{\dots}{objects to which \code{extract} will be applied. Generally
 these will be model components, such as \code{corStruct} and
 \code{varFunc} objects.}
 \item{extract}{an optional extractor function. Defaults to \code{coef}.}
}
\description{
  The extractor function is applied to each object in \code{...}, with
  the result being converted to a vector. A \code{map} attribute is
  included to indicate which pieces of the returned vector correspond to
  the original objects in \code{...}.
}
\value{
  a vector with all elements, generally coefficients, obtained by
  applying \code{extract} to the objects in \code{...}. 
}
\author{Jose' Pinheiro and Douglas Bates}

\seealso{\code{\link{lmeStruct}},\code{\link{nlmeStruct}}}

\examples{
library(nlme)
cs1 <- corAR1(0.1)
vf1 <- varPower(0.5)
allCoef(cs1, vf1)
}
\keyword{models}

