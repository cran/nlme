% $Id: formula.modelStruct.Rd,v 1.7 2002/03/05 14:59:39 bates Exp $
\name{formula.modelStruct}
\title{Extract modelStruct Object Formula}
\usage{
\method{formula}{modelStruct}(x, \dots)
}
\alias{formula.modelStruct}
\arguments{
 \item{x}{an object inheriting from class \code{modelStruct},
   representing a list of model components, such as \code{corStruct} and
   \code{varFunc} objects.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  This method function extracts a formula from each of the 
  components of \code{x}, returning a list of formulas.
}
\value{
  a list with the formulas of each component of \code{x}.
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{formula}}}

\examples{
lms1 <- lmeStruct(reStruct = reStruct(pdDiag(diag(2), ~age)),
   corStruct = corAR1(0.3))
formula(lms1)
}
\keyword{models}
