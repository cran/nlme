% $Id: formula.modelStruct.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{formula.modelStruct}
%\title{Extract Formula from \code{modelStruct} Object}
\title{Extract modelStruct Object Formula}
\usage{
\method{formula}{modelStruct}(object)
}
\alias{formula.modelStruct}
\arguments{
 \item{object}{an object inheriting from class \code{modelStruct},
   representing a list of model components, such as \code{corStruct} and
   \code{varFunc} objects.}
}
\description{
  This method function extracts a formula from each of the 
  components of \code{object}, returning a list of formulas.
}
\value{
  a list with the formulas of each component of \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{formula}}}

\examples{
lms1 <- lmeStruct(reStruct = reStruct(pdDiag(diag(2), ~age)),
   corStruct = corAR1(0.3))
formula(lms1)
}
\keyword{models}
