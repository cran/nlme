% $Id: print.reStruct.Rd,v 1.2 2000/04/17 19:03:40 bates Exp $
\name{print.reStruct}
\title{Print an reStruct Object}
\usage{
print(x, sigma, reEstimates, verbose=F, ...)
}
\alias{print.reStruct}
\arguments{
 \item{x}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
 \item{sigma}{an optional numeric value used as a multiplier for
   the square-root factors of the \code{pdMat} components (usually the
   estimated within-group standard deviation from a mixed-effects
   model). Defaults to 1.} 
 \item{reEstimates}{an optional list with the random effects estimates
   for each level of grouping. Only used when \code{verbose = TRUE}.}
 \item{verbose}{an optional logical value determining if the random
   effects estimates should be printed. Defaults to \code{FALSE}.}
 \item{\dots}{optional arguments passed to \code{print.default}; see
   the documentation on that method function.}
}
\description{
  Each \code{pdMat} component of \code{object} is printed, together with
  its formula and the associated grouping level. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{reStruct}}}

\examples{
data(Pixel)
rs1 <- reStruct(list(Dog = ~day, Side = ~1), data = Pixel)
matrix(rs1) <- list(diag(2), 3)
print(rs1)
}
\keyword{models}
