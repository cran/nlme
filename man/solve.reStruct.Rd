% $Id: solve.reStruct.Rd,v 1.2 1999/11/04 16:38:59 saikat Exp $
\name{solve.reStruct}
\title{Apply Solve to an reStruct Object}
\usage{
solve(a, b, tol)
}
\alias{solve.reStruct}
\arguments{
 \item{a}{an object inheriting from class \code{reStruct},
   representing a random effects structure and consisting of a list of
   \code{pdMat} objects.} 
 \item{b}{this argument is only included for consistency with the
   generic function and is not used in this method function.}
 \item{tol}{an optional numeric value for the tolerance used in the
   numerical algorithm. Defaults to \code{1e-7}.}
}
\description{
  \code{Solve} is applied to each \code{pdMat} component of \code{a},
  which results in inverting the positive-definite matrices they
  represent. 
}
\value{
  an \code{reStruct} object similar to \code{a}, but with the
  \code{pdMat} components representing the inverses of the
  matrices represented by the components of \code{a}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{solve.pdMat}}, \code{\link{reStruct}}}

\examples{
library(nlme)
rs1 <- reStruct(list(A = pdSymm(diag(1:3), form = ~Score),
  B = pdDiag(2 * diag(4), form = ~Educ)))
solve(rs1)
}
\keyword{models}
