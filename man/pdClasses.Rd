% File nlme/man/pdClasses.Rd
% Part of the nlme package for R
% Distributed under GPL 2 or later: see nlme/LICENCE.note

\name{pdClasses}
\title{Positive-Definite Matrix Classes}
\alias{pdClasses}
\description{
  Standard classes of positive-definite matrices (\code{pdMat})
  available in the \pkg{nlme} package.
}
\value{
  Available standard classes:
  \item{pdSymm}{general positive-definite matrix, with no additional
    structure}
  \item{pdLogChol}{general positive-definite matrix, with no additional
    structure, using a log-Cholesky parameterization}
  \item{pdDiag}{diagonal}
  \item{pdIdent}{multiple of an identity}
  \item{pdCompSymm}{compound symmetry structure (constant diagonal and
    constant off-diagonal elements)}
  \item{pdBlocked}{block-diagonal matrix, with diagonal blocks of any
    "atomic" \code{pdMat} class}
  \item{pdNatural}{general positive-definite matrix in natural
    parametrization (i.e. parametrized in terms of standard deviations
    and correlations). The underlying coefficients are not unrestricted,
    so this class should NOT be used for optimization.}
}
\references{
 Pinheiro, J.C., and Bates, D.M. (2000) "Mixed-Effects Models
 in S and S-PLUS", Springer.  
}
\author{José Pinheiro and Douglas Bates \email{bates@stat.wisc.edu}}
\note{
  Users may define their own \code{pdMat} classes by specifying a
  constructor function and, at a minimum, methods for the
  functions \code{\link{pdConstruct}}, \code{\link{pdMatrix}}, and \code{\link{coef}}.
  For examples of these functions, see the methods for classes \code{pdSymm}
  and \code{pdDiag}. 
}
\seealso{
  \code{\link{pdBlocked}},
  \code{\link{pdCompSymm}},
  \code{\link{pdDiag}},
  \code{\link{pdFactor}}, 
  \code{\link{pdIdent}},
  \code{\link{pdMat}},
  \code{\link{pdMatrix}},
  \code{\link{pdNatural}},
  \code{\link{pdSymm}},
  \code{\link{pdLogChol}}
}
\keyword{models}
