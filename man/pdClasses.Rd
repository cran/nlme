% $Id: pdClasses.Rd,v 1.5 2001/06/20 16:55:20 bates Exp $
\name{pdClasses}
\title{Positive-Definite Matrix Classes}
\alias{pdClasses}
\description{
  Standard classes of positive-definite matrices (\code{pdMat})
  structures  available in the \code{nlme} library. 
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
\author{Jose Pinheiro and Douglas Bates}
\note{
  Users may define their own \code{pdMat} classes by specifying a
  \code{constructor} function and, at a minimum, methods for the
  functions \code{pdConstruct}, \code{pdMatrix} and \code{coef}. For
  examples of these functions, see the methods for classes \code{pdSymm}
  and \code{pdDiag}. 
}
\seealso{\code{\link{pdCompSymm}}, \code{\link{pdDiag}},
  \code{\link{pdIdent}}, \code{\link{pdNatural}}, \code{\link{pdSymm}},
  \code{\link{pdLogChol}}
}
\keyword{models}
