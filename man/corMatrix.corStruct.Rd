% $Id: corMatrix.corStruct.Rd,v 1.6 2001/06/18 21:16:50 bates Exp $
\name{corMatrix.corStruct}
\title{Matrix of a corStruct Object}
\usage{
\method{corMatrix}{corStruct}(object, covariate, corr)
}
\alias{corMatrix.corStruct}
\alias{corMatrix.compSymm}
\alias{corMatrix.corAR1}
\alias{corMatrix.corARMA}
\alias{corMatrix.corCAR1}
\alias{corMatrix.corCompSymm}
\alias{corMatrix.corIdent}
\alias{corMatrix.corNatural}
\alias{corMatrix.corSpatial}
\alias{corMatrix.corSymm}
\arguments{
 \item{object}{an object inheriting from class \code{corStruct}
   representing a correlation structure.}
 \item{covariate}{an optional covariate vector (matrix), or list of
   covariate vectors (matrices), at which values the correlation matrix,
   or list of correlation  matrices, are to be evaluated. Defaults to
   \code{getCovariate(object)}.} 
 \item{corr}{a logical value. If \code{TRUE} the function returns the
   correlation matrix, or list of correlation matrices, represented by
   \code{object}. If \code{FALSE} the function returns a transpose
   inverse square-root of the correlation matrix, or a list of transpose
   inverse square-root factors of the correlation matrices.}
}
\description{
  This method function extracts the correlation matrix (or its transpose
  inverse square-root factor), or list of correlation matrices (or their
  transpose inverse square-root factors) corresponding to
  \code{covariate} and \code{object}. Letting \eqn{\Sigma}{S} denote
  a correlation matrix, a square-root factor of \eqn{\Sigma}{S} is
  any square matrix \eqn{L} such that \eqn{\Sigma =  L'L}{S=L'L}. When
  \code{corr = FALSE}, this method extracts \eqn{L^{-t}}{L^(-t)}. 
}
\value{
  If \code{covariate} is a vector (matrix), the returned value will be
  an array with the corresponding correlation matrix (or its transpose
  inverse square-root factor). If the \code{covariate} is a list of
  vectors (matrices), the returned value will be a list with the
  correlation matrices (or their transpose inverse square-root factors)
  corresponding to each component of \code{covariate}.  
}
\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{corFactor.corStruct}}, \code{\link{initialize.corStruct}}}

\examples{
cs1 <- corAR1(0.3)
corMatrix(cs1, covariate = 1:4)
corMatrix(cs1, covariate = 1:4, corr = F)
}
\keyword{models}
