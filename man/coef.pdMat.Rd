% $Id: coef.pdMat.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{coef.pdMat}
\title{pdMat Object Coefficients}
\usage{
\method{coef}{pdMat}(object, unconstrained)
coef(object) <- value
}
\alias{coef.pdMat}
\alias{coef.pdBlocked}
\alias{coef.pdNatural}
\alias{coef.pdSymm}
\alias{coef<-.pdMat}
\alias{coef<-.pdBlocked}
\arguments{
 \item{object}{an object inheriting from class \code{pdMat},
   representing a positive-definite matrix.}
 \item{unconstrained}{a logical value. If \code{TRUE} the coefficients
   are returned in unconstrained form (the same used in the optimization
   algorithm). If \code{FALSE} the upper triangular elements of the
   positive-definite matrix represented by \code{object} are
   returned. Defaults to \code{TRUE}.}
 \item{value}{a vector with the replacement values for the coefficients
   associated with \code{object}. It must be a vector with the same length
   of \code{coef{object}} and must be given in unconstrained form.}
}
\description{
  This method function extracts the coefficients associated with the
  positive-definite matrix represented by \code{object}. 
}
\value{
  a vector with the coefficients corresponding to \code{object}.
}
\references{
 Pinheiro, J.C. and Bates., D.M.  (1996) "Unconstrained Parametrizations
 for Variance-Covariance Matrices", Statistics and Computing, 6, 289-296.
}
\author{Jose Pinheiro and Douglas Bates }
\section{SIDE EFFECTS}{
  On the left side of an assignment, sets the values of the coefficients
  of \code{object} to \code{value}. 
}

\seealso{\code{\link{pdMat}}}

\examples{
coef(pdSymm(diag(3)))
}
\keyword{models}
