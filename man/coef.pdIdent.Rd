% $Id: coef.pdIdent.Rd,v 1.4 2000/07/03 18:22:46 bates Exp $
\name{coef.pdIdent}
\title{pdIdent Object Coefficients}
\usage{
coef(object, unconstrained)
coef(object) <- value
}
\alias{coef.pdIdent}
\arguments{
 \item{object}{an object inheriting from class \code{pdIdent},
   representing a multiple of the identity positive-definite matrix.}
 \item{unconstrained}{a logical value. If \code{TRUE} the logarithm of
   the standard deviation corresponding to the variance-covariance
   matrix represented by \code{object} is returned. If \code{FALSE} the
   standard deviation is returned. Defaults to \code{TRUE}.} 
 \item{value}{a vector with the replacement values for the coefficients
   associated with \code{object}. It must be a numeric value given in
   unconstrained form. 
 }
}
\description{
  This method function extracts the coefficients associated with the
  positive-definite matrix represented by \code{object}. 
}
\value{
  a vector with the coefficients corresponding to \code{object}.
}
\author{Jose Pinheiro and Douglas Bates }
\section{SIDE EFFECTS}{
  On the left side of an assignment, sets the values of the coefficients
  of \code{object} to \code{value}. \code{Object} must be initialized (using
  \code{initialize}) before new values can be assigned to its
  coefficients. 
}

\seealso{\code{\link{coef.pdMat}}, \code{\link{pdMat}}}

\examples{
coef(pdIdent(diag(3)))
}
\keyword{models}
