% $Id: coef.pdIdent.Rd,v 1.6.2.1 2002/08/09 19:45:29 bates Exp $
\name{coef.pdIdent}
\title{pdIdent Object Coefficients}
\usage{
\method{coef}{pdIdent}(object, unconstrained, \dots)
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
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
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
  \code{Initialize}) before new values can be assigned to its
  coefficients. 
}

\seealso{\code{\link{coef.pdMat}}, \code{\link{pdMat}}}

\examples{
coef(pdIdent(diag(3)))
}
\keyword{models}
