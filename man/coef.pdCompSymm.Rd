% $Id: coef.pdCompSymm.Rd,v 1.6 2001/10/30 20:51:15 bates Exp $
\name{coef.pdCompSymm}
\title{pdCompSymm Object Coefficients}
\usage{
\method{coef}{pdCompSymm}(object, unconstrained, \dots)
coef(object) <- value
}
\alias{coef.pdCompSymm}
\arguments{
 \item{object}{an object inheriting from class \code{pdCompSymm},
   representing a positive-definite matrix with compound symmetry
   structure.}
 \item{unconstrained}{a logical value. If \code{TRUE} the coefficients
   are returned in unconstrained form (the same used in the optimization
   algorithm). If \code{FALSE} the standard deviation and the
   correlation coefficient of the compound symmetry  of
   positive-definite matrix represented by \code{object} are
   returned. Defaults to \code{TRUE}.}
 \item{value}{a vector with the replacement values for the coefficients
   associated with \code{object}. It must be a vector of length two
   and must be given in unconstrained form.
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
  \code{initialize}) before new values can be assigned to its
  coefficients. 
}
\seealso{\code{\link{coef.pdMat}}, \code{\link{pdMat}}}

\examples{
coef(pdCompSymm(diag(3)), FALSE)
}
\keyword{models}
