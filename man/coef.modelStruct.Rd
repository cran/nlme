% $Id: coef.modelStruct.Rd,v 1.4 2000/07/03 18:22:46 bates Exp $
\name{coef.modelStruct}
\title{Extract modelStruct Object Coefficients}
\usage{
coef(object, unconstrained)
coef(object) <- value
}
\alias{coef.modelStruct}
\alias{coef<-.modelStruct}
\arguments{
 \item{object}{an object inheriting from class \code{modelStruct},
   representing a list of model components, such as \code{corStruct} and
   \code{varFunc} objects.}
 \item{unconstrained}{a logical value. If \code{TRUE} the coefficients
   are returned in unconstrained form (the same used in the optimization
   algorithm). If \code{FALSE} the coefficients are returned in
   "natural", possibly constrained, form. Defaults to \code{TRUE}.}
 \item{value}{a vector with the replacement values for the coefficients
   associated with \code{object}. It must be a vector with the same length
   of \code{coef{object}} and must be given in unconstrained form.}
}
\description{
  This method function extracts the coefficients associated with each
  component of the \code{modelStruct} list.
}
\value{
  a vector with all coefficients corresponding to the components of
  \code{object}.
}
\author{Jose Pinheiro and Douglas Bates}

\section{SIDE EFFECTS}{
  On the left side of an assignment, sets the values of the coefficients
  of \code{object} to \code{value}. \code{Object} must be initialized (using
  \code{initialize}) before new values can be assigned to its
  coefficients. 
}

\seealso{\code{\link{initialize}}}

\examples{
lms1 <- lmeStruct(reStruct = reStruct(pdDiag(diag(2), ~age)),
   corStruct = corAR1(0.3))
coef(lms1)
}
\keyword{models}
