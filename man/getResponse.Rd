% $Id: getResponse.Rd,v 1.1 1999/10/13 00:50:11 saikat Exp $
\name{getResponse}
\title{Extract Response Variable from an Object}
\usage{
getResponse(object, form, data)
}
\alias{getResponse}
\arguments{
 \item{object}{any object}
 \item{form}{an optional two-sided formula. Defaults to
   \code{formula(object)}.}
 \item{data}{a data frame in which to interpret the variables named in
   \code{form}. Optional for most methods.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include \code{data.frame}, \code{gls}, \code{lme},
  and \code{lmList}.
}
\value{
  will depend on the method function used; see the appropriate documentation.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{getResponseFormula}}}
\examples{
## see the method function documentation
}
\keyword{models}
