% $Id: getResponse.Rd,v 1.4 2001/06/18 21:16:50 bates Exp $
\name{getResponse}
\title{Extract Response Variable from an Object}
\usage{
getResponse(object, form)
}
\alias{getResponse}
\arguments{
 \item{object}{any object}
 \item{form}{an optional two-sided formula. Defaults to
   \code{formula(object)}.}
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
