% $Id: getData.Rd,v 1.1 1999/10/13 00:50:11 saikat Exp $
\name{getData}
\title{Extract Data from an Object}
\usage{
getData(object)
}
\alias{getData}
\arguments{
 \item{object}{an object from which a data.frame can be extracted,
   generally a fitted model object.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include \code{gls}, \code{lme},
  and \code{lmList}.
}
\value{
  will depend on the method function used; see the appropriate documentation.
}

\author{Jose Pinheiro and Douglas Bates}

\examples{
## see the method function documentation
}
\keyword{models}
