% $Id: initialize.Rd,v 1.1 1999/10/13 00:50:11 saikat Exp $
\name{initialize}
\title{Initialize Object}
\usage{
initialize(object, data, ...)
}
\alias{initialize}
\arguments{
 \item{object}{any object requiring initialization, e.g. "plug-in"
   structures such as \code{corStruct} and \code{varFunc} objects. }
 \item{data}{a data frame to be used in the initialization procedure.}
 \item{\dots}{some methods for this generic function require additional
   arguments.}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. Classes which already have methods for
  this function include: \code{corStruct}, \code{lmeStruct},
  \code{reStruct}, and \code{varFunc}.
}
\value{
  an initialized object with the same class as \code{object}. Changes
  introduced by the initialization procedure will depend on the method
  function used; see the appropriate documentation. 
}
\author{Jose Pinheiro and Douglas Bates}
\examples{
## see the method function documentation
}
\keyword{models}
