% $Id: needUpdate.Rd,v 1.5 2002/03/05 14:59:39 bates Exp $
\name{needUpdate}
\title{Check if Update is Needed}
\usage{
needUpdate(object)
}
\alias{needUpdate}
\alias{needUpdate.default}
\alias{needUpdate.varComb}
\alias{needUpdate.varIdent}
\arguments{
 \item{object}{any object}
}
\description{
  This function is generic; method functions can be written to handle
  specific classes of objects. By default, it tries to extract a
  \code{needUpdate} attribute of \code{object}. If this is \code{NULL}
  or \code{FALSE} it returns \code{FALSE}; else it returns \code{TRUE}.
  Updating of objects usually takes place in iterative algorithms in
  which auxiliary quantities associated with the object, and not being
  optimized over, may change.
}
\value{
  a logical value indicating whether \code{object} needs to be updated. 
}
\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\examples{
data(Orthodont)
vf1 <- varExp()
vf1 <- initialize(vf1, data = Orthodont)
needUpdate(vf1)
}
\keyword{models}
