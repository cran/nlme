% $Id: getCovariate.varFunc.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{getCovariate.varFunc}
\title{Extract varFunc Covariate}
\usage{
getCovariate(object)
}
\alias{getCovariate.varFunc}
\arguments{
 \item{object}{an object inheriting from class \code{varFunc},
   representing a variance function structure.} 
}
\description{
  This method function extracts the covariate(s) associated with the
  variance function represented by \code{object}, if any is present.
}
\value{
  if \code{object} has a \code{covariate} attribute, its value is
  returned; else \code{NULL} is returned.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{covariate<-.varFunc}}}

\examples{
data(Orthodont)
vf1 <- varPower(1.1, form = ~age)
covariate(vf1) <- Orthodont[["age"]]
getCovariate(vf1)
}
\keyword{models}
