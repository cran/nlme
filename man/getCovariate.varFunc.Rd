% $Id: getCovariate.varFunc.Rd,v 1.7 2002/03/05 14:59:39 bates Exp $
\name{getCovariate.varFunc}
\title{Extract varFunc Covariate}
\usage{
\method{getCovariate}{varFunc}(object, form, data)
}
\alias{getCovariate.varFunc}
\arguments{
 \item{object}{an object inheriting from class \code{varFunc},
   representing a variance function structure.}
 \item{form}{an optional formula specifying the covariate to be
   evaluated in \code{object}. Defaults to \code{formula(object)}.}
 \item{data}{some methods for this generic require a \code{data}
   object.  Not used in this method.} 
}
\description{
  This method function extracts the covariate(s) associated with the
  variance function represented by \code{object}, if any is present.
}
\value{
  if \code{object} has a \code{covariate} attribute, its value is
  returned; else \code{NULL} is returned.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{covariate<-.varFunc}}}

\examples{
data(Orthodont)
vf1 <- varPower(1.1, form = ~age)
covariate(vf1) <- Orthodont[["age"]]
getCovariate(vf1)
}
\keyword{models}
