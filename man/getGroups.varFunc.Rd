% $Id: getGroups.varFunc.Rd,v 1.4 2000/07/03 18:22:47 bates Exp $
\name{getGroups.varFunc}
\title{Extract varFunc Groups}
\usage{
getGroups(object)
}
\alias{getGroups.varFunc}
\arguments{
 \item{object}{an object inheriting from class \code{varFunc},
   representing a variance function structure.} 
}
\description{
  This method function extracts the grouping factor associated with the
  variance function represented by \code{object}, if any is present.
}
\value{
  if \code{object} has a \code{groups} attribute, its value is
  returned; else \code{NULL} is returned.
}

\author{Jose Pinheiro and Douglas Bates}

\examples{
data(Orthodont)
vf1 <- varPower(form = ~ age | Sex)
vf1 <- initialize(vf1, Orthodont)
getGroups(vf1)
}
\keyword{models}
