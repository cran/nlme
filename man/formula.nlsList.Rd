% $Id: formula.nlsList.Rd,v 1.7 2001/10/30 20:51:15 bates Exp $
\name{formula.nlsList}
\title{Extract nlsList Object Formula}
\usage{
\method{formula}{nlsList}(x, \dots)
}
\alias{formula.nlsList}
\arguments{
 \item{x}{an object inheriting from class \code{nlsList}, representing
   a list of \code{nls} objects with a common model.
 }
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  This method function extracts the common nonlinear model formula
  associated with each \code{nls} component of \code{x}. 
}
\value{
  a two-sided nonlinear formula specifying the model used to
  obtain the \code{nls} components of \code{x}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{nlsList}}}

\examples{
data(Loblolly)
fm1 <- nlsList(SSasymp, Loblolly)
formula(fm1)
}
\keyword{models}
