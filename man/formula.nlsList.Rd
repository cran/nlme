% $Id: formula.nlsList.Rd,v 1.6 2001/06/18 21:16:50 bates Exp $
\name{formula.nlsList}
\title{Extract nlsList Object Formula}
\usage{
\method{formula}{nlsList}(object)
}
\alias{formula.nlsList}
\arguments{
 \item{object}{an object inheriting from class \code{nlsList}, representing
   a list of \code{nls} objects with a common model.
 }
}
\description{
  This method function extracts the common nonlinear model formula
  associated with each \code{nls} component of \code{object}. 
}
\value{
  a two-sided nonlinear formula specifying the model used to
  obtain the \code{nls} components of \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{nlsList}}}

\examples{
data(Loblolly)
fm1 <- nlsList(SSasymp, Loblolly)
formula(fm1)
}
\keyword{models}
