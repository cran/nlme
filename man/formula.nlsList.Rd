% $Id: formula.nlsList.Rd,v 1.4.2.2 2000/12/07 21:54:10 bates Exp $
\name{formula.nlsList}
\title{Extract nlsList Object Formula}
\usage{
formula(object)
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
