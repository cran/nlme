% $Id: formula.nlsList.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
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
data(Soybean)
%\dontrun{
fm1 <- nlsList(weight ~ SSlogis(Time, Asym, xmid, scal)|Plot, data=Soybean)
formula(fm1)
%}
}
\keyword{models}
