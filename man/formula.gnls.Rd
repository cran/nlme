% $Id: formula.gnls.Rd,v 1.5 2001/01/10 19:21:04 bates Exp $
\name{formula.gnls}
\title{Extract gnls Object Formula}
\usage{
formula(object)
}
\alias{formula.gnls}
\arguments{
 \item{object}{an object inheriting from class \code{gnls}, representing
   a generalized nonlinear least squares fitted model.}
}
\description{
  This method function extracts the nonlinear model formula
  associated with \code{object}. 
}
\value{
  a two-sided formula specifying the nonlinear model used to
  obtain \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gnls}}}

\examples{
data(Soybean)
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
formula(fm1)
}
\keyword{models}
