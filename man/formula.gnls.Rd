% $Id: formula.gnls.Rd,v 1.4.2.1 2000/12/02 21:06:08 bates Exp $
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
