% $Id: coef.gnls.Rd,v 1.6 2001/06/18 21:16:50 bates Exp $
\name{coef.gnls}
\title{Extract gnls Coefficients}
\usage{
\method{coef}{gnls}(object)
}
\alias{coef.gnls}
\arguments{
 \item{object}{an object inheriting from class \code{gnls}, representing
   a generalized nonlinear least squares fitted model.}
}
\description{
  The estimated coefficients for the nonlinear model represented by
  \code{object} are extracted. 
}
\value{
  a vector with the estimated coefficients for the nonlinear model
  represented by \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gnls}}} 

\examples{
data(Soybean)
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
coef(fm1)
}
\keyword{models}
