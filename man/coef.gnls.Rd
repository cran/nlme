% $Id: coef.gnls.Rd,v 1.5 2000/07/03 18:22:46 bates Exp $
\name{coef.gnls}
\title{Extract gnls Coefficients}
\usage{
coef(object)
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
