% $Id: fitted.gnls.Rd,v 1.4.2.1 2000/12/02 21:06:08 bates Exp $
\name{fitted.gnls}
\title{Extract gnls Fitted Values}
\usage{
fitted(object)
}
\alias{fitted.gnls}
\arguments{
 \item{object}{an object inheriting from class \code{gnls}, representing
   a generalized nonlienar least squares fitted model.}
}
\description{
  The fitted values for the nonlinear model represented by \code{object}
  are extracted.
}
\value{
  a vector with the fitted values for the nonlinear model represented by
  \code{object}. 
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{gnls}}, \code{\link{residuals.gnls}} }
\examples{
data(Soybean)
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
fitted(fm1)
}
\keyword{models}
