% $Id: residuals.gnls.Rd,v 1.4.2.1 2000/12/02 21:06:08 bates Exp $
\name{residuals.gnls}
\title{Extract gnls Residuals}
\usage{
residuals(object, type)
}
\alias{residuals.gnls}
\arguments{
 \item{object}{an object inheriting from class \code{gnls}, representing
   a generalized least squares fitted linear model.}
  \item{type}{an optional character string specifying the type of
    residuals to be used. If \code{"response"}, the "raw" residuals
    (observed - fitted) are used; else, if \code{"pearson"}, the
    standardized residuals (raw residuals divided by the corresponding
    standard errors) are used; else, if \code{"normalized"}, the
    normalized residuals (standardized residuals pre-multiplied by the
    inverse square-root factor of the estimated error correlation
    matrix) are used. Partial matching of arguments is used, so only the
    first character needs to be provided. Defaults to \code{"pearson"}.} 
}
\description{
  The residuals for the linear model represented by \code{object}
  are extracted.
}
\value{
  a vector with the residuals for the linear model represented by
  \code{object}. 
}
\author{Jose Pinheiro and Douglas Bates}
\seealso{\code{\link{gnls}}, \code{\link{fitted.gnls}} }
\examples{
data(Soybean)
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
fitted(fm1)
}
\keyword{models}
