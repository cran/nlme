% $Id: predict.gnls.Rd,v 1.8 2002/03/05 14:59:39 bates Exp $
\name{predict.gnls}
\title{Predictions from a gnls Object}
\usage{
\method{predict}{gnls}(object, newdata, na.action, naPattern, \dots)
}
\alias{predict.gnls}
\arguments{
 \item{object}{an object inheriting from class \code{gnls}, representing
   a generalized nonlinear least squares fitted model.}
 \item{newdata}{an optional data frame to be used for obtaining the
   predictions. All variables used in the nonlinear model must be present
   in the data frame. If missing, the fitted values are returned.}
 \item{na.action}{a function that indicates what should happen when
   \code{newdata} contains \code{NA}s.  The default action
   (\code{na.fail}) causes the function to print an error message and
   terminate if there are any incomplete observations.}
 \item{naPattern}{an expression or formula object, specifying which returned
   values are to be regarded as missing.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  The predictions for the nonlinear model represented by \code{object} are
  obtained at the covariate values defined in \code{newdata}.
}
\value{
  a vector with the predicted values.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{gnls}}, \code{\link{fitted.gnls}}}

\examples{
data(Soybean)
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
newSoybean <- data.frame(Time = c(10,30,50,80,100))
predict(fm1, newSoybean)
}
\keyword{models}
