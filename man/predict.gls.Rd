% $Id: predict.gls.Rd,v 1.5 2001/01/10 19:21:05 bates Exp $
\name{predict.gls}
\title{Predictions from a gls Object}
\usage{
predict(object, newdata, na.action)
}
\alias{predict.gls}
\arguments{
 \item{object}{an object inheriting from class \code{gls}, representing
   a generalized least squares fitted linear model.}
 \item{newdata}{an optional data frame to be used for obtaining the
   predictions. All variables used in the linear model must be present
   in the data frame. If missing, the fitted values are returned.}
 \item{na.action}{a function that indicates what should happen when
   \code{newdata} contains \code{NA}s.  The default action
   (\code{na.fail}) causes the function to print an error message and
   terminate if there are any incomplete observations.}
}
\description{
  The predictions for the linear model represented by \code{object} are
  obtained at the covariate values defined in \code{newdata}.
}
\value{
  a vector with the predicted values.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{gls}}, \code{\link{fitted.gls}}}

\examples{
data(Ovary)
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
           correlation = corAR1(form = ~ 1 | Mare))
newOvary <- data.frame(Time = c(-0.75, -0.5, 0, 0.5, 0.75))
predict(fm1, newOvary)
}
\keyword{models}
