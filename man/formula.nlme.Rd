% $Id: formula.nlme.Rd,v 1.3 1999/11/18 22:31:04 saikat Exp $
\name{formula.nlme}
\title{Extract nlme Object Formula}
\usage{
formula(object)
}
\alias{formula.nlme}
\arguments{
 \item{object}{an object inheriting from class \code{nlme}, representing
   a fitted nonlinear mixed-effects model.}
}
\description{
  This method function extracts the nonlinear model formula
  associated with \code{object}. 
}
\value{
  a two-sided nonlinear formula specifying the model used to
  obtain \code{object}.
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{nlme}}}

\examples{
library(nlme)
data(Soybean)
%\dontrun{
fm1 <- nlme(weight ~ SSlogis(Time, Asym, xmid, scal), data = Soybean,
            fixed = Asym + xmid + scal ~ 1, start = c(18, 52, 7.5),
            control = list(nlmStepMax = 1.0))
formula(fm1)
%}
}
\keyword{models}
