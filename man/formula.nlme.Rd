% $Id: formula.nlme.Rd,v 1.6 2001/01/10 19:21:04 bates Exp $
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
data(Loblolly)
fm1 <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
            data = Loblolly,
            fixed = Asym + R0 + lrc ~ 1,
            random = Asym ~ 1,
            start = c(Asym = 103, R0 = -8.5, lrc = -3.3))
formula(fm1)
}
\keyword{models}
