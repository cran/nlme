% $Id: formula.nlme.Rd,v 1.8 2001/10/30 20:51:15 bates Exp $
\name{formula.nlme}
\title{Extract nlme Object Formula}
\usage{
\method{formula}{nlme}(x, \dots)
}
\alias{formula.nlme}
\arguments{
 \item{x}{an object inheriting from class \code{nlme}, representing
   a fitted nonlinear mixed-effects model.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  This method function extracts the nonlinear model formula
  associated with \code{x}. 
}
\value{
  a two-sided nonlinear formula specifying the model used to
  obtain \code{x}.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

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
