% $Id: formula.gnls.Rd,v 1.8 2002/03/05 14:59:39 bates Exp $
\name{formula.gnls}
\title{Extract gnls Object Formula}
\usage{
\method{formula}{gnls}(x, \dots)
}
\alias{formula.gnls}
\arguments{
 \item{x}{an object inheriting from class \code{gnls}, representing
   a generalized nonlinear least squares fitted model.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  This method function extracts the nonlinear model formula
  associated with \code{x}. 
}
\value{
  a two-sided formula specifying the nonlinear model used to
  obtain \code{x}.
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{gnls}}}

\examples{
data(Soybean)
fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
            weights = varPower())
formula(fm1)
}
\keyword{models}
