% $Id: random.effects.lmList.Rd,v 1.5 2001/06/18 21:16:50 bates Exp $
\name{random.effects.lmList}
\title{Extract lmList Random Effects}
\usage{
\method{random.effects}{lmList}(object)
}
\alias{random.effects.lmList}
\alias{ranef.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
}
\description{
  The difference between the individual \code{lm} components
  coefficients and their average is calculated.
}
\value{
  a vector with the differences between the individual \code{lm}
  coefficients in \code{object} and their average. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}, \code{\link{fixed.effects.lmList}}}
 
\examples{
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
random.effects(fm1)
}
\keyword{models}
