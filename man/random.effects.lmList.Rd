% $Id: random.effects.lmList.Rd,v 1.4 2000/07/03 18:22:48 bates Exp $
\name{random.effects.lmList}
\title{Extract lmList Random Effects}
\usage{
random.effects(object)
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
