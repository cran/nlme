% $Id: fixed.effects.lmList.Rd,v 1.2 2000/04/17 19:03:39 bates Exp $
\name{fixed.effects.lmList}
\title{Extract lmList Fixed Effects}
\usage{
fixed.effects(object)
}
\alias{fixed.effects.lmList}
\alias{fixef.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
}
\description{
  The average of the coefficients corresponding to the \code{lm}
  components of \code{object} is calculated.
}
\value{
  a vector with the average of the individual \code{lm} coefficients in
  \code{object}. 
}

\author{Jose Pinheiro and Douglas Bates}

\seealso{\code{\link{lmList}}, \code{\link{random.effects.lmList}}}

\examples{
data(Orthodont)
fm1 <- lmList(distance ~ age | Subject, Orthodont)
fixed.effects(fm1)
}
\keyword{models}
