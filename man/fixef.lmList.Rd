% $Id: fixef.lmList.Rd,v 1.3 2001/10/30 20:51:15 bates Exp $
\name{fixef.lmList}
\title{Extract lmList Fixed Effects}
\usage{
\method{fixef}{lmList}(object, \dots)
}
\alias{fixed.effects.lmList}
\alias{fixef.lmList}
\arguments{
 \item{object}{an object inheriting from class \code{lmList}, representing
   a list of \code{lm} objects with a common model.
 }
 \item{\dots}{some methods for this generic require additional
   arguments.  None are used in this method.} 
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
