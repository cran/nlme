% $Id: fixef.lme.Rd,v 1.4 2002/03/05 14:59:39 bates Exp $
\name{fixef.lme}
\title{Extract lme Fixed Effects}
\usage{
\method{fixef}{lme}(object, \dots)
}
\alias{fixed.effects.lme}
\alias{fixef.lme}
\arguments{
 \item{object}{an object inheriting from class \code{lme}, representing
   a fitted linear mixed-effects model.}
  \item{\dots}{some methods for this generic require additional
   arguments.  None are used in this method.} 
}
\description{
  The fixed effects estimates corresponding to the linear mixed-effects
  model represented by \code{object} are returned.
}
\value{
  a vector with the fixed effects estimates corresponding to
  \code{object}. 
}

\author{Jose Pinheiro \email{Jose.Pinheiro@pharma.novartis.com} and Douglas Bates \email{bates@stat.wisc.edu}}

\seealso{\code{\link{coef.lme}}, \code{\link{random.effects.lme}} }

\examples{
data(Orthodont)
fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
fixed.effects(fm1)
fixef(fm1)
}
\keyword{models}
